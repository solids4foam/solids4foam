/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fluidSolidInterface.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"
#include "tetPointFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetPolyPatchInterpolation.H"
#include "tetFemMatrices.H"
#include "fixedValuePointPatchFields.H"
#include "ZoneIDs.H"
#include "TPSFunction.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"
#include "newSubsetMotionSolverFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSolidInterface, 0);
    defineRunTimeSelectionTable(fluidSolidInterface, dictionary);

    template<>
    const char*
    NamedEnum<fluidSolidInterface::interfaceTransferMethod, 3>::names[] =
    {
        "directMap",
        "RBF",
        "GGI"
    };


    const NamedEnum<fluidSolidInterface::interfaceTransferMethod, 3>
        fluidSolidInterface::interfaceTransferMethodNames_;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fluidSolidInterface::updateCoupled()
{
    if (couplingStartTime_ > SMALL && !coupled_)
    {
        if (runTime().value() > (couplingStartTime_ - SMALL))
        {
            InfoIn("fluidSolidInterface::updateCoupled()")
                << "Enabling fluid-solid coupling" << endl;

            // Enable coupling
            coupled_ = true;

            return true;
        }
    }

    return false;
}


void Foam::fluidSolidInterface::calcCurrentSolidZonesPoints() const
{
    const List<tmp<vectorField> > currentFaceZonesPoints
    (
        solid().currentFaceZonesPoints()
    );

    currentSolidZonesPointsPtrList_.setSize(nGlobalPatches());

    forAll(solid().globalPatches(), i)
    {
        // Find global face zones
        if (currentSolidZonesPointsPtrList_.set(i))
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::"
                "calcCurrentSolidZonesPoints() const"
            )   << "Current solid zone points already exists "
                << "for global patch: "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[i].patch().index()
                   ].name()
                << "!" << abort(FatalError);
        }

        currentSolidZonesPointsPtrList_.set
        (
            i,
            new vectorField(currentFaceZonesPoints[i]())
        );
    }
}


void Foam::fluidSolidInterface::calcCurrentSolidZonesPatches() const
{
    currentSolidZonesPatchesPtrList_.setSize(nGlobalPatches());

    forAll(solid().globalPatches(), i)
    {
        // Find global face zones
        if (currentSolidZonesPatchesPtrList_.set(i))
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::"
                "calcCurrentSolidZonesPatches() const"
            )   << "Current solid zone patch already exists "
                << "for global patch: "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[i].patch().index()
                   ].name()
                << "!" << abort(FatalError);
        }

        currentSolidZonesPatchesPtrList_.set
        (
            i,
            new standAlonePatch
            (
                solid().globalPatches()[i].globalPatch().localFaces(),
                currentSolidZonesPoints()[i]
            )
        );
    }
}


void Foam::fluidSolidInterface::calcRbfFluidToSolidInterpolators() const
{
    rbfFluidToSolidPtrList_.setSize(nGlobalPatches());

    forAll(fluid().globalPatches(), i)
    {
        const label patchID = fluid().globalPatches()[i].patch().index();

        if (rbfFluidToSolidPtrList_[i])
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::calcRbfFluidToSolidInterpolators() const"
            )   << "Fluid to solid interpolator already exists for global patch: "
                << fluidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        Info<< "Create RBF fluid-to-solid interpolator for interface patch: "
            << fluidMesh().boundary()[patchID].name() << "!" << endl;

        std::shared_ptr<RBFFunctionInterface> rbfFunction;
        rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

        rbfFluidToSolidPtrList_[i] =
            std::shared_ptr<RBFInterpolation>(new RBFInterpolation(rbfFunction));

        const vectorField solidZoneFaceCentres =
            currentSolidZonesPatches()[i].faceCentres();

        const vectorField fluidZoneFaceCentres =
            fluid().globalPatches()[i].globalPatch().faceCentres();

        matrix fluidX(fluidZoneFaceCentres.size(), 3);
        matrix solidX(solidZoneFaceCentres.size(), 3);

        forAll(fluidZoneFaceCentres, faceI)
        {
            fluidX(faceI, 0) = fluidZoneFaceCentres[faceI].x();
            fluidX(faceI, 1) = fluidZoneFaceCentres[faceI].y();
            fluidX(faceI, 2) = fluidZoneFaceCentres[faceI].z();
        }

        forAll(solidZoneFaceCentres, faceI)
        {
            solidX(faceI, 0) = solidZoneFaceCentres[faceI].x();
            solidX(faceI, 1) = solidZoneFaceCentres[faceI].y();
            solidX(faceI, 2) = solidZoneFaceCentres[faceI].z();
        }

        rbfFluidToSolidPtrList_[i]->compute(fluidX, solidX);

        Info<< "Checking fluid-to-solid interpolator" << endl;

        {
            vectorField fluidPatchFaceCentres =
                vectorField
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndices()[i]].faceCentres()
                );

            vectorField fluidZoneFaceCentres =
                fluid().globalPatches()[i].patchFaceToGlobal(fluidPatchFaceCentres);


            vectorField solidPatchFaceCentres =
                vectorField
                (
                    solidMesh().boundaryMesh()[solidPatchIndices()[i]].faceCentres()
                );

            matrix fluidX(fluidZoneFaceCentres.size(), 3);
            matrix fluidXsolid(solidPatchFaceCentres.size(), 3);

            forAll(fluidZoneFaceCentres, faceI)
            {
                fluidX(faceI, 0) = fluidZoneFaceCentres[faceI].x();
                fluidX(faceI, 1) = fluidZoneFaceCentres[faceI].y();
                fluidX(faceI, 2) = fluidZoneFaceCentres[faceI].z();
            }

            rbfFluidToSolidPtrList_[i]->interpolate(fluidX, fluidXsolid);

            vectorField fluidPatchFaceCentresAtSolid
            (
                solidPatchFaceCentres.size(),
                vector::zero
            );

            forAll(fluidPatchFaceCentresAtSolid, faceI)
            {
                fluidPatchFaceCentresAtSolid[faceI].x() = fluidXsolid(faceI, 0);
                fluidPatchFaceCentresAtSolid[faceI].y() = fluidXsolid(faceI, 1);
                fluidPatchFaceCentresAtSolid[faceI].z() = fluidXsolid(faceI, 2);
            }

            scalar maxDist = gMax
            (
                mag
                (
                    fluidPatchFaceCentresAtSolid
                  - solidPatchFaceCentres
                )
            );

            Info<< "Fluid-to-solid face interpolation error: " << maxDist
                << endl;
        }

        Info<< endl;
    }
}


void Foam::fluidSolidInterface::calcRbfSolidToFluidInterpolators() const
{
    rbfSolidToFluidPtrList_.setSize(nGlobalPatches());

    forAll(solid().globalPatches(), i)
    {
        const label patchID = solid().globalPatches()[i].patch().index();

        if (rbfSolidToFluidPtrList_[i])
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::calcRbfSolidToFluidInterpolators() const"
            )   << "Solid to fluid interpolator already exists for global patch: "
                << solidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        Info<< "Create RBF solid-to-fluid interpolator for interface patch: "
            << solidMesh().boundary()[patchID].name() << endl;

        std::shared_ptr<RBFFunctionInterface> rbfFunction;
        rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

        rbfSolidToFluidPtrList_[i] =
            std::shared_ptr<RBFInterpolation>(new RBFInterpolation(rbfFunction));

        const vectorField solidZonePoints =
            currentSolidZonesPatches()[i].localPoints();

        const vectorField fluidZonePoints =
            fluid().globalPatches()[i].globalPatch().localPoints();

        matrix fluidX(fluidZonePoints.size(), 3);
        matrix solidX(solidZonePoints.size(), 3);

        forAll(fluidZonePoints, faceI)
        {
            fluidX(faceI, 0) = fluidZonePoints[faceI].x();
            fluidX(faceI, 1) = fluidZonePoints[faceI].y();
            fluidX(faceI, 2) = fluidZonePoints[faceI].z();
        }

        forAll(solidZonePoints, faceI)
        {
            solidX(faceI, 0) = solidZonePoints[faceI].x();
            solidX(faceI, 1) = solidZonePoints[faceI].y();
            solidX(faceI, 2) = solidZonePoints[faceI].z();
        }

       rbfSolidToFluidPtrList_[i]->compute(solidX, fluidX);

       Info<< "Checking solid-to-fluid interpolator" << endl;

       {
           matrix fluidPoints(fluidZonePoints.size(), 3);
           matrix solidPoints(solidZonePoints.size(), 3);
           vectorField fluidZonePointsInterp(fluidZonePoints.size(), vector::zero);


           forAll(solidZonePoints, faceI)
           {
               solidPoints(faceI, 0) = solidZonePoints[faceI].x();
               solidPoints(faceI, 1) = solidZonePoints[faceI].y();
               solidPoints(faceI, 2) = solidZonePoints[faceI].z();
           }

           rbfSolidToFluidPtrList_[i]->interpolate(solidPoints, fluidPoints);

           forAll(fluidZonePoints, faceI)
           {
               fluidZonePointsInterp[faceI].x() = fluidPoints(faceI, 0);
               fluidZonePointsInterp[faceI].y() = fluidPoints(faceI, 1);
               fluidZonePointsInterp[faceI].z() = fluidPoints(faceI, 2);
           }

           scalar maxDist = gMax
           (
               mag
               (
                   fluidZonePointsInterp
                 - fluidZonePoints
               )
           );

           Info<< "Solid-to-fluid point interpolation error: " << maxDist
               << endl;
        }

        Info<< endl;
    }
}


void Foam::fluidSolidInterface::calcGgiInterpolators() const
{
    // Note: Remove current solid zone and points so that
    // it will be re-created in the deformed position
    currentSolidZonesPatchesPtrList_.clear();
    currentSolidZonesPointsPtrList_.clear();

    ggiInterpolatorsPtrList_.setSize(nGlobalPatches());

    forAll(fluid().globalPatches(), i)
    {
        const label patchID = fluid().globalPatches()[i].patch().index();

        // Create ggi interpolation
        if (ggiInterpolatorsPtrList_.set(i))
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::"
                "calcGgiInterpolators() const"
            )   << "Ggi interpolator already exists for global patch: "
                << fluidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        Info<< "Create GGI zone-to-zone interpolator for interface patch: "
            << fluidMesh().boundary()[patchID].name() << endl;

        ggiInterpolatorsPtrList_.set
        (
            i,
            new GGIInterpolation<standAlonePatch, standAlonePatch>
            (
                fluid().globalPatches()[i].globalPatch(),
                currentSolidZonesPatches()[i],
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation. Bug fix
                true,           // Patch data is complete on all processors
                SMALL,          // Non-overlapping face tolerances
                SMALL,
                true,           // Rescale weighting factors
                ggiInterpolation::BB_OCTREE
            )
        );

        Info<< "Checking fluid-to-solid face interpolator" << endl;

        {
            const vectorField fluidPatchFaceCentres =
                vectorField
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndices()[i]].faceCentres()
                );

            const vectorField fluidZoneFaceCentres =
                fluid().globalPatches()[i].patchFaceToGlobal(fluidPatchFaceCentres);

            const vectorField solidZoneFaceCentres =
                ggiInterpolatorsPtrList_(i)->masterToSlave
                (
                    fluidZoneFaceCentres
                );

            const vectorField solidPatchFaceCentres =
                solid().globalPatches()[i].globalFaceToPatch(solidZoneFaceCentres);

            scalar maxDist = gMax
            (
                mag
                (
                    solidPatchFaceCentres
                  - solidMesh().boundaryMesh()[solidPatchIndices()[i]].faceCentres()
                )
            );

            Info<< "Fluid-to-solid face interpolation error: " << maxDist
                << endl;
        }

        Info<< "Checking solid-to-fluid point interpolator" << endl;

        {
            const vectorField solidZonePointsAtFluid =
                ggiInterpolatorsPtrList_(i)->slaveToMasterPointInterpolate
                (
                    currentSolidZonesPoints()[i]
                );

            const vectorField fluidZonePoints =
                fluid().globalPatches()[i].globalPatch().localPoints();

            const scalar maxDist = gMax
            (
                mag
                (
                    fluidZonePoints
                  - solidZonePointsAtFluid
                )
            );

            Info<< "Solid-to-fluid point interpolation error: " << maxDist
                << endl;
        }

        Info<< "Number of uncovered master faces: "
            << ggiInterpolatorsPtrList_(i)->uncoveredMasterFaces().size() << endl;

        Info<< "Number of uncovered slave faces: "
            << ggiInterpolatorsPtrList_(i)->uncoveredSlaveFaces().size() << endl;

        ggiInterpolatorsPtrList_(i)->slavePointDistanceToIntersection();
        ggiInterpolatorsPtrList_(i)->masterPointDistanceToIntersection();

        Info<< endl;
    }
}


void Foam::fluidSolidInterface::calcFluidToSolidFaceMaps() const
{
    fluidToSolidFaceMapsPtrList_.setSize(nGlobalPatches());

    forAll(fluid().globalPatches(), i)
    {
        const label patchID = fluid().globalPatches()[i].patch().index();

        if (fluidToSolidFaceMapsPtrList_.set(i))
        {
            FatalErrorIn(type() + "::calcFluidToSolidFaceMaps() const")
                << "pointer already set for global patch: "
                << fluidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        if
        (
            solid().globalPatches()[i].globalPatch().size()
         != fluid().globalPatches()[i].globalPatch().size()
        )
        {
            FatalErrorIn(type() + "::calcFluidToSolidFaceMaps() const")
                << "Fluid and solid interfaces are not conformal for "
                << "global patch: "
                << fluidMesh().boundary()[patchID].name() << "!" << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word fluidToSolidFaceMapName
        (
            "fluidToSolidFaceMap"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
        );

        IOobject mapHeader
        (
            fluidToSolidFaceMapName,
            fluid().runTime().timeName(),
            fluidMesh(),
            IOobject::MUST_READ
        );

        if (mapHeader.headerOk())
        {
            // Read map
            Info<< "Reading fluidToSolidFaceMap for global patch "
                << fluidMesh().boundary()[patchID].name()
                << " from disk..." << endl;

            fluidToSolidFaceMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        fluidToSolidFaceMapName,
                        fluid().runTime().timeName(),
                        fluidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            // Initialise map
            fluidToSolidFaceMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        fluidToSolidFaceMapName,
                        fluid().runTime().timeName(),
                        fluidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    labelList(solid().globalPatches()[i].globalPatch().size(), -1)
                )
            );
            labelList& fluidToSolidMap = fluidToSolidFaceMapsPtrList_[i];

            // Perform N^2 search for corresponding faces
            // We will take 0.1% of the minEdgeLength as the exact match tolerance

            const vectorField& fluidCf =
                fluid().globalPatches()[i].globalPatch().faceCentres();

            const vectorField& solidCf =
                solid().globalPatches()[i].globalPatch().faceCentres();

            const scalar tol = 0.001*gMin(minEdgeLengths()[i]);

            forAll(solidCf, solidFaceI)
            {
                const vector& curSolidCf = solidCf[solidFaceI];

                forAll(fluidCf, fluidFaceI)
                {
                    if (mag(curSolidCf - fluidCf[fluidFaceI]) < tol)
                    {
                        fluidToSolidMap[solidFaceI] = fluidFaceI;
                        break;
                    }
                }
            }
        }

        if (gMin(fluidToSolidFaceMapsPtrList_[i]) == -1)
        {
            FatalErrorIn(type() + "::calcFluidToSolidFaceMaps() const")
                << "Cannot calculate the map between interfaces for global patch "
                << fluidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::fluidToSolidFaceMaps() const
{
    if (fluidToSolidFaceMapsPtrList_.empty())
    {
        calcFluidToSolidFaceMaps();
    }

    return fluidToSolidFaceMapsPtrList_;
}


void Foam::fluidSolidInterface::calcSolidToFluidFaceMaps() const
{
    solidToFluidFaceMapsPtrList_.setSize(nGlobalPatches());

    forAll(solid().globalPatches(), i)
    {
        const label patchID = solid().globalPatches()[i].patch().index();

        if (solidToFluidFaceMapsPtrList_.set(i))
        {
            FatalErrorIn(type() + "::calcSolidToFluidFaceMaps() const")
                << "pointer already set for global patch: "
                << solidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        if
        (
            solid().globalPatches()[i].globalPatch().size()
         != fluid().globalPatches()[i].globalPatch().size()
        )
        {
            FatalErrorIn(type() + "::calcSolidToFluidFaceMaps() const")
                << "Solid and fluid interfaces are not conformal for "
                << "global patch: "
                << solidMesh().boundary()[patchID].name() << "!" << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word solidToFluidFaceMapName
        (
            "solidToFluidFaceMap"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
        );

        IOobject mapHeader
        (
            solidToFluidFaceMapName,
            solid().runTime().timeName(),
            solidMesh(),
            IOobject::MUST_READ
        );

        if (mapHeader.headerOk())
        {
            // Read map
            Info<< "Reading solidToFluidFaceMap for global patch "
                << solidMesh().boundary()[patchID].name()
                << " from disk..." << endl;

            solidToFluidFaceMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        solidToFluidFaceMapName,
                        solid().runTime().timeName(),
                        solidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            // Initialise map
            solidToFluidFaceMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        solidToFluidFaceMapName,
                        solid().runTime().timeName(),
                        solidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    labelList(fluid().globalPatches()[i].globalPatch().size(), -1)
                )
            );
            labelList& solidToFluidMap = solidToFluidFaceMapsPtrList_[i];

            // Perform N^2 search for corresponding faces
            // We will take 0.1% of the minEdgeLength as the exact match tolerance

            const vectorField& fluidCf =
                fluid().globalPatches()[i].globalPatch().faceCentres();

            const vectorField& solidCf =
                solid().globalPatches()[i].globalPatch().faceCentres();

            const scalar tol = 0.001*gMin(minEdgeLengths()[i]);

            forAll(fluidCf, fluidFaceI)
            {
                const vector& curFluidCf = fluidCf[fluidFaceI];

                forAll(solidCf, solidFaceI)
                {
                    if (mag(curFluidCf - solidCf[solidFaceI]) < tol)
                    {
                        solidToFluidMap[fluidFaceI] = solidFaceI;
                        break;
                    }
                }
            }
        }

        if (gMin(solidToFluidFaceMapsPtrList_[i]) == -1)
        {
            FatalErrorIn(type() + "::calcSolidToFluidFaceMaps() const")
                << "Cannot calculate the map between interfaces for global patch "
                << solidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::solidToFluidFaceMaps() const
{
    if (solidToFluidFaceMapsPtrList_.empty())
    {
        calcSolidToFluidFaceMaps();
    }

    return solidToFluidFaceMapsPtrList_;
}


void Foam::fluidSolidInterface::calcFluidToSolidPointMaps() const
{
    fluidToSolidPointMapsPtrList_.setSize(nGlobalPatches());

    forAll(fluid().globalPatches(), i)
    {
        const label patchID = fluid().globalPatches()[i].patch().index();

        if (fluidToSolidPointMapsPtrList_.set(i))
        {
            FatalErrorIn(type() + "::calcFluidToSolidPointMaps() const")
                << "pointer already set for global patch: "
                << fluidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        if
        (
            solid().globalPatches()[i].globalPatch().size()
         != fluid().globalPatches()[i].globalPatch().size()
        )
        {
            FatalErrorIn(type() + "::calcFluidToSolidPointMaps() const")
                << "Fluid and solid interfaces are not conformal for "
                << "global patch: "
                << fluidMesh().boundary()[patchID].name() << "!" << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word fluidToSolidPointMapName
        (
            "fluidToSolidPointMap"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
        );

        IOobject mapHeader
        (
            fluidToSolidPointMapName,
            fluid().runTime().timeName(),
            fluidMesh(),
            IOobject::MUST_READ
        );

        if (mapHeader.headerOk())
        {
            // Read map
            Info<< "Reading fluidToSolidPointMap for global patch "
                << fluidMesh().boundary()[patchID].name()
                << " from disk..." << endl;

            fluidToSolidPointMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        fluidToSolidPointMapName,
                        fluid().runTime().timeName(),
                        fluidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            // Initialise map
            fluidToSolidPointMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        fluidToSolidPointMapName,
                        fluid().runTime().timeName(),
                        fluidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    labelList(solid().globalPatches()[i].globalPatch().nPoints(), -1)
                )
            );
            labelList& fluidToSolidMap = fluidToSolidPointMapsPtrList_[i];

            // Perform N^2 search for corresponding points
            // We will take 0.1% of the minEdgeLength as the exact match tolerance

            const vectorField& fluidLP =
                fluid().globalPatches()[i].globalPatch().localPoints();

            const vectorField& solidLP =
                solid().globalPatches()[i].globalPatch().localPoints();

            const scalar tol = 0.001*gMin(minEdgeLengths()[i]);

            forAll(solidLP, solidPointI)
            {
                const vector& curSolidLP = solidLP[solidPointI];

                forAll(fluidLP, fluidPointI)
                {
                    if (mag(curSolidLP - fluidLP[fluidPointI]) < tol)
                    {
                        fluidToSolidMap[solidPointI] = fluidPointI;
                        break;
                    }
                }
            }
        }

        if (gMin(fluidToSolidPointMapsPtrList_[i]) == -1)
        {
            FatalErrorIn(type() + "::calcFluidToSolidPointMaps() const")
                << "Cannot calculate the map between interfaces for global patch "
                << fluidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::fluidToSolidPointMaps() const
{
    if (fluidToSolidPointMapsPtrList_.empty())
    {
        calcFluidToSolidPointMaps();
    }

    return fluidToSolidPointMapsPtrList_;
}


void Foam::fluidSolidInterface::calcSolidToFluidPointMaps() const
{
    solidToFluidPointMapsPtrList_.setSize(nGlobalPatches());

    forAll(solid().globalPatches(), i)
    {
        const label patchID = solid().globalPatches()[i].patch().index();

        if (solidToFluidPointMapsPtrList_.set(i))
        {
            FatalErrorIn(type() + "::calcSolidToFluidPointMaps() const")
                << "pointer already set for global patch: "
                << solidMesh().boundary()[patchID].name() << "!"
                << abort(FatalError);
        }

        if
        (
            solid().globalPatches()[i].globalPatch().size()
         != fluid().globalPatches()[i].globalPatch().size()
        )
        {
            FatalErrorIn(type() + "::calcSolidToFluidPointMaps() const")
                << "Solid and fluid interfaces are not conformal for "
                << "global patch: "
                << solidMesh().boundary()[patchID].name() << "!" << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word solidToFluidPointMapName
        (
            "solidToFluidPointMap"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
        );

        IOobject mapHeader
        (
            solidToFluidPointMapName,
            solid().runTime().timeName(),
            solidMesh(),
            IOobject::MUST_READ
        );

        if (mapHeader.headerOk())
        {
            // Read map
            Info<< "Reading solidToFluidPointMap for global patch "
                << solidMesh().boundary()[patchID].name()
                << " from disk..." << endl;

            solidToFluidPointMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        solidToFluidPointMapName,
                        solid().runTime().timeName(),
                        solidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            // Initialise map
            solidToFluidPointMapsPtrList_.set
            (
                i,
                new labelIOList
                (
                    IOobject
                    (
                        solidToFluidPointMapName,
                        solid().runTime().timeName(),
                        solidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    labelList(fluid().globalPatches()[i].globalPatch().nPoints(), -1)
                )
            );
            labelList& solidToFluidMap = solidToFluidPointMapsPtrList_[i];

            // Perform N^2 search for corresponding points
            // We will take 0.1% of the minEdgeLength as the exact match tolerance

            const vectorField& fluidLP =
                fluid().globalPatches()[i].globalPatch().localPoints();

            const vectorField& solidLP =
                solid().globalPatches()[i].globalPatch().localPoints();

            const scalar tol = 0.001*gMin(minEdgeLengths()[i]);

            forAll(fluidLP, fluidPointI)
            {
                const vector& curFluidLP = fluidLP[fluidPointI];

                forAll(solidLP, solidPointI)
                {
                    if (mag(curFluidLP - solidLP[solidPointI]) < tol)
                    {
                        solidToFluidMap[fluidPointI] = solidPointI;
                        break;
                    }
                }
            }
        }

        if (gMin(solidToFluidPointMapsPtrList_[i]) == -1)
        {
            FatalErrorIn(type() + "::calcSolidToFluidPointMaps() const")
                << "Cannot calculate the map between interfaces for global patch "
                << solidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::solidToFluidPointMaps() const
{
    if (solidToFluidPointMapsPtrList_.empty())
    {
        calcSolidToFluidPointMaps();
    }

    return solidToFluidPointMapsPtrList_;

}


void Foam::fluidSolidInterface::
calcAccumulatedFluidInterfacesDisplacements() const
{
    accumulatedFluidInterfacesDisplacementsPtrList_.setSize
    (
        nGlobalPatches()
    );

    forAll(fluid().globalPatches(), i)
    {
        const label patchID = fluid().globalPatches()[i].patch().index();

        // Read accumulated displacement
        if (accumulatedFluidInterfacesDisplacementsPtrList_.set(i))
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::"
                "calcAccumulatedFluidInterfacesDisplacements() const"
            )   << "Accumulated displacement field already exists for global patch: "
                << fluidMesh().boundary()[patchID].name()
                << "!" << abort(FatalError);
        }

        const word accumulatedFluidInterfaceDisplacementName
        (
            "accumulatedFluidInterfaceDisplacement"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
        );

        // Accumulated fluid interface displacement
        IOobject accumulatedFluidInterfaceDisplacementHeader
        (
            accumulatedFluidInterfaceDisplacementName,
            fluid().runTime().timeName(),
            fluidMesh(),
            IOobject::MUST_READ
        );

        if (accumulatedFluidInterfaceDisplacementHeader.headerOk())
        {
            Pout<< "Reading accumulated fluid interface "
                << "displacement for global patch "
                << fluidMesh().boundary()[patchID].name()
                << " from disk..." << endl;

            accumulatedFluidInterfacesDisplacementsPtrList_.set
            (
                i,
                new vectorIOField
                (
                    IOobject
                    (
                        accumulatedFluidInterfaceDisplacementName,
                        fluid().runTime().timeName(),
                        fluidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            Pout<< "Creating accumulated fluid interface "
                << "displacement for global patch "
                << fluidMesh().boundary()[patchID].name()
                << " from disk..." << endl;

            accumulatedFluidInterfacesDisplacementsPtrList_.set
            (
                i,
                new vectorIOField
                (
                    IOobject
                    (
                        accumulatedFluidInterfaceDisplacementName,
                        fluid().runTime().timeName(),
                        fluidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    vectorField
                    (
                        fluidMesh().boundaryMesh()[fluidPatchIndices()[i]].nPoints(),
                        vector::zero
                    )
                )
            );
        }
    }
}


Foam::PtrList<Foam::vectorIOField>&
Foam::fluidSolidInterface::accumulatedFluidInterfacesDisplacements()
{
    if (accumulatedFluidInterfacesDisplacementsPtrList_.empty())
    {
        calcAccumulatedFluidInterfacesDisplacements();
    }

    return accumulatedFluidInterfacesDisplacementsPtrList_;
}


void Foam::fluidSolidInterface::calcMinEdgeLengths() const
{
    minEdgeLengthsPtrList_.setSize(nGlobalPatches());

    forAll(fluid().globalPatches(), i)
    {
        // Read accumulated displacement
        if (minEdgeLengthsPtrList_.set(i))
        {
            FatalErrorIn
            (
                "void fluidSolidInterface::"
                "calcMinEdgeLengths() const"
            )
                << "Minimal edge length already exists for global patch: "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << "!" << abort(FatalError);
        }

        minEdgeLengthsPtrList_.set
        (
            i,
            new scalarField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                0
            )
        );
        scalarField& minEdgeLength = minEdgeLengthsPtrList_[i];


        const edgeList& edges =
            fluid().globalPatches()[i].globalPatch().edges();

        const vectorField& points =
            fluid().globalPatches()[i].globalPatch().localPoints();

        const labelListList& pointEdges =
            fluid().globalPatches()[i].globalPatch().pointEdges();

        forAll(points, pointI)
        {
            const labelList& curPointEdges = pointEdges[pointI];

            scalar minLength = GREAT;

            forAll(curPointEdges, edgeI)
            {
                const edge& curEdge = edges[curPointEdges[edgeI]];

                scalar Le = curEdge.mag(points);

                if (Le < minLength)
                {
                    minLength = Le;
                }
            }

            minEdgeLength[pointI] = minLength;
        }

        // Pout << "Min edge length: " << min(minEdgeLength) << endl;
        // Pout << "gMin edge length: " << gMin(minEdgeLength) << endl;
    }
}


const Foam::PtrList<Foam::scalarField>&
Foam::fluidSolidInterface::minEdgeLengths() const
{
    if (minEdgeLengthsPtrList_.empty())
    {
        calcMinEdgeLengths();
    }

    return minEdgeLengthsPtrList_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::fluidSolidInterface
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime),
    IOdictionary
    (
        IOobject
        (
            "fsiProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    fsiProperties_(subDict(type + "Coeffs")),
    fluid_(fluidModel::New(runTime, "fluid")),
    solid_(solidModel::New(runTime, "solid")),
    solidPatchIndices_(),
    fluidPatchIndices_(),
    nGlobalPatches_(label(1)),
    currentSolidZonesPointsPtrList_(),
    currentSolidZonesPatchesPtrList_(),
    rbfFluidToSolidPtrList_(),
    rbfSolidToFluidPtrList_(),
    ggiInterpolatorsPtrList_(),
    fluidToSolidFaceMapsPtrList_(),
    solidToFluidFaceMapsPtrList_(),
    fluidToSolidPointMapsPtrList_(),
    solidToFluidPointMapsPtrList_(),
    outerCorrTolerance_
    (
        fsiProperties_.lookupOrDefault<scalar>("outerCorrTolerance", 1e-06)
    ),
    nOuterCorr_
    (
        fsiProperties_.lookupOrDefault<int>("nOuterCorr", 30)
    ),
    coupled_
    (
        fsiProperties_.lookupOrDefault<Switch>("coupled", true)
    ),
    couplingStartTime_
    (
        fsiProperties_.lookupOrDefault<scalar>("couplingStartTime", -1.0)
    ),
    predictor_(fsiProperties_.lookupOrDefault<Switch>("predictor", false)),
    interfaceDeformationLimit_
    (
        fsiProperties_.lookupOrDefault<scalar>("interfaceDeformationLimit", 0.0)
    ),
    fluidZonesPointsDispls_(),
    fluidZonesPointsDisplsRef_(),
    fluidZonesPointsDisplsPrev_(),
    solidZonesPointsDispls_(),
    solidZonesPointsDisplsRef_(),
    interfacesPointsDispls_(),
    interfacesPointsDisplsPrev_(),
    residuals_(),
    residualsPrev_(),
    maxResidualsNorm_(),
    maxIntsDisplsNorm_(),
    outerCorr_(0),
    writeResidualsToFile_
    (
        fsiProperties_.lookupOrDefault<Switch>("writeResidualsToFile", false)
    ),
    residualFilePtr_(),
    interpolatorUpdateFrequency_
    (
        fsiProperties_.lookupOrDefault<int>("interpolatorUpdateFrequency", 0)
    ),
    accumulatedFluidInterfacesDisplacementsPtrList_(),
    minEdgeLengthsPtrList_(),
    transferMethod_
    (
        interfaceTransferMethodNames_
        [
            fsiProperties_.lookupOrDefault<word>
            (
                "interfaceTransferMethod", "GGI"
            )
        ]
    )
{
    // Check if couplingStartTime is specified
    if (couplingStartTime_ > SMALL)
    {
        if (coupled_)
        {
            WarningIn(type + "::fsiProperties(...)")
                << "When using the coupilngStartTime option, the coupled "
                << "option should be set to off: resetting coupled to off"
                << endl;

            coupled_ = false;
        }
    }

    // Read interface patches names for regions
    const wordList solidPatchNames(fsiProperties_.lookup("solidPatch"));
    const wordList fluidPatchNames(fsiProperties_.lookup("fluidPatch"));

    if
    (
        solidPatchNames.size() != fluidPatchNames.size()
    )
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Defined number of coupled fluid and solid patches "
            << "must be equal!"
            << abort(FatalError);
    }

    solidPatchIndices_.setSize(solidPatchNames.size(), label(-1));
    fluidPatchIndices_.setSize(fluidPatchNames.size(), label(-1));

    // loop over all coupled patches
    forAll(solidPatchNames, i)
    {
        // Solid patch index
        const polyPatchID solidPatch
        (
            solidPatchNames[i],
            solidMesh().boundaryMesh()
        );

        if (!solidPatch.active())
        {
            FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
                << "Solid patch name " << solidPatchNames[i] << " not found."
                << abort(FatalError);
        }

        solidPatchIndices_[i] = solidPatch.index();
    }

    // Create solid global patches
    solid().makeGlobalPatches(solidPatchNames);

    // loop over all coupled patches
    forAll(fluidPatchNames, i)
    {
        // Fluid patch index
        const polyPatchID fluidPatch
        (
            fluidPatchNames[i],
            fluidMesh().boundaryMesh()
        );

        if (!fluidPatch.active())
        {
            FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
                << "Fluid patch name " << fluidPatchNames[i] << " not found."
                << abort(FatalError);
        }

        fluidPatchIndices_[i] = fluidPatch.index();
    }

    // Create fluid global patches
    fluid().makeGlobalPatches(fluidPatchNames);

    // Note: perform one last check before going any further
    // in case something goes wrong
    if
    (
        fluid().globalPatches().size()
     != solid().globalPatches().size()
    )
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "The number of assigned global poly patches "
            << "differs on fluid and solid regions!"
            << abort(FatalError);
    }

    // Set the number of global poly patches: solid or fluid
    nGlobalPatches_ = fluid().globalPatches().size();

    // Set interface fields list size and initialize residual
    fluidZonesPointsDispls_.setSize(nGlobalPatches());
    fluidZonesPointsDisplsRef_.setSize(nGlobalPatches());
    fluidZonesPointsDisplsPrev_.setSize(nGlobalPatches());
    solidZonesPointsDispls_.setSize(nGlobalPatches());
    solidZonesPointsDisplsRef_.setSize(nGlobalPatches());
    interfacesPointsDispls_.setSize(nGlobalPatches());
    interfacesPointsDisplsPrev_.setSize(nGlobalPatches());
    residuals_.setSize(nGlobalPatches());
    residualsPrev_.setSize(nGlobalPatches());
    maxResidualsNorm_.setSize(nGlobalPatches());
    maxIntsDisplsNorm_.setSize(nGlobalPatches());

    forAll(residuals_, i)
    {
        residuals_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(), vector::zero
            );
    }

    // Check if deprecated option rbfInterpolation is specified
    if (fsiProperties_.found("rbfInterpolation"))
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "The 'rbfInterpolation' is deprecated: instead please use the "
            << "'transferMethod' to specify the approach" << abort(FatalError);
    }

    if (transferMethod_ == directMap)
    {
        // Force maps to be created, in case of restart, they need to be read
        fluidToSolidFaceMaps();
        solidToFluidFaceMaps();
        fluidToSolidPointMaps();
        solidToFluidPointMaps();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::~fluidSolidInterface()
{
    currentSolidZonesPointsPtrList_.clear();
    currentSolidZonesPatchesPtrList_.clear();
    ggiInterpolatorsPtrList_.clear();
    accumulatedFluidInterfacesDisplacementsPtrList_.clear();
    minEdgeLengthsPtrList_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidSolidInterface> Foam::fluidSolidInterface::New
(
    Time& runTime,
    const word& region
)
{
    word fsiTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary fsiProperties
        (
            IOobject
            (
                "fsiProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        fsiProperties.lookup("fluidSolidInterface")
            >> fsiTypeName;
    }

    Info<< "Selecting fluidSolidInterface " << fsiTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fsiTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fluidSolidInterface::New(Time&, const word&)"
        )   << "Unknown fluidSolidInterface type " << fsiTypeName
            << endl << endl
            << "Valid fluidSolidInterface types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<fluidSolidInterface>(cstrIter()(runTime, region));
}


Foam::OFstream& Foam::fluidSolidInterface::residualFile()
{
    if (Pstream::parRun())
    {
        if (!Pstream::master())
        {
            FatalErrorIn
            (
                "Foam::OFstream& Foam::fluidSolidInterface::residualFile()"
            )   << "Only the master processor can call this functon!"
                << abort(FatalError);
        }
    }

    if (residualFilePtr_.empty())
    {
        const fileName historyDir = runTime().path()/"residuals";
        mkDir(historyDir);
        residualFilePtr_.set(new OFstream(historyDir/"fsiResiduals.dat"));
        residualFilePtr_()
            << "Time outerCorrector residual" << endl;
    }

    return residualFilePtr_();
}


const Foam::PtrList<Foam::vectorField>&
Foam::fluidSolidInterface::currentSolidZonesPoints() const
{
    if (currentSolidZonesPointsPtrList_.empty())
    {
        calcCurrentSolidZonesPoints();
    }

    return currentSolidZonesPointsPtrList_;
}


const Foam::PtrList<Foam::standAlonePatch>&
Foam::fluidSolidInterface::currentSolidZonesPatches() const
{
    if (currentSolidZonesPatchesPtrList_.empty())
    {
        calcCurrentSolidZonesPatches();
    }

    return currentSolidZonesPatchesPtrList_;
}


const Foam::List<std::shared_ptr<RBFInterpolation> >&
Foam::fluidSolidInterface::rbfFluidToSolid() const
{
    if (rbfFluidToSolidPtrList_.empty())
    {
        calcRbfFluidToSolidInterpolators();
    }

    return rbfFluidToSolidPtrList_;
}


const Foam::List<std::shared_ptr<RBFInterpolation> >&
Foam::fluidSolidInterface::rbfSolidToFluid() const
{
    if (rbfSolidToFluidPtrList_.empty())
    {
        calcRbfSolidToFluidInterpolators();
    }

    return rbfSolidToFluidPtrList_;
}


const Foam::PtrList<Foam::GGIInterpolation<standAlonePatch, standAlonePatch> >&
Foam::fluidSolidInterface::ggiInterpolators() const
{
    if (ggiInterpolatorsPtrList_.empty())
    {
        calcGgiInterpolators();
    }

    return ggiInterpolatorsPtrList_;
}


void Foam::fluidSolidInterface::setDeltaT(Time& runTime)
{
    // For now, the fluid sets the time-step
    fluid().setDeltaT(runTime);
}


void Foam::fluidSolidInterface::initializeFields()
{
    outerCorr_ = 0;

    forAll(fluid().globalPatches(), i)
    {
        fluidZonesPointsDispls_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        fluidZonesPointsDisplsRef_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        fluidZonesPointsDisplsPrev_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        solidZonesPointsDispls_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        solidZonesPointsDisplsRef_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        residualsPrev_[i] = residuals_[i];

        residuals_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        maxResidualsNorm_[i] = scalar(0);

        interfacesPointsDispls_[i] =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );

        interfacesPointsDisplsPrev_ =
            vectorField
            (
                fluid().globalPatches()[i].globalPatch().nPoints(),
                vector::zero
            );
    }
}


void Foam::fluidSolidInterface::updateInterpolatorAndGlobalPatches()
{
    if (ggiInterpolatorsPtrList_.empty())
    {
        ggiInterpolators();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
            ggiInterpolatorsPtrList_.clear();
            fluid().clearGlobalPatches();
            solid().clearGlobalPatches();

            const wordList fluidPatchNames(fsiProperties_.lookup("fluidPatch"));
            fluid().makeGlobalPatches(fluidPatchNames);

            const wordList solidPatchNames(fsiProperties_.lookup("solidPatch"));
            solid().makeGlobalPatches(solidPatchNames);

            // update ggi interpolators
            ggiInterpolators();
        }
    }
}


void Foam::fluidSolidInterface::moveFluidMesh()
{
    // Get fluid patch displacement from fluid zone displacement
    // Take care: these are local patch fields not global patch fields

    List<vectorField> fluidPatchesPointsDispls
    (
        nGlobalPatches(), vectorField()
    );

    List<vectorField> fluidPatchesPointsDisplsPrev
    (
        nGlobalPatches(), vectorField()
    );

    scalar maxDelta = 0;

    forAll(fluid().globalPatches(), i)
    {
        fluidPatchesPointsDispls[i] =
            fluid().globalPatches()[i].globalPointToPatch
            (
                fluidZonesPointsDispls()[i]
            );

        fluidPatchesPointsDisplsPrev[i] =
            fluid().globalPatches()[i].globalPointToPatch
            (
                fluidZonesPointsDisplsPrev()[i]
            );

        // Patch point normals
        const vectorField& n =
            fluid().mesh().boundaryMesh()
            [
                fluid().globalPatches()[i].patch().index()
            ].pointNormals();

        // Patch deltaCoeffs
        const scalarField fluidZoneDeltaCoeffs =
            fluid().globalPatches()[i].patchFaceToGlobal
            (
                fluidMesh().boundary()
                [
                    fluid().globalPatches()[i].patch().index()
                ].deltaCoeffs()
            );

        // Zone deltaCoeffs at points
        const scalarField fluidZonePointDeltaCoeffs =
            fluid().globalPatches()[i].interpolator().faceToPointInterpolate
            (
                fluidZoneDeltaCoeffs
            );

        // Patch deltaCoeffs at points
        const scalarField fluidPatchPointDeltaCoeffs =
            fluid().globalPatches()[i].globalPointToPatch
            (
                fluidZonePointDeltaCoeffs
            );

        const scalar delta =
            gMax
            (
                mag
                (
                    n
                  & (
                        accumulatedFluidInterfacesDisplacements()[i]
                      + fluidPatchesPointsDispls[i]
                      - fluidPatchesPointsDisplsPrev[i]
                    )
                )*fluidPatchPointDeltaCoeffs
            );

        if (delta > maxDelta)
        {
            maxDelta = delta;
        }
    }

    Info<< "Maximal accumulated displacement of interface points: "
        << maxDelta << endl;

    if (maxDelta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().allPoints();

        forAll(fluid().globalPatches(), i)
        {
            const labelList& meshPoints =
                fluid().globalPatches()[i].globalPatch().meshPoints();

            forAll(fluidPatchesPointsDispls[i], pointI)
            {
                newPoints[meshPoints[pointI]] +=
                    fluidPatchesPointsDispls[i][pointI]
                  - fluidPatchesPointsDisplsPrev[i][pointI];
            }

            twoDPointCorrector twoDCorrector(fluidMesh());
            twoDCorrector.correctPoints(newPoints);

            fluidMesh().movePoints(newPoints);

            // Accumulate interface points displacement
            accumulatedFluidInterfacesDisplacements()[i] +=
                fluidPatchesPointsDispls[i]
              - fluidPatchesPointsDisplsPrev[i];
        }
    }
    else
    {
        // Move whole fluid mesh

        // Check mesh motion solver type

        // If the motionU field is in the object registry then we assume that
        // the fe motion solver is being used
        const bool feMotionSolver =
            fluidMesh().foundObject<tetPointVectorField>("motionU");

        // If the pointMotionU field is in the object registry then we assume
        // that the fv motion solver is being used
        const bool fvMotionSolver =
            fluidMesh().foundObject<pointVectorField>("pointMotionU");

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            forAll(fluid().globalPatches(), i)
            {
                fixedValueTetPolyPatchVectorField& motionUFluidPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchIndices()[i]]
                    );

                tetPolyPatchInterpolation tppi
                (
                    refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
                );

                motionUFluidPatch ==
                    tppi.pointToPointInterpolate
                    (
                        accumulatedFluidInterfacesDisplacements()[i]
                      / fluid().runTime().deltaT().value()
                    );
            }
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            forAll(fluid().globalPatches(), i)
            {
                fixedValuePointPatchVectorField& motionUFluidPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchIndices()[i]]
                    );

                motionUFluidPatch ==
                    (
                        fluidPatchesPointsDispls[i]
                      - fluidPatchesPointsDisplsPrev[i]
                    )/fluid().runTime().deltaT().value();
            }
        }
        else if (isA<newSubsetMotionSolverFvMesh>(fluidMesh()))
        {
            newSubsetMotionSolverFvMesh& dynMesh =
                refCast<newSubsetMotionSolverFvMesh>
                (
                    fluidMesh()
                );

            const fvMesh& subMesh = dynMesh.subsetMesh().subMesh();
            
            const bool fvMotionSolver =
                subMesh.foundObject<pointVectorField>("pointMotionU");

            // Info << subMesh.boundaryMesh() << endl;
            
            if (fvMotionSolver)
            {
                pointVectorField& motionU =
                    const_cast<pointVectorField&>
                    (
                        subMesh.objectRegistry::
                        lookupObject<pointVectorField>
                        (
                            "pointMotionU"
                        )
                    );

                forAll(fluid().globalPatches(), i)
                {
                    fixedValuePointPatchVectorField& motionUFluidPatch =
                        refCast<fixedValuePointPatchVectorField>
                        (
                            motionU.boundaryField()[fluidPatchIndices()[i]]
                        );

                    motionUFluidPatch ==
                    (
                        fluidPatchesPointsDispls[i]
                      - fluidPatchesPointsDisplsPrev[i]
                    )/fluid().runTime().deltaT().value();

                    // FatalErrorIn("fluidSolidInterface::moveFluidMesh()")
                    //   << "subset fvMotionSolver"
                    //   << abort(FatalError);
                }
            }
        }
        else
        {
            FatalErrorIn("fluidSolidInterface::moveFluidMesh()")
                << "Problem with fluid mesh motion solver selection"
                << abort(FatalError);
        }

        bool meshChanged = fluidMesh().update();
        reduce(meshChanged, orOp<bool>());
        fluid().fsiMeshUpdate() = true;
        fluid().fsiMeshUpdateChanged() = meshChanged;

        forAll(fluid().globalPatches(), i)
        {
            accumulatedFluidInterfacesDisplacements()[i] =
                vectorField
                (
                    accumulatedFluidInterfacesDisplacements()[i].size(),
                    vector::zero
                );
        }
    }


    // Move unused fluid mesh points
    // Not needed anymore as globalFaceZones are not used
    // {
    //     vectorField newPoints = fluidMesh().allPoints();

    //     const labelList& fluidZoneMeshPoints =
    //         fluid().globalPatch().globalPatch().meshPoints();

    //     forAll(fluidZonePointsDispl(), pointI)
    //     {
    //         if (fluidZoneMeshPoints[pointI] >= fluidMesh().nPoints())
    //         {
    //             newPoints[fluidZoneMeshPoints[pointI]] +=
    //                 fluidZonePointsDispl()[pointI]
    //               - fluidZonePointsDisplPrev()[pointI];
    //         }
    //     }

    //     twoDPointCorrector twoDCorrector(fluidMesh());

    //     twoDCorrector.correctPoints(newPoints);

    //     fluidMesh().movePoints(newPoints);
    // }
}


void Foam::fluidSolidInterface::updateForce()
{
    Info<< "Setting traction on solid patch/patches" << endl;

    const List<tmp<vectorField> > faceZonesViscousForce
    (
        fluid().faceZonesViscousForce()
    );

    const List<tmp<scalarField> > faceZonesPressureForce
    (
        fluid().faceZonesPressureForce()
    );

    List<vectorField> solidZonesTotalTraction
    (
        nGlobalPatches(), vectorField()
    );

    List<vectorField> fluidZonesTotalTraction
    (
        nGlobalPatches(), vectorField()
    );

    forAll(fluid().globalPatches(), i)
    {
        const vectorField& fluidZoneTraction = faceZonesViscousForce[i]();

        const scalarField& fluidZonePressure = faceZonesPressureForce[i]();

        // Fluid zone face normals
        const vectorField& n =
            fluid().globalPatches()[i].globalPatch().faceNormals();

        // Fluid zone total traction
        fluidZonesTotalTraction[i] =
            fluidZoneTraction - fluidZonePressure*n;

        // Solid zone traction is to be interpolated from the fluid zone
        solidZonesTotalTraction[i] =
            vectorField
            (
                solid().globalPatches()[i].globalPatch().size(),
                vector::zero
            );
    }

    transferFacesZoneToZone
    (
        "fluid",                 // from region name
        "solid",                 // to region name
        fluid().globalPatches(), // from zones
        solid().globalPatches(), // to zones
        fluidZonesTotalTraction, // from fields
        solidZonesTotalTraction  // to fields
    );

    // Flip traction sign after transferring from fluid to solid
    forAll(fluid().globalPatches(), i)
    {
        solidZonesTotalTraction[i] = -solidZonesTotalTraction[i];
    }

    // Update coupling
    if (!coupled_)
    {
        updateCoupled();
    }

    // Set traction on solid
    if (coupled())
    {
        solid().setTraction
        (
            solidPatchIndices(),
            solidZonesTotalTraction
        );

        // Set interface pressure for elasticWallPressure
        // boundary condition
        forAll(fluid().globalPatches(), i)
        {
            if
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()[fluidPatchIndices()[i]]
                )
            )
            {
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().p().boundaryField()[fluidPatchIndices()[i]]
                    )
                ).prevPressure() = fluid().patchPressureForce(fluidPatchIndices()[i]);
            }
        }
    }
    else
    {
        // Set interface pressure for elasticWallPressure
        // boundary condition
        forAll(fluid().globalPatches(), i)
        {
            if
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()[fluidPatchIndices()[i]]
                )
            )
            {
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().p().boundaryField()[fluidPatchIndices()[i]]
                    )
                ).prevPressure() = 0;
            }
        }
    }

    // Total force at the fluid side of the interface
    {
        forAll(fluid().globalPatches(), i)
        {
            const vectorField& p =
                fluid().globalPatches()[i].globalPatch().localPoints();
            const faceList& f =
                fluid().globalPatches()[i].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);
            vectorField C(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
                C[faceI] = f[faceI].centre(p);
            }

            const vector totalTractionForce = sum(fluidZonesTotalTraction[i]*mag(S));

            Info<< "Total force on interface patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << " (fluid) = " << totalTractionForce << endl;
        }
    }

    // Total force at the solid side of the interface
    {
        forAll(solid().globalPatches(), i)
        {
            const vectorField& p =
                solid().globalPatches()[i].globalPatch().localPoints();
            const faceList& f =
                solid().globalPatches()[i].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);
            vectorField C(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
                C[faceI] = f[faceI].centre(p);
            }

            const vector totalTractionForce = sum(solidZonesTotalTraction[i]*mag(S));

            Info<< "Total force on interface patch "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[i].patch().index()
                   ].name()
                << " (solid) = " << totalTractionForce << endl;
        }
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    List<scalar> minResidual(nGlobalPatches(), scalar(0));

    const List<tmp<vectorField> > faceZonesPointDisplacementIncrement
    (
        solid().faceZonesPointDisplacementIncrement()
    );

    const List<tmp<vectorField> > faceZonesPointDisplacementOld
    (
        solid().faceZonesPointDisplacementOld()
    );

    List<vectorField> solidZonesPointsDisplsAtSolid
    (
        nGlobalPatches(),
        vectorField()
    );

    List<vectorField> solidZonesPointsTotDisplsAtSolid
    (
        nGlobalPatches(),
        vectorField()
    );

    List<vectorField> solidZonesPointsTotDispls
    (
        nGlobalPatches(),
        vectorField()
    );

    forAll(solid().globalPatches(), i)
    {
        solidZonesPointsDisplsAtSolid[i] =
            faceZonesPointDisplacementIncrement[i]();

        solidZonesPointsTotDisplsAtSolid[i] =
            faceZonesPointDisplacementOld[i]();

        solidZonesPointsTotDispls[i] =
            vectorField
            (
                solidZonesPointsDispls()[i].size(),
                vector::zero
            );
    }

    // Transfer displacement fields from the solid to the fluid
    transferPointsZoneToZone
    (
        "solid",                       // from region name
        "fluid",                       // to region name
        solid().globalPatches(),       // from zones
        fluid().globalPatches(),       // to zones
        solidZonesPointsDisplsAtSolid, // from fields
        solidZonesPointsDispls()       // to fields
    );

    transferPointsZoneToZone
    (
        "solid",                          // from region name
        "fluid",                          // to region name
        solid().globalPatches(),          // from zones
        fluid().globalPatches(),          // to zones
        solidZonesPointsTotDisplsAtSolid, // from fields
        solidZonesPointsTotDispls         // to fields
    );

    forAll(solid().globalPatches(), i)
    {
        // Update interface residuals
        residualsPrev()[i] = residuals()[i];

        residuals()[i] =
            solidZonesPointsDispls()[i] - fluidZonesPointsDispls()[i];

        scalar residualNorm = Foam::sqrt(gSum(magSqr(residuals()[i])));
        scalar residualNorm_2 = residualNorm;

        if (residualNorm > maxResidualsNorm_[i])
        {
            maxResidualsNorm_[i] = residualNorm;
        }

        residualNorm /= maxResidualsNorm_[i] + SMALL;

        Info<< "Current fsi relative residual norm ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[i].patch().index()
               ].name()
            << "): " << residualNorm << endl;

        interfacesPointsDisplsPrev_[i] = interfacesPointsDispls_[i];

        interfacesPointsDispls_[i] = solidZonesPointsDispls()[i];

        const vectorField intTotDispl =
            interfacesPointsDispls_[i] + solidZonesPointsTotDispls[i];

        const scalar intTotDisplNorm = Foam::sqrt(gSum(magSqr(intTotDispl)));

        if (intTotDisplNorm > maxIntsDisplsNorm_[i])
        {
            maxIntsDisplsNorm_[i] = intTotDisplNorm;
        }

        residualNorm_2 /= maxIntsDisplsNorm_[i] + SMALL;

        Info<< "Alternative fsi residual ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[i].patch().index()
               ].name()
            << "): " << residualNorm_2 << endl;

        minResidual[i] = min(residualNorm_2, residualNorm);
    }

    return max(minResidual);
}


void Foam::fluidSolidInterface::updateMovingWallPressureAcceleration()
{
    const List<tmp<vectorField> > faceZonesAcceleration
    (
        solid().faceZonesAcceleration()
    );

    forAll(fluid().globalPatches(), i)
    {
        if
        (
            isA<movingWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndices()[i]]
            )
        )
        {
            Info<< "Setting acceleration at fluid side of the interface" << endl;

            const vectorField& solidZoneAcceleration = faceZonesAcceleration[i]();

            const vectorField fluidZoneAcceleration =
                ggiInterpolators()[i].slaveToMaster(solidZoneAcceleration);

            const vectorField fluidPatchAcceleration =
                fluid().globalPatches()[i].globalFaceToPatch(fluidZoneAcceleration);

            const_cast<movingWallPressureFvPatchScalarField&>
            (
                refCast<const movingWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()
                    [
                        fluid().globalPatches()[i].patch().index()
                    ]
                )
            ).prevAcceleration() = fluidPatchAcceleration;
        }
    }
}


void Foam::fluidSolidInterface::updateElasticWallPressureAcceleration()
{
    const List<tmp<vectorField> > faceZonesAcceleration
    (
        solid().faceZonesAcceleration()
    );

    forAll(fluid().globalPatches(), i)
    {
        // Set interface acceleration
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndices()[i]]
            )
        )
        {
            const vectorField& solidZoneAcceleration = faceZonesAcceleration[i]();

            const vectorField fluidZoneAcceleration =
                ggiInterpolators()[i].slaveToMaster(solidZoneAcceleration);

            const vectorField fluidPatchAcceleration =
                fluid().globalPatches()[i].globalFaceToPatch(fluidZoneAcceleration);

            const_cast<elasticWallPressureFvPatchScalarField&>
            (
                refCast<const elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()
                    [
                        fluid().globalPatches()[i].patch().index()
                    ]
                )
            ).prevAcceleration() = fluidPatchAcceleration;
        }
    }
}


void Foam::fluidSolidInterface::syncFluidZonePointsDispl
(
    List<vectorField>& fluidZonesPointsDispls
)
{
    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        forAll(fluid().globalPatches(), i)
        {
            if (!Pstream::master())
            {
                fluidZonesPointsDispls[i] = vector::zero;
            }

            // pass to all procs
            reduce(fluidZonesPointsDispls[i], sumOp<vectorField>());

            const labelList& map =
                fluid().globalPatches()[i].globalMasterToCurrentProcPointAddr();

            if (!Pstream::master())
            {
                const vectorField fluidZonePointsDisplGlobal = fluidZonesPointsDispls[i];

                forAll(fluidZonePointsDisplGlobal, globalPointI)
                {
                    const label localPoint = map[globalPointI];

                    fluidZonesPointsDispls[i][localPoint] =
                        fluidZonePointsDisplGlobal[globalPointI];
                }
            }
        }
    }
}

void Foam::fluidSolidInterface::writeFields(const Time& runTime)
{
    // solid calls runTime.write() to write both solid and fluid fields
    //fluid().writeFields(runTime);
    solid().writeFields(runTime);
}

// ************************************************************************* //
