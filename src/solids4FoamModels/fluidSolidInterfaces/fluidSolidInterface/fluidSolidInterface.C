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
    if (currentSolidZonesPointsList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentSolidZonesPoints() const"
        )   << "List already exists"
            << abort(FatalError);
    }

    currentSolidZonesPointsList_.setSize(nGlobalPatches_);

    forAll(solid().globalPatches(), interfaceI)
    {
        currentSolidZonesPointsList_.set
        (
            interfaceI,
            new vectorField(solid().currentFaceZonesPoints(interfaceI)())
        );
    }
}


void Foam::fluidSolidInterface::calcCurrentSolidZonesPatches() const
{
    if (currentSolidZonesPatchesList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "currentSolidZonesPatchesPtrList() const"
        )   << "List already exists"
            << abort(FatalError);
    }

    currentSolidZonesPatchesList_.setSize(nGlobalPatches_);

    forAll(solid().globalPatches(), interfaceI)
    {
        currentSolidZonesPatchesList_.set
        (
            interfaceI,
            new standAlonePatch
            (
                solid().globalPatches()[interfaceI].globalPatch().localFaces(),
                currentSolidZonesPoints()[interfaceI]
            )
        );
    }
}


void Foam::fluidSolidInterface::calcRbfFluidToSolidInterpolators() const
{
    if (rbfFluidToSolidList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcRbfFluidToSolidInterpolators() const"
        )   << "Fluid to solid interpolators already exist!"
            << abort(FatalError);
    }

    Info<< "Create RBF fluid-to-solid interpolators" << endl;

    rbfFluidToSolidList_.setSize(nGlobalPatches_);

    forAll(fluid().globalPatches(), interfaceI)
    {
        std::shared_ptr<RBFFunctionInterface> rbfFunction;
        rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

        rbfFluidToSolidList_[interfaceI] =
            std::shared_ptr<RBFInterpolation>
            (
                new RBFInterpolation(rbfFunction)
            );

        const vectorField solidZoneFaceCentres =
            currentSolidZonesPatches()[interfaceI].faceCentres();

        const vectorField fluidZoneFaceCentres =
            fluid().globalPatches()[interfaceI].globalPatch().faceCentres();

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

        rbfFluidToSolidList_[interfaceI]->compute(fluidX, solidX);

        // Check interpolation error

        matrix fluidXatSolid(solidZoneFaceCentres.size(), 3);

        rbfFluidToSolidList_[interfaceI]->interpolate(fluidX, fluidXatSolid);

        vectorField fluidFaceCentresAtSolid
        (
            solidZoneFaceCentres.size(), vector::zero
        );

        forAll(fluidFaceCentresAtSolid, faceI)
        {
            fluidFaceCentresAtSolid[faceI].x() = fluidXatSolid(faceI, 0);
            fluidFaceCentresAtSolid[faceI].y() = fluidXatSolid(faceI, 1);
            fluidFaceCentresAtSolid[faceI].z() = fluidXatSolid(faceI, 2);
        }

        const scalar maxDist = gMax
        (
            mag(fluidFaceCentresAtSolid - solidZoneFaceCentres)
        );

        Info<< "    face interpolation error for interface "
            << interfaceI << ": " << maxDist
            << endl;
    }
}


void Foam::fluidSolidInterface::calcRbfSolidToFluidInterpolators() const
{
    if (rbfSolidToFluidList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcRbfSolidToFluidInterpolators() const"
        )   << "Solid to fluid interpolators already exist!"
            << abort(FatalError);
    }

    Info<< "Create RBF solid-to-fluid interpolators" << endl;

    rbfSolidToFluidList_.setSize(nGlobalPatches_);

    forAll(solid().globalPatches(), interfaceI)
    {
        std::shared_ptr<RBFFunctionInterface> rbfFunction;
        rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

        rbfSolidToFluidList_[interfaceI] =
            std::shared_ptr<RBFInterpolation>
            (
                new RBFInterpolation(rbfFunction)
            );

        const vectorField solidZonePoints =
            currentSolidZonesPatches()[interfaceI].localPoints();

        const vectorField fluidZonePoints =
            fluid().globalPatches()[interfaceI].globalPatch().localPoints();

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

        rbfSolidToFluidList_[interfaceI]->compute(solidX, fluidX);

        // Check interpolation error

        matrix solidXatFluid(fluidZonePoints.size(), 3);

        rbfSolidToFluidList_[interfaceI]->interpolate
        (
            solidX, solidXatFluid
        );

        vectorField fluidZonePointsInterp(fluidZonePoints.size(), vector::zero);

        forAll(fluidZonePoints, faceI)
        {
            fluidZonePointsInterp[faceI].x() = solidXatFluid(faceI, 0);
            fluidZonePointsInterp[faceI].y() = solidXatFluid(faceI, 1);
            fluidZonePointsInterp[faceI].z() = solidXatFluid(faceI, 2);
        }

        const scalar maxDist = gMax
            (
                mag(fluidZonePointsInterp - fluidZonePoints)
            );

        Info<< "    point interpolation error for interface "
            << interfaceI << ": " << maxDist
            << endl;
    }
}


void Foam::fluidSolidInterface::calcGgiInterpolators() const
{
    if (ggiInterpolatorsList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcGgiInterpolator() const"
        )   << "GGI interpolators already exists"
            << abort(FatalError);
    }

    // Note: Remove current solid zone and points so that
    // it will be re-created in the deformed position
    currentSolidZonesPatchesList_.clear();
    currentSolidZonesPointsList_.clear();

    Info<< "Create GGI zone-to-zone interpolators"<< endl;

    ggiInterpolatorsList_.setSize(nGlobalPatches_);

    forAll(fluid().globalPatches(), interfaceI)
    {

        ggiInterpolatorsList_.set
        (
            interfaceI,
            new GGIInterpolation<standAlonePatch, standAlonePatch>
            (
                fluid().globalPatches()[interfaceI].globalPatch(),
                currentSolidZonesPatches()[interfaceI],
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

        {
            const vectorField fluidPatchFaceCentres =
                vectorField
                (
                    fluidMesh().boundaryMesh()
                    [
                        fluidPatchIndices()[interfaceI]
                    ].faceCentres()
                );

            const vectorField fluidZoneFaceCentres =
                fluid().globalPatches()[interfaceI].patchFaceToGlobal
                (
                    fluidPatchFaceCentres
                );

            const vectorField solidZoneFaceCentres =
                ggiInterpolatorsList_(interfaceI)->masterToSlave
                (
                    fluidZoneFaceCentres
                );

            const vectorField solidPatchFaceCentres =
                solid().globalPatches()[interfaceI].globalFaceToPatch
                (
                    solidZoneFaceCentres
                );

            scalar maxDist = gMax
            (
                mag
                (
                    solidPatchFaceCentres
                  - solidMesh().boundaryMesh()
                    [
                        solidPatchIndices()[interfaceI]
                    ].faceCentres()
                )
            );

            Info<< "    interface " << interfaceI
                << " fluid-to-solid error: " << maxDist
                << endl;
        }

        {
            const vectorField solidZonePointsAtFluid =
                ggiInterpolatorsList_(interfaceI)->slaveToMasterPointInterpolate
                (
                    currentSolidZonesPoints()[interfaceI]
                );

            const vectorField fluidZonePoints =
                fluid().globalPatches()[interfaceI].globalPatch().localPoints();

            const scalar maxDist = gMax
            (
                mag
                (
                    fluidZonePoints
                  - solidZonePointsAtFluid
                )
            );

            Info<< "    interface " << interfaceI
                << " solid-to-fluid error: " << maxDist
                << endl;
        }

        Info<< "        number of uncovered master faces: "
            << ggiInterpolatorsList_(interfaceI)->uncoveredMasterFaces().size() << endl;

        Info<< "        number of uncovered slave faces: "
            << ggiInterpolatorsList_(interfaceI)->uncoveredSlaveFaces().size() << endl;

        ggiInterpolatorsList_(interfaceI)->slavePointDistanceToIntersection();
        ggiInterpolatorsList_(interfaceI)->masterPointDistanceToIntersection();

        Info<< endl;
    }
}


void Foam::fluidSolidInterface::calcFluidToSolidFaceMaps() const
{
    if (fluidToSolidFaceMapsList_.size())
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMaps() const")
            << "List already set!" << abort(FatalError);
    }

    fluidToSolidFaceMapsList_.setSize(nGlobalPatches_);

    forAll(fluid().globalPatches(), interfaceI)
    {
        const label patchID =
            fluid().globalPatches()[interfaceI].patch().index();

        if
        (
            solid().globalPatches()[interfaceI].globalPatch().size()
         != fluid().globalPatches()[interfaceI].globalPatch().size()
        )
        {
            FatalErrorIn(type() + "::calcFluidToSolidFaceMaps() const")
                << "Fluid and solid interfaces are not conformal (fluid patch ="
                << " " << fluidMesh().boundary()[patchID].name() << ", solid "
                << "patch = " << solidMesh().boundary()[patchID].name() << ")"
                << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word fluidToSolidFaceMapName
        (
            "fluidToSolidFaceMap" + Foam::name(interfaceI)
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
                << " from disk" << endl;

            fluidToSolidFaceMapsList_.set
            (
                interfaceI,
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
            fluidToSolidFaceMapsList_.set
            (
                interfaceI,
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
                    labelList
                    (
                        solid().globalPatches()
                        [
                            interfaceI
                        ].globalPatch().size(),
                        -1
                    )
                )
            );
            labelList& fluidToSolidMap = fluidToSolidFaceMapsList_[interfaceI];

            // Perform N^2 search for corresponding faces
            // We will take 0.1% of the minEdgeLength as the exact match
            // tolerance

            const vectorField& fluidCf =
                fluid().globalPatches()[interfaceI].globalPatch().faceCentres();

            const vectorField& solidCf =
                solid().globalPatches()[interfaceI].globalPatch().faceCentres();

            const scalar tol = 0.001*gMin(minEdgeLengths()[interfaceI]);

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

        if (gMin(fluidToSolidFaceMapsList_[interfaceI]) == -1)
        {
            FatalErrorIn(type() + "::calcFluidToSolidFaceMaps() const")
                << "Cannot calculate the map between interfaces for global "
                << "patch "
                << fluidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::fluidToSolidFaceMaps() const
{
    if (fluidToSolidFaceMapsList_.empty())
    {
        calcFluidToSolidFaceMaps();
    }

    return fluidToSolidFaceMapsList_;
}


void Foam::fluidSolidInterface::calcSolidToFluidFaceMaps() const
{
    if (solidToFluidFaceMapsList_.size())
    {
        FatalErrorIn(type() + "::calcSolidToFluidFaceMaps() const")
            << "Maps already set!"
            << abort(FatalError);
    }

    solidToFluidFaceMapsList_.setSize(nGlobalPatches_);

    forAll(solid().globalPatches(), interfaceI)
    {
        const label patchID =
            solid().globalPatches()[interfaceI].patch().index();

        if
        (
            solid().globalPatches()[interfaceI].globalPatch().size()
         != fluid().globalPatches()[interfaceI].globalPatch().size()
        )
        {
            FatalErrorIn(type() + "::calcSolidToFluidFaceMaps() const")
                << "Fluid and solid interfaces are not conformal (fluid patch ="
                << " " << fluidMesh().boundary()[patchID].name() << ", solid "
                << "patch = " << solidMesh().boundary()[patchID].name() << ")"
                << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word solidToFluidFaceMapName
        (
            "solidToFluidFaceMap" + Foam::name(interfaceI)
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
                << " from disk" << endl;

            solidToFluidFaceMapsList_.set
            (
                interfaceI,
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
            solidToFluidFaceMapsList_.set
            (
                interfaceI,
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
                    labelList
                    (
                        fluid().globalPatches()
                        [
                            interfaceI
                        ].globalPatch().size(),
                        -1
                    )
                )
            );
            labelList& solidToFluidMap = solidToFluidFaceMapsList_[interfaceI];

            // Perform N^2 search for corresponding faces
            // We will take 0.1% of the minEdgeLength as the exact match
            // tolerance

            const vectorField& fluidCf =
                fluid().globalPatches()[interfaceI].globalPatch().faceCentres();

            const vectorField& solidCf =
                solid().globalPatches()[interfaceI].globalPatch().faceCentres();

            const scalar tol = 0.001*gMin(minEdgeLengths()[interfaceI]);

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

        if (gMin(solidToFluidFaceMapsList_[interfaceI]) == -1)
        {
            FatalErrorIn(type() + "::calcSolidToFluidFaceMaps() const")
                << "Cannot calculate the map between interfaces for global "
                << "patch "
                << solidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::solidToFluidFaceMaps() const
{
    if (solidToFluidFaceMapsList_.empty())
    {
        calcSolidToFluidFaceMaps();
    }

    return solidToFluidFaceMapsList_;
}


void Foam::fluidSolidInterface::calcFluidToSolidPointMaps() const
{
    if (fluidToSolidPointMapsList_.size())
    {
        FatalErrorIn(type() + "::calcFluidToSolidPointMaps() const")
            << "Maps already set!"
            << abort(FatalError);
    }

    fluidToSolidPointMapsList_.setSize(nGlobalPatches_);

    forAll(fluid().globalPatches(), interfaceI)
    {
        const label patchID =
            fluid().globalPatches()[interfaceI].patch().index();

        if
        (
            solid().globalPatches()[interfaceI].globalPatch().nPoints()
         != fluid().globalPatches()[interfaceI].globalPatch().nPoints()
        )
        {
            FatalErrorIn(type() + "::calcFluidToSolidPointMaps() const")
                << "Fluid and solid interfaces are not conformal (fluid patch ="
                << " " << fluidMesh().boundary()[patchID].name() << ", solid "
                << "patch = " << solidMesh().boundary()[patchID].name() << ")"
                << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word fluidToSolidPointMapName
        (
            "fluidToSolidPointMap" + Foam::name(interfaceI)
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
                << " from disk" << endl;

            fluidToSolidPointMapsList_.set
            (
                interfaceI,
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
            fluidToSolidPointMapsList_.set
            (
                interfaceI,
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
                    labelList
                    (
                        solid().globalPatches()
                        [
                            interfaceI
                        ].globalPatch().nPoints(),
                        -1
                    )
                )
            );
            labelList& fluidToSolidMap = fluidToSolidPointMapsList_[interfaceI];

            // Perform N^2 search for corresponding points
            // We will take 0.1% of the minEdgeLength as the exact match
            // tolerance

            const vectorField& fluidLP =
                fluid().globalPatches()[interfaceI].globalPatch().localPoints();

            const vectorField& solidLP =
                solid().globalPatches()[interfaceI].globalPatch().localPoints();

            const scalar tol = 0.001*gMin(minEdgeLengths()[interfaceI]);

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

        if (gMin(fluidToSolidPointMapsList_[interfaceI]) == -1)
        {
            FatalErrorIn(type() + "::calcFluidToSolidPointMaps() const")
                << "Cannot calculate the map between interfaces for global "
                << "patch "
                << fluidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::fluidToSolidPointMaps() const
{
    if (fluidToSolidPointMapsList_.empty())
    {
        calcFluidToSolidPointMaps();
    }

    return fluidToSolidPointMapsList_;
}


void Foam::fluidSolidInterface::calcSolidToFluidPointMaps() const
{
    if (solidToFluidPointMapsList_.size())
    {
        FatalErrorIn(type() + "::calcSolidToFluidPointMaps() const")
            << "Map already set!"
            << abort(FatalError);
    }

    solidToFluidPointMapsList_.setSize(nGlobalPatches_);

    forAll(solid().globalPatches(), interfaceI)
    {
        const label patchID =
            solid().globalPatches()[interfaceI].patch().index();

        if
        (
            solid().globalPatches()[interfaceI].globalPatch().nPoints()
         != fluid().globalPatches()[interfaceI].globalPatch().nPoints()
        )
        {
            FatalErrorIn(type() + "::calcSolidToFluidPointMaps() const")
                << "Fluid and solid interfaces are not conformal (fluid patch ="
                << " " << fluidMesh().boundary()[patchID].name() << ", solid "
                << "patch = " << solidMesh().boundary()[patchID].name() << ")"
                << nl
                << "directMap method requires conformal interfaces."
                << abort(FatalError);
        }

        const word solidToFluidPointMapName
        (
            "solidToFluidPointMap" + Foam::name(interfaceI)
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
                << " from disk" << endl;

            solidToFluidPointMapsList_.set
            (
                interfaceI,
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
            solidToFluidPointMapsList_.set
            (
                interfaceI,
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
                    labelList
                    (
                        fluid().globalPatches()
                        [
                            interfaceI
                        ].globalPatch().nPoints(),
                        -1
                    )
                )
            );
            labelList& solidToFluidMap = solidToFluidPointMapsList_[interfaceI];

            // Perform N^2 search for corresponding points
            // We will take 0.1% of the minEdgeLength as the exact match
            // tolerance

            const vectorField& fluidLP =
                fluid().globalPatches()[interfaceI].globalPatch().localPoints();

            const vectorField& solidLP =
                solid().globalPatches()[interfaceI].globalPatch().localPoints();

            const scalar tol = 0.001*gMin(minEdgeLengths()[interfaceI]);

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

        if (gMin(solidToFluidPointMapsList_[interfaceI]) == -1)
        {
            FatalErrorIn(type() + "::calcSolidToFluidPointMaps() const")
                << "Cannot calculate the map between interfaces for global "
                << "patch "
                << solidMesh().boundary()[patchID].name() << "." << nl
                << "Direct mapping can only be used with conformal interfaces!"
                << abort(FatalError);
        }
    }
}


const Foam::PtrList<Foam::labelIOList>&
Foam::fluidSolidInterface::solidToFluidPointMaps() const
{
    if (solidToFluidPointMapsList_.empty())
    {
        calcSolidToFluidPointMaps();
    }

    return solidToFluidPointMapsList_;

}


void Foam::fluidSolidInterface::
calcAccumulatedFluidInterfacesDisplacements() const
{
    if (accumulatedFluidInterfacesDisplacementsList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcAccumulatedFluidInterfacesDisplacements() const"
        )   << "List already exists!" << abort(FatalError);
    }

    accumulatedFluidInterfacesDisplacementsList_.setSize
    (
        nGlobalPatches_
    );

    forAll(fluid().globalPatches(), interfaceI)
    {
        const label patchID =
            fluid().globalPatches()[interfaceI].patch().index();

        const word accumulatedFluidInterfaceDisplacementName
        (
            "accumulatedFluidInterfaceDisplacement" + Foam::name(interfaceI)
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
                << " from disk" << endl;

            accumulatedFluidInterfacesDisplacementsList_.set
            (
                interfaceI,
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
            accumulatedFluidInterfacesDisplacementsList_.set
            (
                interfaceI,
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
                        fluidMesh().boundaryMesh()
                        [
                            fluidPatchIndices()[interfaceI]
                        ].nPoints(),
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
    if (accumulatedFluidInterfacesDisplacementsList_.empty())
    {
        calcAccumulatedFluidInterfacesDisplacements();
    }

    return accumulatedFluidInterfacesDisplacementsList_;
}


void Foam::fluidSolidInterface::calcMinEdgeLengths() const
{
    if (minEdgeLengthsList_.size())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcMinEdgeLengths() const"
        )   << "List a;ready exists!" << abort(FatalError);
    }

    minEdgeLengthsList_.setSize(nGlobalPatches_);

    forAll(fluid().globalPatches(), interfaceI)
    {
        minEdgeLengthsList_.set
        (
            interfaceI,
            new scalarField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                0
            )
        );
        scalarField& minEdgeLength = minEdgeLengthsList_[interfaceI];

        const edgeList& edges =
            fluid().globalPatches()[interfaceI].globalPatch().edges();

        const vectorField& points =
            fluid().globalPatches()[interfaceI].globalPatch().localPoints();

        const labelListList& pointEdges =
            fluid().globalPatches()[interfaceI].globalPatch().pointEdges();

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
    if (minEdgeLengthsList_.empty())
    {
        calcMinEdgeLengths();
    }

    return minEdgeLengthsList_;
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
    solidPatchNames_(),
    fluidPatchNames_(),
    solidPatchIndices_(),
    fluidPatchIndices_(),
    nGlobalPatches_(-1),
    currentSolidZonesPointsList_(),
    currentSolidZonesPatchesList_(),
    rbfFluidToSolidList_(),
    rbfSolidToFluidList_(),
    ggiInterpolatorsList_(),
    fluidToSolidFaceMapsList_(),
    solidToFluidFaceMapsList_(),
    fluidToSolidPointMapsList_(),
    solidToFluidPointMapsList_(),
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
    accumulatedFluidInterfacesDisplacementsList_(),
    minEdgeLengthsList_(),
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
    // To maintain backwards compatibility, we will first check if a single
    // interface pair are defined, then we will check for multiple interface
    // pairs

    if
    (
        fsiProperties_.found("solidPatch")
     && fsiProperties_.found("fluidPatch")
     && !fsiProperties_.found("solidPatches")
     && !fsiProperties_.found("fluidPatches")
    )
    {
        solidPatchNames_.setSize(1, word(fsiProperties_.lookup("solidPatch")));
        fluidPatchNames_.setSize(1, word(fsiProperties_.lookup("fluidPatch")));
    }
    else if
    (
        !fsiProperties_.found("solidPatch")
     && !fsiProperties_.found("fluidPatch")
     && fsiProperties_.found("solidPatches")
     && fsiProperties_.found("fluidPatches")
    )
    {
        solidPatchNames_ = wordList(fsiProperties_.lookup("solidPatches"));
        fluidPatchNames_ = wordList(fsiProperties_.lookup("fluidPatches"));

        if (solidPatchNames_.size() != fluidPatchNames_.size())
        {
            FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
                << "Defined number of coupled fluid and solid patches "
                << "must be equal!" << nl
                << "Currently, there are " << solidPatchNames_.size()
                << " solid interface patches and " << fluidPatchNames_.size()
                << " patches!" << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Either 'solidPatch' and 'fluidPatch' should be defined OR "
            << "'solidPatches' and 'fluidPatches' but not both or neither!"
            << abort(FatalError);
    }

    solidPatchIndices_.setSize(solidPatchNames_.size(), label(-1));
    fluidPatchIndices_.setSize(fluidPatchNames_.size(), label(-1));

    // Set the number of global poly patches: solid or fluid
    nGlobalPatches_ = fluid().globalPatches().size();

    // loop over all coupled patches
    forAll(solidPatchNames_, interfaceI)
    {
        // Solid patch index
        const polyPatchID solidPatch
        (
            solidPatchNames_[interfaceI],
            solidMesh().boundaryMesh()
        );

        if (!solidPatch.active())
        {
            FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
                << "Solid patch name " << solidPatchNames_[interfaceI]
                << " not found!" << abort(FatalError);
        }

        solidPatchIndices_[interfaceI] = solidPatch.index();
    }

    // Create solid global patches
    solid().makeGlobalPatches(solidPatchNames_);

    // loop over all coupled patches
    forAll(fluidPatchNames_, interfaceI)
    {
        // Fluid patch index
        const polyPatchID fluidPatch
        (
            fluidPatchNames_[interfaceI],
            fluidMesh().boundaryMesh()
        );

        if (!fluidPatch.active())
        {
            FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
                << "Fluid patch name " << fluidPatchNames_[interfaceI]
                << " not found!" << abort(FatalError);
        }

        fluidPatchIndices_[interfaceI] = fluidPatch.index();
    }

    // Create fluid global patches
    fluid().makeGlobalPatches(fluidPatchNames_);

    // Set interface fields list size and initialize residual
    fluidZonesPointsDispls_.setSize(nGlobalPatches_);
    fluidZonesPointsDisplsRef_.setSize(nGlobalPatches_);
    fluidZonesPointsDisplsPrev_.setSize(nGlobalPatches_);
    solidZonesPointsDispls_.setSize(nGlobalPatches_);
    solidZonesPointsDisplsRef_.setSize(nGlobalPatches_);
    interfacesPointsDispls_.setSize(nGlobalPatches_);
    interfacesPointsDisplsPrev_.setSize(nGlobalPatches_);
    residuals_.setSize(nGlobalPatches_);
    residualsPrev_.setSize(nGlobalPatches_);
    maxResidualsNorm_.setSize(nGlobalPatches_);
    maxIntsDisplsNorm_.setSize(nGlobalPatches_);

    forAll(residuals_, interfaceI)
    {
        residuals_[interfaceI] = vectorField
        (
            fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
            vector::zero
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
    currentSolidZonesPointsList_.clear();
    currentSolidZonesPatchesList_.clear();
    ggiInterpolatorsList_.clear();
    accumulatedFluidInterfacesDisplacementsList_.clear();
    minEdgeLengthsList_.clear();
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
    if (currentSolidZonesPointsList_.empty())
    {
        calcCurrentSolidZonesPoints();
    }

    return currentSolidZonesPointsList_;
}


const Foam::PtrList<Foam::standAlonePatch>&
Foam::fluidSolidInterface::currentSolidZonesPatches() const
{
    if (currentSolidZonesPatchesList_.empty())
    {
        calcCurrentSolidZonesPatches();
    }

    return currentSolidZonesPatchesList_;
}


const Foam::List<std::shared_ptr<RBFInterpolation> >&
Foam::fluidSolidInterface::rbfFluidToSolid() const
{
    if (rbfFluidToSolidList_.empty())
    {
        calcRbfFluidToSolidInterpolators();
    }

    return rbfFluidToSolidList_;
}


const Foam::List<std::shared_ptr<RBFInterpolation> >&
Foam::fluidSolidInterface::rbfSolidToFluid() const
{
    if (rbfSolidToFluidList_.empty())
    {
        calcRbfSolidToFluidInterpolators();
    }

    return rbfSolidToFluidList_;
}


const Foam::PtrList<Foam::GGIInterpolation<standAlonePatch, standAlonePatch> >&
Foam::fluidSolidInterface::ggiInterpolators() const
{
    if (ggiInterpolatorsList_.empty())
    {
        calcGgiInterpolators();
    }

    return ggiInterpolatorsList_;
}


void Foam::fluidSolidInterface::setDeltaT(Time& runTime)
{
    // For now, the fluid sets the time-step
    fluid().setDeltaT(runTime);
}


void Foam::fluidSolidInterface::initializeFields()
{
    outerCorr_ = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        fluidZonesPointsDispls_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        fluidZonesPointsDisplsRef_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        fluidZonesPointsDisplsPrev_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        solidZonesPointsDispls_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        solidZonesPointsDisplsRef_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        residualsPrev_[interfaceI] = residuals_[interfaceI];

        residuals_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        maxResidualsNorm_[interfaceI] = scalar(0);

        interfacesPointsDispls_[interfaceI] =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );

        interfacesPointsDisplsPrev_ =
            vectorField
            (
                fluid().globalPatches()[interfaceI].globalPatch().nPoints(),
                vector::zero
            );
    }
}


void Foam::fluidSolidInterface::updateInterpolatorAndGlobalPatches()
{
    if (ggiInterpolatorsList_.empty())
    {
        ggiInterpolators();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
            // Clear current interpolators
            ggiInterpolatorsList_.clear();

            // Re-create global patches
            fluid().clearGlobalPatches();
            solid().clearGlobalPatches();
            fluid().makeGlobalPatches(fluidPatchNames_);
            solid().makeGlobalPatches(solidPatchNames_);

            // Re-create interpolators
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
        nGlobalPatches_, vectorField()
    );

    List<vectorField> fluidPatchesPointsDisplsPrev
    (
        nGlobalPatches_, vectorField()
    );

    scalar maxDelta = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        fluidPatchesPointsDispls[interfaceI] =
            fluid().globalPatches()[interfaceI].globalPointToPatch
            (
                fluidZonesPointsDispls()[interfaceI]
            );

        fluidPatchesPointsDisplsPrev[interfaceI] =
            fluid().globalPatches()[interfaceI].globalPointToPatch
            (
                fluidZonesPointsDisplsPrev()[interfaceI]
            );

        // Patch point normals
        const vectorField& n =
            fluid().mesh().boundaryMesh()
            [
                fluid().globalPatches()[interfaceI].patch().index()
            ].pointNormals();

        // Patch deltaCoeffs
        const scalarField fluidZoneDeltaCoeffs =
            fluid().globalPatches()[interfaceI].patchFaceToGlobal
            (
                fluidMesh().boundary()
                [
                    fluid().globalPatches()[interfaceI].patch().index()
                ].deltaCoeffs()
            );

        // Zone deltaCoeffs at points
        const scalarField fluidZonePointDeltaCoeffs =
            fluid().globalPatches()
            [
                interfaceI
            ].interpolator().faceToPointInterpolate(fluidZoneDeltaCoeffs);

        // Patch deltaCoeffs at points
        const scalarField fluidPatchPointDeltaCoeffs =
            fluid().globalPatches()[interfaceI].globalPointToPatch
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
                        accumulatedFluidInterfacesDisplacements()[interfaceI]
                      + fluidPatchesPointsDispls[interfaceI]
                      - fluidPatchesPointsDisplsPrev[interfaceI]
                    )
                )*fluidPatchPointDeltaCoeffs
            );

        Info<< "Maximal accumulated displacement of interface " << interfaceI
            << ": " << delta << endl;

        if (delta > maxDelta)
        {
            maxDelta = delta;
        }
    }

    Info<< nl << "Maximal accumulated displacement of all interfaces: "
        << maxDelta << endl;

    if (maxDelta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().allPoints();

        forAll(fluid().globalPatches(), interfaceI)
        {
            const labelList& meshPoints =
                fluid().globalPatches()[interfaceI].globalPatch().meshPoints();

            forAll(fluidPatchesPointsDispls[interfaceI], pointI)
            {
                newPoints[meshPoints[pointI]] +=
                    fluidPatchesPointsDispls[interfaceI][pointI]
                  - fluidPatchesPointsDisplsPrev[interfaceI][pointI];
            }

            twoDPointCorrector twoDCorrector(fluidMesh());
            twoDCorrector.correctPoints(newPoints);

            fluidMesh().movePoints(newPoints);

            // Accumulate interface points displacement
            accumulatedFluidInterfacesDisplacements()[interfaceI] +=
                fluidPatchesPointsDispls[interfaceI]
              - fluidPatchesPointsDisplsPrev[interfaceI];
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

            forAll(fluid().globalPatches(), interfaceI)
            {
                fixedValueTetPolyPatchVectorField& motionUFluidPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchIndices()[interfaceI]]
                    );

                tetPolyPatchInterpolation tppi
                (
                    refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
                );

                motionUFluidPatch ==
                    tppi.pointToPointInterpolate
                    (
                        accumulatedFluidInterfacesDisplacements()[interfaceI]
                       /fluid().runTime().deltaT().value()
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

            forAll(fluid().globalPatches(), interfaceI)
            {
                fixedValuePointPatchVectorField& motionUFluidPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchIndices()[interfaceI]]
                    );

                motionUFluidPatch ==
                    (
                        fluidPatchesPointsDispls[interfaceI]
                      - fluidPatchesPointsDisplsPrev[interfaceI]
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

                forAll(fluid().globalPatches(), interfaceI)
                {
                    fixedValuePointPatchVectorField& motionUFluidPatch =
                        refCast<fixedValuePointPatchVectorField>
                        (
                            motionU.boundaryField()
                            [
                                fluidPatchIndices()[interfaceI]
                            ]
                        );

                    motionUFluidPatch ==
                    (
                        fluidPatchesPointsDispls[interfaceI]
                      - fluidPatchesPointsDisplsPrev[interfaceI]
                    )/fluid().runTime().deltaT().value();
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

        forAll(fluid().globalPatches(), interfaceI)
        {
            accumulatedFluidInterfacesDisplacements()[interfaceI] =
                vectorField
                (
                    accumulatedFluidInterfacesDisplacements()
                    [
                        interfaceI
                    ].size(),
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
    // Check if coupling switch needs to be updated
    if (!coupled_)
    {
        updateCoupled();
    }

    // PhilipC: we should get rid if the use of List< tmp< > >
    // Here, I think it owuld make sense to have a loop over interfaces
    // and then perform interpolation for each interface: there is no need to
    // create lists of fields to be interpolated
    // Also, I would rather leave the transferFacesZoneToZone function as it was
    // previously where it takes a pair of patches as arguments, rather than
    // forcing it to take lists.

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
        nGlobalPatches_, vectorField()
    );

    List<vectorField> fluidZonesTotalTraction
    (
        nGlobalPatches_, vectorField()
    );

    forAll(fluid().globalPatches(), interfaceI)
    {
        const vectorField& fluidZoneTraction =
            faceZonesViscousForce[interfaceI]();

        const scalarField& fluidZonePressure =
            faceZonesPressureForce[interfaceI]();

        // Fluid zone face normals
        const vectorField& n =
            fluid().globalPatches()[interfaceI].globalPatch().faceNormals();

        // Fluid zone total traction
        fluidZonesTotalTraction[interfaceI] =
            fluidZoneTraction - fluidZonePressure*n;

        // Solid zone traction is to be interpolated from the fluid zone
        solidZonesTotalTraction[interfaceI] =
            vectorField
            (
                solid().globalPatches()[interfaceI].globalPatch().size(),
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
    forAll(fluid().globalPatches(), interfaceI)
    {
        solidZonesTotalTraction[interfaceI] =
            -solidZonesTotalTraction[interfaceI];
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
        forAll(fluid().globalPatches(), interfaceI)
        {
            if
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()[fluidPatchIndices()[interfaceI]]
                )
            )
            {
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().p().boundaryField()[fluidPatchIndices()[interfaceI]]
                    )
                ).prevPressure() = fluid().patchPressureForce(fluidPatchIndices()[interfaceI]);
            }
        }
    }
    else
    {
        // Set interface pressure for elasticWallPressure
        // boundary condition
        forAll(fluid().globalPatches(), interfaceI)
        {
            if
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()[fluidPatchIndices()[interfaceI]]
                )
            )
            {
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().p().boundaryField()
                        [
                            fluidPatchIndices()[interfaceI]
                        ]
                    )
                ).prevPressure() = 0;
            }
        }
    }

    // Total force at the fluid side of the interface
    {
        forAll(fluid().globalPatches(), interfaceI)
        {
            const vectorField& p =
                fluid().globalPatches()[interfaceI].globalPatch().localPoints();
            const faceList& f =
                fluid().globalPatches()[interfaceI].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);
            vectorField C(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
                C[faceI] = f[faceI].centre(p);
            }

            const vector totalTractionForce =
                sum(fluidZonesTotalTraction[interfaceI]*mag(S));

            Info<< "Total force on interface patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << " (fluid) = " << totalTractionForce << endl;
        }
    }

    // Total force at the solid side of the interface
    {
        forAll(solid().globalPatches(), interfaceI)
        {
            const vectorField& p =
                solid().globalPatches()[interfaceI].globalPatch().localPoints();
            const faceList& f =
                solid().globalPatches()[interfaceI].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);
            vectorField C(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
                C[faceI] = f[faceI].centre(p);
            }

            const vector totalTractionForce =
                sum(solidZonesTotalTraction[interfaceI]*mag(S));

            Info<< "Total force on interface patch "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << " (solid) = " << totalTractionForce << endl;
        }
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    List<scalar> minResidual(nGlobalPatches_, scalar(0));

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
        nGlobalPatches_,
        vectorField()
    );

    List<vectorField> solidZonesPointsTotDisplsAtSolid
    (
        nGlobalPatches_,
        vectorField()
    );

    List<vectorField> solidZonesPointsTotDispls
    (
        nGlobalPatches_,
        vectorField()
    );

    forAll(solid().globalPatches(), interfaceI)
    {
        solidZonesPointsDisplsAtSolid[interfaceI] =
            faceZonesPointDisplacementIncrement[interfaceI]();

        solidZonesPointsTotDisplsAtSolid[interfaceI] =
            faceZonesPointDisplacementOld[interfaceI]();

        solidZonesPointsTotDispls[interfaceI] =
            vectorField
            (
                solidZonesPointsDispls()[interfaceI].size(),
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

    forAll(solid().globalPatches(), interfaceI)
    {
        // Update interface residuals
        residualsPrev()[interfaceI] = residuals()[interfaceI];

        residuals()[interfaceI] =
            solidZonesPointsDispls()[interfaceI]
          - fluidZonesPointsDispls()[interfaceI];

        scalar residualNorm = Foam::sqrt(gSum(magSqr(residuals()[interfaceI])));
        scalar residualNorm_2 = residualNorm;

        if (residualNorm > maxResidualsNorm_[interfaceI])
        {
            maxResidualsNorm_[interfaceI] = residualNorm;
        }

        residualNorm /= maxResidualsNorm_[interfaceI] + SMALL;

        Info<< "Current fsi relative residual norm ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[interfaceI].patch().index()
               ].name()
            << "): " << residualNorm << endl;

        interfacesPointsDisplsPrev_[interfaceI] =
            interfacesPointsDispls_[interfaceI];

        interfacesPointsDispls_[interfaceI] =
            solidZonesPointsDispls()[interfaceI];

        const vectorField intTotDispl =
            interfacesPointsDispls_[interfaceI]
          + solidZonesPointsTotDispls[interfaceI];

        const scalar intTotDisplNorm = Foam::sqrt(gSum(magSqr(intTotDispl)));

        if (intTotDisplNorm > maxIntsDisplsNorm_[interfaceI])
        {
            maxIntsDisplsNorm_[interfaceI] = intTotDisplNorm;
        }

        residualNorm_2 /= maxIntsDisplsNorm_[interfaceI] + SMALL;

        Info<< "Alternative fsi residual ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[interfaceI].patch().index()
               ].name()
            << "): " << residualNorm_2 << endl;

        minResidual[interfaceI] = min(residualNorm_2, residualNorm);
    }

    return max(minResidual);
}


void Foam::fluidSolidInterface::updateMovingWallPressureAcceleration()
{
    const List<tmp<vectorField> > faceZonesAcceleration
    (
        solid().faceZonesAcceleration()
    );

    forAll(fluid().globalPatches(), interfaceI)
    {
        if
        (
            isA<movingWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndices()[interfaceI]]
            )
        )
        {
            Info<< "Setting acceleration at fluid side of the interface" << endl;

            const vectorField& solidZoneAcceleration = faceZonesAcceleration[interfaceI]();

            const vectorField fluidZoneAcceleration =
                ggiInterpolators()[interfaceI].slaveToMaster(solidZoneAcceleration);

            const vectorField fluidPatchAcceleration =
                fluid().globalPatches()[interfaceI].globalFaceToPatch(fluidZoneAcceleration);

            const_cast<movingWallPressureFvPatchScalarField&>
            (
                refCast<const movingWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()
                    [
                        fluid().globalPatches()[interfaceI].patch().index()
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

    forAll(fluid().globalPatches(), interfaceI)
    {
        // Set interface acceleration
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndices()[interfaceI]]
            )
        )
        {
            const vectorField& solidZoneAcceleration = faceZonesAcceleration[interfaceI]();

            const vectorField fluidZoneAcceleration =
                ggiInterpolators()[interfaceI].slaveToMaster(solidZoneAcceleration);

            const vectorField fluidPatchAcceleration =
                fluid().globalPatches()[interfaceI].globalFaceToPatch(fluidZoneAcceleration);

            const_cast<elasticWallPressureFvPatchScalarField&>
            (
                refCast<const elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()
                    [
                        fluid().globalPatches()[interfaceI].patch().index()
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
        forAll(fluid().globalPatches(), interfaceI)
        {
            if (!Pstream::master())
            {
                fluidZonesPointsDispls[interfaceI] = vector::zero;
            }

            // pass to all procs
            reduce(fluidZonesPointsDispls[interfaceI], sumOp<vectorField>());

            const labelList& map =
                fluid().globalPatches()[interfaceI].globalMasterToCurrentProcPointAddr();

            if (!Pstream::master())
            {
                const vectorField fluidZonePointsDisplGlobal = fluidZonesPointsDispls[interfaceI];

                forAll(fluidZonePointsDisplGlobal, globalPointI)
                {
                    const label localPoint = map[globalPointI];

                    fluidZonesPointsDispls[interfaceI][localPoint] =
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
