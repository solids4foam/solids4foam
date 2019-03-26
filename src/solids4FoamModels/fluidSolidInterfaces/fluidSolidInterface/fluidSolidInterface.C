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


void Foam::fluidSolidInterface::calcCurrentSolidZonePoints() const
{
    // Find global face zones
    if (currentSolidZonePointsPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentSolidZonePoints() const"
        )   << "Current solid zone points already exist"
            << abort(FatalError);
    }

    currentSolidZonePointsPtr_ =
        new vectorField(solid().currentFaceZonePoints());
}


void Foam::fluidSolidInterface::calcCurrentSolidZonePatch() const
{
    // Find global face zones
    if (currentSolidZonePatchPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentSolidZonePatch() const"
        )   << "Current solid zone patch already exists"
            << abort(FatalError);
    }

    currentSolidZonePatchPtr_ =
        new standAlonePatch
        (
            solid().globalPatch().globalPatch().localFaces(),
            currentSolidZonePoints()
        );
}


void Foam::fluidSolidInterface::calcRbfFluidToSolidInterpolator() const
{
    if (rbfFluidToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcRbfFluidToSolidInterpolator() const"
        )   << "Fluid to solid interpolator already exists"
            << abort(FatalError);
    }

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    rbfFluidToSolidPtr_ =
        std::shared_ptr<RBFInterpolation>(new RBFInterpolation(rbfFunction));

    const vectorField solidZoneFaceCentres =
        currentSolidZonePatch().faceCentres();
    const vectorField fluidZoneFaceCentres =
        fluid().globalPatch().globalPatch().faceCentres();

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

    rbfFluidToSolidPtr_->compute(fluidX, solidX);

    Info<< "Checking fluid-to-solid interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres =
            fluid().globalPatch().patchFaceToGlobal(fluidPatchFaceCentres);


        vectorField solidPatchFaceCentres =
            vectorField
            (
                solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            );

        matrix fluidX(fluidZoneFaceCentres.size(), 3);
        matrix fluidXsolid(solidPatchFaceCentres.size(), 3);

        forAll(fluidZoneFaceCentres, faceI)
        {
            fluidX(faceI, 0) = fluidZoneFaceCentres[faceI].x();
            fluidX(faceI, 1) = fluidZoneFaceCentres[faceI].y();
            fluidX(faceI, 2) = fluidZoneFaceCentres[faceI].z();
        }

        rbfFluidToSolidPtr_->interpolate(fluidX, fluidXsolid);

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
}


void Foam::fluidSolidInterface::calcRbfSolidToFluidInterpolator() const
{
    if (rbfSolidToFluidPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcRbfSolidToFluidInterpolator() const"
        )   << "Solid to fluid interpolator already exists"
            << abort(FatalError);
    }

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    rbfSolidToFluidPtr_ =
        std::shared_ptr<RBFInterpolation>(new RBFInterpolation(rbfFunction));

    const vectorField solidZonePoints = currentSolidZonePatch().localPoints();
    const vectorField fluidZonePoints =
        fluid().globalPatch().globalPatch().localPoints();

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

    rbfSolidToFluidPtr_->compute(solidX, fluidX);

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

        rbfSolidToFluidPtr_->interpolate(solidPoints, fluidPoints);

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
}


void Foam::fluidSolidInterface::calcGgiInterpolator() const
{
    // Create ggi interpolation
    if (ggiInterpolatorPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcGgiInterpolator() const"
        )   << "Ggi interpolator already exists"
            << abort(FatalError);
    }

    // Remove current solid zone and points so that it will be re-created in the
    // deformed position
    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentSolidZonePointsPtr_);

    Info<< "Create GGI zone-to-zone interpolator" << endl;

    ggiInterpolatorPtr_ =
        new GGIInterpolation<standAlonePatch, standAlonePatch>
        (
            fluid().globalPatch().globalPatch(),
            currentSolidZonePatch(),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            true,           // Patch data is complete on all processors
            SMALL,          // Non-overlapping face tolerances
            SMALL,
            true,           // Rescale weighting factors
            ggiInterpolation::BB_OCTREE
        );


    Info<< "Checking fluid-to-solid face interpolator" << endl;

    {
        const vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        const vectorField fluidZoneFaceCentres =
            fluid().globalPatch().patchFaceToGlobal(fluidPatchFaceCentres);

        const vectorField solidZoneFaceCentres =
            ggiInterpolatorPtr_->masterToSlave
            (
                fluidZoneFaceCentres
            );

        const vectorField solidPatchFaceCentres =
            solid().globalPatch().globalFaceToPatch(solidZoneFaceCentres);

        scalar maxDist = gMax
        (
            mag
            (
                solidPatchFaceCentres
              - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            )
        );

        Info<< "Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }

    Info<< "Checking solid-to-fluid point interpolator" << endl;
    {
        const vectorField solidZonePointsAtFluid =
            ggiInterpolatorPtr_->slaveToMasterPointInterpolate
            (
                currentSolidZonePoints()
            );

        const vectorField fluidZonePoints =
            fluid().globalPatch().globalPatch().localPoints();

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
        << ggiInterpolatorPtr_->uncoveredMasterFaces().size() << endl;

    Info<< "Number of uncovered slave faces: "
        << ggiInterpolatorPtr_ ->uncoveredSlaveFaces().size() << endl;

    ggiInterpolatorPtr_->slavePointDistanceToIntersection();
    ggiInterpolatorPtr_->masterPointDistanceToIntersection();
}


void Foam::fluidSolidInterface::calcFluidToSolidFaceMap() const
{
    if (fluidToSolidFaceMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMap() const")
            << "pointer already set!" << abort(FatalError);
    }

    if
    (
        solid().globalPatch().globalPatch().size()
     != fluid().globalPatch().globalPatch().size()
    )
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMap() const")
            << "Fluid and solid interfaces are not conformal!" << nl
            << "directMap method requires conformal interfaces"
            << abort(FatalError);
    }

    IOobject mapHeader
    (
        "fluidToSolidFaceMap",
        fluid().runTime().timeName(),
        fluidMesh(),
        IOobject::MUST_READ
    );

    if (mapHeader.headerOk())
    {
        // Read map
        Info<< "Reading fluidToSolidFaceMap from disk" << endl;
        fluidToSolidFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "fluidToSolidFaceMap",
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
        fluidToSolidFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "fluidToSolidFaceMap",
                    fluid().runTime().timeName(),
                    fluidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(solid().globalPatch().globalPatch().size(), -1)
            )
        );
        labelList& fluidToSolidMap = fluidToSolidFaceMapPtr_();

        // Perform N^2 search for corresponding faces
    // We will take 0.1% of the minEdgeLength as the exact match tolerance

        const vectorField& fluidCf =
            fluid().globalPatch().globalPatch().faceCentres();
        const vectorField& solidCf =
            solid().globalPatch().globalPatch().faceCentres();

        const scalar tol = 0.001*gMin(minEdgeLength());

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

    if (gMin(fluidToSolidFaceMapPtr_()) == -1)
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMap() const")
            << "Cannot calculate the map between interfaces" << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


const Foam::labelList& Foam::fluidSolidInterface::fluidToSolidFaceMap() const
{
    if (fluidToSolidFaceMapPtr_.empty())
    {
        calcFluidToSolidFaceMap();
    }

    return fluidToSolidFaceMapPtr_();
}


void Foam::fluidSolidInterface::calcSolidToFluidFaceMap() const
{
    if (solidToFluidFaceMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcSolidToFluidFaceMap() const")
            << "pointer already set!" << abort(FatalError);
    }

    if
    (
        solid().globalPatch().globalPatch().size()
     != fluid().globalPatch().globalPatch().size()
    )
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMap() const")
            << "Fluid and solid interfaces are not conformal!" << nl
            << "directMap method requires conformal interfaces"
            << abort(FatalError);
    }

    IOobject mapHeader
    (
        "solidToFluidFaceMap",
        solid().runTime().timeName(),
        solidMesh(),
        IOobject::MUST_READ
    );

    if (mapHeader.headerOk())
    {
        // Read map
        Info<< "Reading solidToFluidFaceMap from disk" << endl;
        solidToFluidFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "solidToFluidFaceMap",
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
        solidToFluidFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "solidToFluidFaceMap",
                    solid().runTime().timeName(),
                    solidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(fluid().globalPatch().globalPatch().size(), -1)
            )
        );
        labelList& solidToFluidMap = solidToFluidFaceMapPtr_();

        // Perform N^2 search for corresponding faces
        // We will take 0.1% of the minEdgeLength as the exact match tolerance

        const vectorField& fluidCf =
            fluid().globalPatch().globalPatch().faceCentres();
        const vectorField& solidCf =
            solid().globalPatch().globalPatch().faceCentres();

        const scalar tol = 0.001*gMin(minEdgeLength());

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

    if (gMin(solidToFluidFaceMapPtr_()) == -1)
    {
        FatalErrorIn(type() + "::calcSolidToFluidFaceMap() const")
            << "Cannot calculate the map between interfaces" << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


const Foam::labelList& Foam::fluidSolidInterface::solidToFluidFaceMap() const
{
    if (solidToFluidFaceMapPtr_.empty())
    {
        calcSolidToFluidFaceMap();
    }

    return solidToFluidFaceMapPtr_();
}


void Foam::fluidSolidInterface::calcFluidToSolidPointMap() const
{
    if (fluidToSolidPointMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcFluidToSolidPointMap() const")
            << "pointer already set!" << abort(FatalError);
    }

    if
    (
        solid().globalPatch().globalPatch().size()
     != fluid().globalPatch().globalPatch().size()
    )
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMap() const")
            << "Fluid and solid interfaces are not conformal!" << nl
            << "directMap method requires conformal interfaces"
            << abort(FatalError);
    }

    IOobject mapHeader
    (
        "fluidToSolidPointMap",
        fluid().runTime().timeName(),
        fluidMesh(),
        IOobject::MUST_READ
    );

    if (mapHeader.headerOk())
    {
        // Read map
        Info<< "Reading fluidToSolidPointMap from disk" << endl;
        fluidToSolidPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "fluidToSolidPointMap",
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
        fluidToSolidPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "fluidToSolidPointMap",
                    fluid().runTime().timeName(),
                    fluidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(solid().globalPatch().globalPatch().nPoints(), -1)
            )
        );
        labelList& fluidToSolidMap = fluidToSolidPointMapPtr_();

        // Perform N^2 search for corresponding points
        // We will take 0.1% of the minEdgeLength as the exact match tolerance

        const vectorField& fluidLP =
            fluid().globalPatch().globalPatch().localPoints();
        const vectorField& solidLP =
            solid().globalPatch().globalPatch().localPoints();

        const scalar tol = 0.001*gMin(minEdgeLength());

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

    if (gMin(fluidToSolidPointMapPtr_()) == -1)
    {
        FatalErrorIn(type() + "::calcFluidToSolidPointMap() const")
            << "Cannot calculate the map between interfaces" << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


const Foam::labelList& Foam::fluidSolidInterface::fluidToSolidPointMap() const
{
    if (fluidToSolidPointMapPtr_.empty())
    {
        calcFluidToSolidPointMap();
    }

    return fluidToSolidPointMapPtr_();
}


void Foam::fluidSolidInterface::calcSolidToFluidPointMap() const
{
    if (solidToFluidPointMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcSolidToFluidPointMap() const")
            << "pointer already set!" << abort(FatalError);
    }

    if
    (
        solid().globalPatch().globalPatch().size()
     != fluid().globalPatch().globalPatch().size()
    )
    {
        FatalErrorIn(type() + "::calcFluidToSolidFaceMap() const")
            << "Fluid and solid interfaces are not conformal!" << nl
            << "directMap method requires conformal interfaces"
            << abort(FatalError);
    }

    IOobject mapHeader
    (
        "solidToFluidPointMap",
        solid().runTime().timeName(),
        solidMesh(),
        IOobject::MUST_READ
    );

    if (mapHeader.headerOk())
    {
        // Read map
        Info<< "Reading solidToFluidPointMap from disk" << endl;
        solidToFluidPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "solidToFluidPointMap",
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
        solidToFluidPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    "solidToFluidPointMap",
                    solid().runTime().timeName(),
                    solidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(fluid().globalPatch().globalPatch().nPoints(), -1)
            )
        );
        labelList& solidToFluidMap = solidToFluidPointMapPtr_();

        // Perform N^2 search for corresponding points
        // We will take 0.1% of the minEdgeLength as the exact match tolerance

        const vectorField& fluidLP =
            fluid().globalPatch().globalPatch().localPoints();
        const vectorField& solidLP =
            solid().globalPatch().globalPatch().localPoints();

        const scalar tol = 0.001*gMin(minEdgeLength());

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

    if (gMin(solidToFluidPointMapPtr_()) == -1)
    {
        FatalErrorIn(type() + "::calcSolidToFluidPointMap() const")
            << "Cannot calculate the map between interpoints" << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


const Foam::labelList& Foam::fluidSolidInterface::solidToFluidPointMap() const
{
    if (solidToFluidPointMapPtr_.empty())
    {
        calcSolidToFluidPointMap();
    }

    return solidToFluidPointMapPtr_();

}


void Foam::fluidSolidInterface::
calcAccumulatedFluidInterfaceDisplacement() const
{
    // Read accumulated displacement
    if (accumulatedFluidInterfaceDisplacementPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcAccumulatedFluidInterfaceDisplacement() const"
        )   << "Accumulated displacement field already exists"
            << abort(FatalError);
    }

    // Accumulated fluid interface displacement
    IOobject accumulatedFluidInterfaceDisplacementHeader
    (
        "accumulatedFluidInterfaceDisplacement",
        fluid().runTime().timeName(),
        fluidMesh(),
        IOobject::MUST_READ
    );

    if (accumulatedFluidInterfaceDisplacementHeader.headerOk())
    {
        Pout << "Reading accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    fluid().runTime().timeName(),
                    fluidMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        Pout<< "Creating accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    fluid().runTime().timeName(),
                    fluidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vectorField
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
                    vector::zero
                )
            );
    }
}


void Foam::fluidSolidInterface::calcMinEdgeLength() const
{
    // Read accumulated displacement
    if (minEdgeLengthPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcMinEdgeLength() const"
        )
            << "Minimal edge lengths already exist"
                << abort(FatalError);
    }

    minEdgeLengthPtr_ =
        new scalarField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            0
        );
    scalarField& minEdgeLength = *minEdgeLengthPtr_;


    const edgeList& edges =
        fluid().globalPatch().globalPatch().edges();

    const vectorField& points =
        fluid().globalPatch().globalPatch().localPoints();

    const labelListList& pointEdges =
        fluid().globalPatch().globalPatch().pointEdges();

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

//     Pout << "Min edge length: " << min(minEdgeLength) << endl;
//     Pout << "gMin edge length: " << gMin(minEdgeLength) << endl;
}


Foam::vectorIOField&
Foam::fluidSolidInterface::accumulatedFluidInterfaceDisplacement()
{
    if (!accumulatedFluidInterfaceDisplacementPtr_)
    {
        calcAccumulatedFluidInterfaceDisplacement();
    }

    return *accumulatedFluidInterfaceDisplacementPtr_;
}


const Foam::scalarField& Foam::fluidSolidInterface::minEdgeLength() const
{
    if (!minEdgeLengthPtr_)
    {
        calcMinEdgeLength();
    }

    return *minEdgeLengthPtr_;
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
    solidPatchIndex_(-1),
    fluidPatchIndex_(-1),
    currentSolidZonePointsPtr_(NULL),
    currentSolidZonePatchPtr_(NULL),
    rbfFluidToSolidPtr_(NULL),
    rbfSolidToFluidPtr_(NULL),
    ggiInterpolatorPtr_(NULL),
    fluidToSolidFaceMapPtr_(),
    solidToFluidFaceMapPtr_(),
    fluidToSolidPointMapPtr_(),
    solidToFluidPointMapPtr_(),
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
    fluidZonePointsDispl_(),
    fluidZonePointsDisplRef_(),
    fluidZonePointsDisplPrev_(),
    solidZonePointsDispl_(),
    solidZonePointsDisplRef_(),
    interfacePointsDispl_(),
    interfacePointsDisplPrev_(),
    residual_(),
    residualPrev_(),
    maxResidualNorm_(0),
    maxIntDisplNorm_(0),
    outerCorr_(0),
    interpolatorUpdateFrequency_
    (
        fsiProperties_.lookupOrDefault<int>("interpolatorUpdateFrequency", 0)
    ),
    accumulatedFluidInterfaceDisplacementPtr_(NULL),
    minEdgeLengthPtr_(NULL),
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

    // Solid patch index
    const word solidPatchName(fsiProperties_.lookup("solidPatch"));
    const polyPatchID solidPatch
    (
        solidPatchName,
        solidMesh().boundaryMesh()
    );

    if (!solidPatch.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Solid patch name " << solidPatchName << " not found."
            << abort(FatalError);
    }

    solidPatchIndex_ = solidPatch.index();

    // Create solid global patch
    solid().makeGlobalPatch(solidPatchName);

    // Fluid patch index
    const word fluidPatchName(fsiProperties_.lookup("fluidPatch"));
    const polyPatchID fluidPatch
    (
        fluidPatchName,
        fluidMesh().boundaryMesh()
    );

    if (!fluidPatch.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Fluid patch name " << fluidPatchName << " not found."
            << abort(FatalError);
    }

    fluidPatchIndex_ = fluidPatch.index();

    // Create fluid global patch
    fluid().makeGlobalPatch(fluidPatchName);

    // Initialize residual
    residual_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(), vector::zero
        );

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
        fluidToSolidFaceMap();
        solidToFluidFaceMap();
        fluidToSolidPointMap();
        solidToFluidPointMap();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::~fluidSolidInterface()
{
    deleteDemandDrivenData(currentSolidZonePointsPtr_);
    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(ggiInterpolatorPtr_);
    deleteDemandDrivenData(accumulatedFluidInterfaceDisplacementPtr_);
    deleteDemandDrivenData(minEdgeLengthPtr_);
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


const Foam::vectorField&
Foam::fluidSolidInterface::currentSolidZonePoints() const
{
    if (!currentSolidZonePointsPtr_)
    {
        calcCurrentSolidZonePoints();
    }

    return *currentSolidZonePointsPtr_;
}


const Foam::standAlonePatch&
Foam::fluidSolidInterface::currentSolidZonePatch() const
{
    if (!currentSolidZonePatchPtr_)
    {
        calcCurrentSolidZonePatch();
    }

    return *currentSolidZonePatchPtr_;
}


const std::shared_ptr<RBFInterpolation>&
Foam::fluidSolidInterface::rbfFluidToSolid() const
{
    if (!rbfFluidToSolidPtr_)
    {
        calcRbfFluidToSolidInterpolator();
    }

    return rbfFluidToSolidPtr_;
}


const Foam::GGIInterpolation<standAlonePatch, standAlonePatch>&
Foam::fluidSolidInterface::ggiInterpolator() const
{
    if (!ggiInterpolatorPtr_)
    {
        calcGgiInterpolator();
    }

    return *ggiInterpolatorPtr_;
}


const std::shared_ptr<RBFInterpolation>&
Foam::fluidSolidInterface::rbfSolidToFluid() const
{
    if (!rbfSolidToFluidPtr_)
    {
        calcRbfSolidToFluidInterpolator();
    }

    return rbfSolidToFluidPtr_;
}


void Foam::fluidSolidInterface::setDeltaT(Time& runTime)
{
    // For now, the fluid sets the time-step
    fluid().setDeltaT(runTime);
}


void Foam::fluidSolidInterface::initializeFields()
{
    fluidZonePointsDispl_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplRef_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplPrev_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    solidZonePointsDispl_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    solidZonePointsDisplRef_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    residualPrev_ = residual_;

    residual_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    maxResidualNorm_ = 0;

    outerCorr_ = 0;

    interfacePointsDispl_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );

    interfacePointsDisplPrev_ =
        vectorField
        (
            fluid().globalPatch().globalPatch().nPoints(),
            vector::zero
        );
}


void Foam::fluidSolidInterface::updateInterpolatorAndGlobalPatches()
{
    if (!ggiInterpolatorPtr_)
    {
        ggiInterpolator();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
            deleteDemandDrivenData(ggiInterpolatorPtr_);
            fluid().clearGlobalPatch();
            solid().clearGlobalPatch();
            fluid().makeGlobalPatch
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].name()
            );
            solid().makeGlobalPatch
            (
                solidMesh().boundaryMesh()[solidPatchIndex_].name()
            );
            ggiInterpolator();
        }
    }
}


void Foam::fluidSolidInterface::moveFluidMesh()
{
    // Get fluid patch displacement from fluid zone displacement

    // Take care: these are local patch fields not global patch fields

    const vectorField fluidPatchPointsDispl =
        fluid().globalPatch().globalPointToPatch(fluidZonePointsDispl());

    const vectorField fluidPatchPointsDisplPrev =
        fluid().globalPatch().globalPointToPatch(fluidZonePointsDisplPrev());

    // Patch point normals
    const vectorField& n =
        fluid().mesh().boundaryMesh()
        [
            fluid().globalPatch().patch().index()
        ].pointNormals();

    // Patch deltaCoeffs
    const scalarField fluidZoneDeltaCoeffs =
        fluid().globalPatch().patchFaceToGlobal
        (
            fluidMesh().boundary()
            [
                fluid().globalPatch().patch().index()
            ].deltaCoeffs()
        );

    // Zone deltaCoeffs at points
    const scalarField fluidZonePointDeltaCoeffs =
        fluid().globalPatch().interpolator().faceToPointInterpolate
        (
            fluidZoneDeltaCoeffs
        );

    // Patch deltaCoeffs at points
    const scalarField fluidPatchPointDeltaCoeffs =
        fluid().globalPatch().globalPointToPatch(fluidZonePointDeltaCoeffs);

    const scalar delta =
        gMax
        (
            mag
            (
                n
              & (
                    accumulatedFluidInterfaceDisplacement()
                  + fluidPatchPointsDispl
                  - fluidPatchPointsDisplPrev
                )
            )*fluidPatchPointDeltaCoeffs
        );

    Info<< "Maximal accumulated displacement of interface points: "
        << delta << endl;

    if (delta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().allPoints();

        const labelList& meshPoints =
            fluid().globalPatch().globalPatch().meshPoints();

        forAll(fluidPatchPointsDispl, pointI)
        {
            newPoints[meshPoints[pointI]] +=
                fluidPatchPointsDispl[pointI]
              - fluidPatchPointsDisplPrev[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());
        twoDCorrector.correctPoints(newPoints);

        fluidMesh().movePoints(newPoints);

        // Accumulate interface points displacement
        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;
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

            fixedValueTetPolyPatchVectorField& motionUFluidPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            tetPolyPatchInterpolation tppi
            (
                refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
            );

            motionUFluidPatch ==
                tppi.pointToPointInterpolate
                (
                    accumulatedFluidInterfaceDisplacement()
                   /fluid().runTime().deltaT().value()
                );
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

            fixedValuePointPatchVectorField& motionUFluidPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            motionUFluidPatch ==
                (
                    fluidPatchPointsDispl
                  - fluidPatchPointsDisplPrev
                )/fluid().runTime().deltaT().value();
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

                fixedValuePointPatchVectorField& motionUFluidPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchIndex()]
                    );

                motionUFluidPatch ==
                (
                    fluidPatchPointsDispl
                  - fluidPatchPointsDisplPrev
                )/fluid().runTime().deltaT().value();

                // FatalErrorIn("fluidSolidInterface::moveFluidMesh()")
                //   << "subset fvMotionSolver"
                //   << abort(FatalError);
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

        accumulatedFluidInterfaceDisplacement() =
            vectorField
            (
                accumulatedFluidInterfaceDisplacement().size(),
                vector::zero
            );
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
    Info<< "Setting traction on solid patch" << endl;

    const vectorField fluidZoneTraction = fluid().faceZoneViscousForce();

    const scalarField fluidZonePressure = fluid().faceZonePressureForce();

    // Fluid zone face normals
    const vectorField& n = fluid().globalPatch().globalPatch().faceNormals();

    // Fluid zone total traction
    const vectorField fluidZoneTotalTraction =
        fluidZoneTraction - fluidZonePressure*n;

    // Solid zone traction is interpolated from the fluid zone

    vectorField solidZoneTotalTraction
    (
        solid().globalPatch().globalPatch().size(),
        vector::zero
    );

    transferFacesZoneToZone
    (
        "fluid",                             // from region name
        "solid",                             // to region name
        fluid().globalPatch().globalPatch(), // from zone
        solid().globalPatch().globalPatch(), // to zone
        fluidZoneTotalTraction,              // from field
        solidZoneTotalTraction               // to field
    );

    // Flip traction sign after transferring from fluid to solid
    solidZoneTotalTraction = -solidZoneTotalTraction;

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
            solidPatchIndex(),
            solidZoneTotalTraction
        );

        // Set interface pressure for elasticWallPressure
        // boundary condition
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndex_]
            )
        )
        {
            const_cast<elasticWallPressureFvPatchScalarField&>
            (
                refCast<const elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()[fluidPatchIndex_]
                )
            ).prevPressure() = fluid().patchPressureForce(fluidPatchIndex_);
        }
    }
    else
    {
        // Set interface pressure for elasticWallPressure
        // boundary condition
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndex_]
            )
        )
        {
            const_cast<elasticWallPressureFvPatchScalarField&>
            (
                refCast<const elasticWallPressureFvPatchScalarField>
                (
                    fluid().p().boundaryField()[fluidPatchIndex_]
                )
            ).prevPressure() = 0;
        }
    }

    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluid().globalPatch().globalPatch().localPoints();
        const faceList& f =
            fluid().globalPatch().globalPatch().localFaces();

        vectorField S(f.size(), vector::zero);
        vectorField C(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
            C[faceI] = f[faceI].centre(p);
        }

        const vector totalTractionForce = sum(fluidZoneTotalTraction*mag(S));

        Info<< "Total force (fluid) = " << totalTractionForce << endl;
    }

    // Total force at the solid side of the interface
    {
        const vectorField& p =
            solid().globalPatch().globalPatch().localPoints();
        const faceList& f =
            solid().globalPatch().globalPatch().localFaces();

        vectorField S(f.size(), vector::zero);
        vectorField C(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
            C[faceI] = f[faceI].centre(p);
        }

        const vector totalTractionForce = sum(solidZoneTotalTraction*mag(S));

        Info<< "Total force (solid) = " << totalTractionForce << endl;
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    const vectorField solidZonePointsDisplAtSolid =
        solid().faceZonePointDisplacementIncrement();

    const vectorField solidZonePointsTotDisplAtSolid =
        solid().faceZonePointDisplacementOld();

    vectorField solidZonePointsTotDispl
    (
        solidZonePointsDispl().size(),
        vector::zero
    );

    // Transfer displacement fields from the solid to the fluid

    transferPointsZoneToZone
    (
        "solid",                             // from region name
        "fluid",                             // to region name
        solid().globalPatch().globalPatch(), // from zone
        fluid().globalPatch().globalPatch(), // to zone
        solidZonePointsDisplAtSolid,         // from field
        solidZonePointsDispl()               // to field
    );

    transferPointsZoneToZone
    (
        "solid",                             // from region name
        "fluid",                             // to region name
        solid().globalPatch().globalPatch(), // from zone
        fluid().globalPatch().globalPatch(), // to zone
        solidZonePointsTotDisplAtSolid,      // from field
        solidZonePointsTotDispl              // to field
    );


    // Update interface residuals

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

    scalar residualNorm = Foam::sqrt(gSum(magSqr(residual())));
    scalar residualNorm_2 = residualNorm;

    if (residualNorm > maxResidualNorm_)
    {
        maxResidualNorm_ = residualNorm;
    }

    residualNorm /= maxResidualNorm_ + SMALL;

    Info<< "Current fsi relative residual norm: " << residualNorm << endl;

    interfacePointsDisplPrev_ = interfacePointsDispl_;

    interfacePointsDispl_ = solidZonePointsDispl();

    const vectorField intTotDispl =
        interfacePointsDispl_ + solidZonePointsTotDispl;
    const scalar intTotDisplNorm = ::sqrt(gSum(magSqr(intTotDispl)));
    if (intTotDisplNorm > maxIntDisplNorm_)
    {
        maxIntDisplNorm_ = intTotDisplNorm;
    }

    residualNorm_2 /= maxIntDisplNorm_ + SMALL;

    Info<< "Alternative fsi residual: " << residualNorm_2 << endl;

    return min(residualNorm_2, residualNorm);
}


void Foam::fluidSolidInterface::updateMovingWallPressureAcceleration()
{
    if
    (
        isA<movingWallPressureFvPatchScalarField>
        (
            fluid().p().boundaryField()[fluidPatchIndex()]
        )
    )
    {
        Info<< "Setting acceleration at fluid side of the interface" << endl;

        const vectorField solidZoneAcceleration =
            solid().faceZoneAcceleration();

        const vectorField fluidZoneAcceleration =
            ggiInterpolator().slaveToMaster(solidZoneAcceleration);

        const vectorField fluidPatchAcceleration =
            fluid().globalPatch().globalFaceToPatch(fluidZoneAcceleration);

        const_cast<movingWallPressureFvPatchScalarField&>
        (
            refCast<const movingWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()
                [
                    fluid().globalPatch().patch().index()
                ]
            )
        ).prevAcceleration() = fluidPatchAcceleration;
    }
}


void Foam::fluidSolidInterface::updateElasticWallPressureAcceleration()
{
    // Set interface acceleration
    if
    (
        isA<elasticWallPressureFvPatchScalarField>
        (
            fluid().p().boundaryField()[fluidPatchIndex()]
        )
    )
    {
        const vectorField solidZoneAcceleration =
            solid().faceZoneAcceleration();

        const vectorField fluidZoneAcceleration =
            ggiInterpolator().slaveToMaster(solidZoneAcceleration);

        const vectorField fluidPatchAcceleration =
            fluid().globalPatch().globalFaceToPatch(fluidZoneAcceleration);

        const_cast<elasticWallPressureFvPatchScalarField&>
        (
            refCast<const elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()
                [
                    fluid().globalPatch().patch().index()
                ]
            )
        ).prevAcceleration() = fluidPatchAcceleration;
    }
}


void Foam::fluidSolidInterface::syncFluidZonePointsDispl
(
    vectorField& fluidZonePointsDispl
)
{
    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if (!Pstream::master())
        {
            fluidZonePointsDispl = vector::zero;
        }

        // pass to all procs
        reduce(fluidZonePointsDispl, sumOp<vectorField>());

        const labelList& map =
            fluid().globalPatch().globalMasterToCurrentProcPointAddr();

        if (!Pstream::master())
        {
            const vectorField fluidZonePointsDisplGlobal = fluidZonePointsDispl;

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                const label localPoint = map[globalPointI];

                fluidZonePointsDispl[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
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
