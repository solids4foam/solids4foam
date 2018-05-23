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
#include "ggiInterpolation.H"
#include "TPSFunction.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSolidInterface, 0);
    defineRunTimeSelectionTable(fluidSolidInterface, dictionary);
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
        new vectorField(solid().currentFaceZonePoints(solidZoneIndex()));
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
        new PrimitivePatch<face, List, const pointField&>
        (
            solidMesh().faceZones()[solidZoneIndex_]().localFaces(),
            currentSolidZonePoints()
        );
}


void Foam::fluidSolidInterface::calcFluidToSolidInterpolator() const
{
    if (fluidToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcFluidToSolidInterpolator() const"
        )   << "Fluid to solid interpolator already exists"
            << abort(FatalError);
    }

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    fluidToSolidPtr_ =
        std::shared_ptr<RBFInterpolation>(new RBFInterpolation(rbfFunction));

    vectorField solidZoneFaceCentres = currentSolidZonePatch().faceCentres();
    vectorField fluidZoneFaceCentres =
        fluidMesh().faceZones()[fluidZoneIndex_]().faceCentres();

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

    fluidToSolidPtr_->compute(fluidX, solidX);

    Info<< "Checking fluid-to-solid interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

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

        fluidToSolidPtr_->interpolate(fluidX, fluidXsolid);

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

    // Create copy of solid face zone primitive patch in current configuration

    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentSolidZonePointsPtr_);

    //currentSolidZonePatch().movePoints(currentSolidZonePoints());

    Info<< "Create GGI zone-to-zone interpolator" << endl;

    ggiInterpolatorPtr_ =
        new ggiZoneInterpolation
        (
            fluidMesh().faceZones()[fluidZoneIndex_](),
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
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] = fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres =
            ggiInterpolatorPtr_->masterToSlave
            (
                fluidZoneFaceCentres
            );

        vectorField solidPatchFaceCentres
        (
            solidMesh().boundaryMesh()[solidPatchIndex_].size(),
            vector::zero
        );

        const label solidPatchStart =
            solidMesh().boundaryMesh()[solidPatchIndex_].start();

        forAll(solidPatchFaceCentres, i)
        {
            solidPatchFaceCentres[i] =
                solidZoneFaceCentres
                [
                    solidMesh().faceZones()[solidZoneIndex_]
                   .whichFace(solidPatchStart + i)
                ];
        }

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

    Info<< "Checking solid-to-fluid point interpolator (GGI)" << endl;
    {
        vectorField solidZonePoints_ = currentSolidZonePoints();
        // solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        vectorField solidZonePoints =
            ggiInterpolatorPtr_->slaveToMasterPointInterpolate
            (
                solidZonePoints_
            );

        vectorField fluidZonePoints =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        scalar maxDist = gMax
        (
            mag
            (
                fluidZonePoints
              - solidZonePoints
            )
        );

        Info<< "Solid-to-fluid point interpolation error (GGI): " << maxDist
            << endl;
    }

    Info<< "Number of uncovered master faces: "
        << ggiInterpolatorPtr_->uncoveredMasterFaces().size() << endl;

    Info<< "Number of uncovered slave faces: "
        << ggiInterpolatorPtr_ ->uncoveredSlaveFaces().size() << endl;

    ggiInterpolatorPtr_->slavePointDistanceToIntersection();
    ggiInterpolatorPtr_->masterPointDistanceToIntersection();
}


void Foam::fluidSolidInterface::calcSolidToFluidInterpolator() const
{
    if (solidToFluidPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcSolidToFluidInterpolator() const"
        )   << "Solid to fluid interpolator already exists"
            << abort(FatalError);
    }

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    solidToFluidPtr_ =
        std::shared_ptr<RBFInterpolation>(new RBFInterpolation(rbfFunction));

    vectorField solidZonePoints = currentSolidZonePatch().localPoints();
    vectorField fluidZonePoints =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

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

    solidToFluidPtr_->compute(solidX, fluidX);

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

        solidToFluidPtr_->interpolate(solidPoints, fluidPoints);

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
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            0
        );
    scalarField& minEdgeLength = *minEdgeLengthPtr_;


    const edgeList& edges =
        fluidMesh().faceZones()[fluidZoneIndex_]().edges();

    const vectorField& points =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

    const labelListList& pointEdges =
        fluidMesh().faceZones()[fluidZoneIndex_]().pointEdges();

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
    solidZoneIndex_(-1),
    fluidPatchIndex_(-1),
    fluidZoneIndex_(-1),
    currentSolidZonePointsPtr_(NULL),
    currentSolidZonePatchPtr_(NULL),
    fluidToSolidPtr_(NULL),
    ggiInterpolatorPtr_(NULL),
    solidToFluidPtr_(NULL),
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
    predictor_(lookupOrDefault<Switch>("predictor", false)),
    rbfInterpolation_(lookupOrDefault<Switch>("rbfInterpolation", false)),
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
    solidZonePressure_(),
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
    minEdgeLengthPtr_(NULL)
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


    // Solid face zone index

    const word solidZoneName(fsiProperties_.lookup("solidZone"));

    const faceZoneID solidZone
    (
        solidZoneName,
        solidMesh().faceZones()
    );

    if (!solidZone.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Solid face zone name " << solidZoneName
            << " not found.  Please check your face zone definition."
            << abort(FatalError);
    }

    solidZoneIndex_ = solidZone.index();


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


    // Fluid face zone index

    const word fluidZoneName(fsiProperties_.lookup("fluidZone"));

    const faceZoneID fluidZone
    (
        fluidZoneName,
        fluidMesh().faceZones()
    );

    if (!fluidZone.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Fluid face zone name " << fluidZoneName
            << " not found.  Please check your face zone definition."
            << abort(FatalError);
    }

    fluidZoneIndex_ = fluidZone.index();

    // Initialize solid zone pressure field
    solidZonePressure_ =
        scalarField(solidMesh().faceZones()[solidZoneIndex()].size(), 0.0);

    // Initialize residual
    residual_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );
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


const Foam::PrimitivePatch<Foam::face, Foam::List, const Foam::pointField&>&
Foam::fluidSolidInterface::currentSolidZonePatch() const
{
    if (!currentSolidZonePatchPtr_)
    {
        calcCurrentSolidZonePatch();
    }

    return *currentSolidZonePatchPtr_;
}


const std::shared_ptr<RBFInterpolation>&
Foam::fluidSolidInterface::fluidToSolid() const
{
    if (!fluidToSolidPtr_)
    {
        calcFluidToSolidInterpolator();
    }

    return fluidToSolidPtr_;
}


const Foam::ggiZoneInterpolation&
Foam::fluidSolidInterface::ggiInterpolator() const
{
    if (!ggiInterpolatorPtr_)
    {
        calcGgiInterpolator();
    }

    return *ggiInterpolatorPtr_;
}


const std::shared_ptr<RBFInterpolation>&
Foam::fluidSolidInterface::solidToFluid() const
{
    if (!solidToFluidPtr_)
    {
        calcSolidToFluidInterpolator();
    }

    return solidToFluidPtr_;
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
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplPrev_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDispl_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    residualPrev_ = residual_;

    residual_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    maxResidualNorm_ = 0;

    outerCorr_ = 0;

    nOuterCorr_ = fsiProperties_.lookupOrDefault<int>("nOuterCorr", 30);

    outerCorrTolerance_ =
        fsiProperties_.lookupOrDefault<scalar>("outerCorrTolerance", 1e-6);

    coupled_ = fsiProperties_.lookupOrDefault<Switch>("coupled", true);

    interfacePointsDispl_ =
        vectorField
        (
            fluid().mesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    interfacePointsDisplPrev_ =
        vectorField
        (
            fluid().mesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );
}


void Foam::fluidSolidInterface::updateInterpolator()
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
            ggiInterpolator();
        }
    }
}


void Foam::fluidSolidInterface::moveFluidMesh()
{
    // Get fluid patch displacement from fluid zone displacement

    vectorField fluidPatchPointsDispl
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
        vector::zero
    );

    vectorField fluidPatchPointsDisplPrev
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
        vector::zero
    );

    const labelList& fluidPatchMeshPoints =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

    forAll(fluidPatchPointsDispl, pointI)
    {
        label curMeshPointID = fluidPatchMeshPoints[pointI];

        label curFluidZonePointID =
            fluidMesh().faceZones()[fluidZoneIndex()]()
           .whichPoint(curMeshPointID);

        fluidPatchPointsDispl[pointI] =
            fluidZonePointsDispl()[curFluidZonePointID];

        fluidPatchPointsDisplPrev[pointI] =
            fluidZonePointsDisplPrev()[curFluidZonePointID];
    }

    // Move fluid mesh
    const vectorField& n =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].pointNormals();

    primitivePatchInterpolation patchInterpolator
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()]
    );

    scalarField pointDeltaCoeffs =
        patchInterpolator.faceToPointInterpolate
        (
            fluidMesh().boundary()[fluidPatchIndex()].deltaCoeffs()
        );

    scalar delta =
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
            )*pointDeltaCoeffs
        );

    Info<< "Maximal accumulated displacement of interface points: "
        << delta << endl;

    if (delta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().allPoints();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll (fluidPatchPointsDispl, pointI)
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
        // pointField newPoints = fluidMesh().allPoints();

        // const labelList& meshPoints =
        //     fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        // forAll (accumulatedFluidInterfaceDisplacement(), pointI)
        // {
        //     newPoints[meshPoints[pointI]] -=
        //         accumulatedFluidInterfaceDisplacement()[pointI];
        // }

        // twoDPointCorrector twoDCorrector(fluidMesh());

        // twoDCorrector.correctPoints(newPoints);

        // fluidMesh().movePoints(newPoints);

        // accumulatedFluidInterfaceDisplacement() +=
        //     fluidPatchPointsDispl
        //   - fluidPatchPointsDisplPrev;

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
    {
        vectorField newPoints = fluidMesh().allPoints();

        const labelList& fluidZoneMeshPoints =
            fluidMesh().faceZones()[fluidZoneIndex()]().meshPoints();

        forAll(fluidZonePointsDispl(), pointI)
        {
            if (fluidZoneMeshPoints[pointI] >= fluidMesh().nPoints())
            {
                newPoints[fluidZoneMeshPoints[pointI]] +=
                    fluidZonePointsDispl()[pointI]
                  - fluidZonePointsDisplPrev()[pointI];
            }
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh().movePoints(newPoints);
    }
}


void Foam::fluidSolidInterface::updateForce()
{
    Info<< "Setting traction on solid patch" << endl;

    const vectorField fluidZoneTraction =
        fluid().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        );

    const scalarField fluidZonePressure =
        fluid().faceZonePressureForce(fluidZoneIndex(), fluidPatchIndex());


    // Fluid zone face normals

    vectorField p = fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();
    const labelList& mp =
        fluidMesh().faceZones()[fluidZoneIndex_]().meshPoints();
    const vectorField& allPoints = fluidMesh().allPoints();
    forAll(mp, pI)
    {
        p[pI] = allPoints[mp[pI]];
    }

    const faceList& f =
        fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    // Fluid zone total traction
    const vectorField fluidZoneTotalTraction =
        fluidZoneTraction - fluidZonePressure*n;

    // Solid zone traction is interpolated from the fluid zone

    vectorField solidZoneTotalTraction
    (
        solidMesh().faceZones()[solidZoneIndex()].size(),
        vector::zero
    );

    if (rbfInterpolation_)
    {
        Info << "... using RBF interpolation" << endl;

        matrix fluidForce(fluidZoneTotalTraction.size(), 3);
        matrix solidForce(solidZoneTotalTraction.size(), 3);

        forAll(fluidZoneTotalTraction, faceI)
        {
            fluidForce(faceI, 0) = fluidZoneTotalTraction[faceI].x();
            fluidForce(faceI, 1) = fluidZoneTotalTraction[faceI].y();
            fluidForce(faceI, 2) = fluidZoneTotalTraction[faceI].z();
        }

        fluidToSolid()->interpolate(fluidForce, solidForce);

        forAll(solidZoneTotalTraction, faceI)
        {
            solidZoneTotalTraction[faceI].x() = -solidForce(faceI, 0);
            solidZoneTotalTraction[faceI].y() = -solidForce(faceI, 1);
            solidZoneTotalTraction[faceI].z() = -solidForce(faceI, 2);
        }
    }
    else
    {
        solidZoneTotalTraction =
            ggiInterpolator().masterToSlave
            (
              - fluidZoneTotalTraction
            );
    }

    solidZonePressure_ =
        ggiInterpolator().masterToSlave
        (
            fluidZonePressure
        );

    // Debugging: print traction fields
    // if (false)
    // {
    //     volVectorField fluidTraction
    //     (
    //         IOobject
    //         (
    //             "fluidTraction",
    //             runTime().timeName(),
    //             fluidMesh(),
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         fluidMesh(),
    //         dimensionedVector("0", dimPressure, vector::zero)
    //     );
    //     fluidTraction.boundaryField()[fluidPatchIndex_] =
    //       fluidZoneTotalTraction;
    //     fluidTraction.write();

    //     volVectorField solidTraction
    //     (
    //         IOobject
    //         (
    //             "solidTraction",
    //             runTime().timeName(),
    //             solidMesh(),
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         solidMesh(),
    //         dimensionedVector("0", dimPressure, vector::zero)
    //     );
    //     solidTraction.boundaryField()[solidPatchIndex_] =
    //         solidZoneTotalTraction;
    //     solidTraction.write();
    // }

    if (!coupled_)
    {
        updateCoupled();
    }

    if (coupled())
    {
        solid().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
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
    //if (true)
    {
        vectorField p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const labelList& mp =
            fluidMesh().faceZones()[fluidZoneIndex_]().meshPoints();
        const vectorField& allPoints = fluidMesh().allPoints();
        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

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

    // Totla force at the solid side of the interface
    //if (true)
    {
        vectorField p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const labelList& mp =
            solidMesh().faceZones()[solidZoneIndex_]().meshPoints();
        const vectorField& allPoints = solidMesh().allPoints();
        forAll(mp, pI)
        {
            p[pI] = allPoints[mp[pI]];
        }

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

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
        solid().faceZonePointDisplacementIncrement(solidZoneIndex());

    const vectorField solidZonePointsTotDisplAtSolid =
        solid().faceZonePointDisplacement(solidZoneIndex());

    vectorField solidZonePointsTotDispl
    (
        solidZonePointsDispl().size(),
        vector::zero
    );

    if (rbfInterpolation_)
    {
        Info<< "Displacement interpolation using RBF interpolation" << endl;

        matrix fluidDispl(solidZonePointsDispl().size(), 3);
        matrix solidDispl(solidZonePointsDisplAtSolid.size(), 3);

        forAll(solidZonePointsDisplAtSolid, pointI)
        {
            solidDispl(pointI, 0) = solidZonePointsDisplAtSolid[pointI].x();
            solidDispl(pointI, 1) = solidZonePointsDisplAtSolid[pointI].y();
            solidDispl(pointI, 2) = solidZonePointsDisplAtSolid[pointI].z();
        }

        solidToFluid()->interpolate(solidDispl, fluidDispl);

        forAll(solidZonePointsDispl(), pointI)
        {
            solidZonePointsDispl()[pointI].x() = fluidDispl(pointI, 0);
            solidZonePointsDispl()[pointI].y() = fluidDispl(pointI, 1);
            solidZonePointsDispl()[pointI].z() = fluidDispl(pointI, 2);
        }

        // Total displacement
        forAll(solidZonePointsTotDisplAtSolid, pointI)
        {
            solidDispl(pointI, 0) = solidZonePointsTotDisplAtSolid[pointI].x();
            solidDispl(pointI, 1) = solidZonePointsTotDisplAtSolid[pointI].y();
            solidDispl(pointI, 2) = solidZonePointsTotDisplAtSolid[pointI].z();
        }

        solidToFluid()->interpolate(solidDispl, fluidDispl);

        forAll(solidZonePointsTotDispl, pointI)
        {
            solidZonePointsTotDispl[pointI].x() = fluidDispl(pointI, 0);
            solidZonePointsTotDispl[pointI].y() = fluidDispl(pointI, 1);
            solidZonePointsTotDispl[pointI].z() = fluidDispl(pointI, 2);
        }
    }
    else
    {
        solidZonePointsDispl() =
            ggiInterpolator().slaveToMasterPointInterpolate
            (
                solidZonePointsDisplAtSolid
            );

        solidZonePointsTotDispl =
            ggiInterpolator().slaveToMasterPointInterpolate
            (
                solidZonePointsTotDisplAtSolid
            );
    }

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

    scalar residualNorm = ::sqrt(gSum(magSqr(residual())));
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
            solid().faceZoneAcceleration(solidZoneIndex(), solidPatchIndex());

        const vectorField fluidZoneAcceleration =
            ggiInterpolator().slaveToMaster(solidZoneAcceleration);

        vectorField fluidPatchAcceleration
        (
            fluidMesh().boundary()[fluidPatchIndex()].size(),
            vector::zero
        );

        const label patchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].start();

        forAll(fluidPatchAcceleration, I)
        {
            fluidPatchAcceleration[I] =
                fluidZoneAcceleration
                [
                    fluidMesh().faceZones()[fluidZoneIndex()].whichFace
                    (
                        patchStart + I
                    )
                ];
        }

        const_cast<movingWallPressureFvPatchScalarField&>
        (
            refCast<const movingWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndex()]
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
            solid().faceZoneAcceleration(solidZoneIndex(), solidPatchIndex());

        const vectorField fluidZoneAcceleration =
            ggiInterpolator().slaveToMaster(solidZoneAcceleration);

        vectorField fluidPatchAcceleration
        (
            fluidMesh().boundary()[fluidPatchIndex()].size(),
            vector::zero
        );

        const label patchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].start();

        forAll(fluidPatchAcceleration, I)
        {
            fluidPatchAcceleration[I] =
                fluidZoneAcceleration
                [
                    fluidMesh().faceZones()[fluidZoneIndex()].whichFace
                    (
                        patchStart + I
                    )
                ];
        }

        const_cast<elasticWallPressureFvPatchScalarField&>
        (
            refCast<const elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchIndex()]
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

        //- pass to all procs
        reduce(fluidZonePointsDispl, sumOp<vectorField>());

        const label globalFluidZoneIndex =
            findIndex(fluid().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "void Foam::fluidSolidInterface::syncFluidZonePointsDispl\n"
                "(\n"
                "    vectorField& fluidZonePointsDispl\n"
                ")"
            )   << "global zone point map is not available"
                << abort(FatalError);
        }

        const labelList& map =
            fluid().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

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
    fluid().writeFields(runTime);
    solid().writeFields(runTime);
}

// ************************************************************************* //
