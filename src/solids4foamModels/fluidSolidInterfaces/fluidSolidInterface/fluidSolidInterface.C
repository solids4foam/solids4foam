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
//#include "tetPointFields.H"
//#include "fixedValueTetPolyPatchFields.H"
//#include "tetPolyPatchInterpolation.H"
//#include "tetFemMatrices.H"
#include "fixedValuePointPatchFields.H"
#include "ZoneIDs.H"
//#include "newGgiInterpolation.H"
//#include "RBFMotionSolver.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

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
        )   << "Current solid zone points alarady exist"
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
//        new PrimitivePatch<face, List, const pointField&>
        new PrimitivePatch<face, SubList, const pointField&>//SubList is required for AMI!!!
        (
//            solidMesh().faceZones()[solidZoneIndex_]().localFaces(),
	    SubList<face>
            (
                solidMesh().faceZones()[solidZoneIndex_]().localFaces(),
                solidMesh().faceZones()[solidZoneIndex_]().size()
            ),
            currentSolidZonePoints()
        );
}

void Foam::fluidSolidInterface::calcCurrentFluidZonePoints() const
{
    // Find global face zones
    if (currentFluidZonePointsPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentFluidZonePoints() const"
        )   << "Current fluid zone points alarady exist"
            << abort(FatalError);
    }

    currentFluidZonePointsPtr_ =
        new vectorField(fluid().currentFaceZonePoints(fluidZoneIndex()));
}

void Foam::fluidSolidInterface::calcCurrentFluidZonePatch() const
{
    // Find global face zones
    if (currentFluidZonePatchPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcCurrentFluidZonePatch() const"
        )   << "Current fluid zone patch already exists"
            << abort(FatalError);
    }

    currentFluidZonePatchPtr_ =
        new PrimitivePatch<face, SubList, const pointField&>//SubList is required for AMI!!!
        (
	    SubList<face>
            (
                fluidMesh().faceZones()[solidZoneIndex_]().localFaces(),
                fluidMesh().faceZones()[solidZoneIndex_]().size()
            ),
            currentFluidZonePoints()
        );
}



void Foam::fluidSolidInterface::calcAMIInterpolator() const
{
    // Create extended ggi interpolation
////    if (ggiInterpolatorPtr_)
////    {
////        FatalErrorIn
////        (
////            "void fluidSolidInterface::"
////            "calcGgiInterpolator() const"
////        )   << "Ggi interpolator already exists"
////            << abort(FatalError);
////    }

    // Create AMI interpolation
////    if (AMIPtr_)//Do we need this?
////    {
////        FatalErrorIn
////        (
////            "void fluidSolidInterface::"
////            "calcAMIInterpolator() const"
////        )   << "AMI interpolator already exists"
////            << abort(FatalError);
////    }

    // Create copy of solid face zone primitive patch in current configuration

    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentSolidZonePointsPtr_);

    //currentSolidZonePatch().movePoints(currentSolidZonePoints());//Do we need this?

    Info<< "Create AMI" << endl;
    const bool flip = false;//good question

	    AMIPtr_.set
            (
                new AMIPatchToPatchInterpolation
                (
                    currentFluidZonePatch(),//fluidMesh().faceZones()[fluidZoneIndex_](),//p,
                    currentSolidZonePatch(),//solidMesh().boundaryMesh(),//nbrP,
                    faceAreaIntersect::tmMesh,
                    true,
                    AMIPatchToPatchInterpolation::imFaceAreaWeight,
                    -1,
                    flip
                )
            );
    Info<< "AMI Created" << endl;

    currentSolidZonePatch(),
    currentFluidZonePatch(),

////    ggiInterpolatorPtr_ =
////        new extendedGgiZoneInterpolation
////        (
////            fluidMesh().faceZones()[fluidZoneIndex_](),
////            currentSolidZonePatch(),
////            tensorField(0),
////            tensorField(0),
////            vectorField(0), // Slave-to-master separation. Bug fix
////            0,              // Non-overlapping face tolerances
////            0,              // HJ, 24/Oct/2008
////            true,           // Rescale weighting factors.  Bug fix, MB.
////            newGgiInterpolation::AABB
            // N_SQUARED BB_OCTREE AABB THREE_D_DISTANCE
            // Octree search, MB.
////        );


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
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres = AMI().interpolateToSource(fluidZoneFaceCentres);

////        vectorField solidZoneFaceCentres =
////            ggiInterpolatorPtr_->masterToSlave
////            (
////                fluidZoneFaceCentres
////            );

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

    Info<< "Checking solid-to-fluid point interpolator (AMI)" << endl;
    {
        vectorField solidZonePoints_ =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

//////        vectorField solidZonePoints = AMI().interpolateToSource(solidZonePoints_);
        vectorField solidZonePoints = solidZonePoints_;//Is this enough??? Otherwise we'll need an AMI for point interpolation!

////        vectorField solidZonePoints =
////            ggiInterpolatorPtr_->slaveToMasterPointInterpolate
////            (
////                solidZonePoints_
////            );

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

        Info<< "Solid-to-fluid point interpolation error (AMI): " << maxDist
            << endl;
    }

////    Info<< "Number of uncovered master faces: "
////        << ggiInterpolatorPtr_->uncoveredMasterFaces().size() << endl;

////    Info<< "Number of uncovered slave faces: "
////        << ggiInterpolatorPtr_ ->uncoveredSlaveFaces().size() << endl;
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
    dynamicFvMesh& fluidMesh,
    dynamicFvMesh& solidMesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fsiProperties",
            fluidMesh.time().constant(),
            fluidMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    fsiProperties_(subDict(type + "Coeffs")),
    fluidMesh_(fluidMesh),
    fluid_(fluidModel::New(fluidMesh_)),
    solidMesh_(solidMesh),
    solid_(solidModel::New(solidMesh_)),
    solidPatchIndex_(-1),
    solidZoneIndex_(-1),
    fluidPatchIndex_(-1),
    fluidZoneIndex_(-1),
    currentSolidZonePointsPtr_(NULL),
    currentSolidZonePatchPtr_(NULL),
    currentFluidZonePointsPtr_(NULL),
    currentFluidZonePatchPtr_(NULL),
    AMIPtr_(NULL),
////    ggiInterpolatorPtr_(NULL),
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
    interfaceDeformationLimit_
    (
        fsiProperties_.lookupOrDefault<scalar>("interfaceDeformationLimit", 0.0)
    ),
    fluidZonePointsDispl_(),
    fluidZonePointsDisplRef_(),
    fluidZonePointsDisplPrev_(),
    solidZonePointsDispl_(),
    solidZonePointsDisplRef_(),
    solidZonePressure_(),
    residual_(),
    residualPrev_(),
    maxResidualNorm_(0),
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
        solidMesh.boundaryMesh()
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
        solidMesh.faceZones()
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
        fluidMesh.boundaryMesh()
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
        fluidMesh.faceZones()
    );

    if (!fluidZone.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Fluid face zone name " << fluidZoneName
            << " not found.  Please check your face zone definition."
            << abort(FatalError);
    }

    fluidZoneIndex_ = fluidZone.index();

    // Initialize solid zone pressure
    solidZonePressure_ =
        scalarField(solidMesh.faceZones()[solidZoneIndex()].size(), 0.0);

    // Initialize residual
    residual_ =
        vectorField
        (
            fluidMesh.faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::~fluidSolidInterface()
{
    deleteDemandDrivenData(currentSolidZonePointsPtr_);
    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentFluidZonePointsPtr_);
    deleteDemandDrivenData(currentFluidZonePatchPtr_);
    AMIPtr_.clear();
////    deleteDemandDrivenData(AMIPtr_);
////    deleteDemandDrivenData(ggiInterpolatorPtr_);
    deleteDemandDrivenData(accumulatedFluidInterfaceDisplacementPtr_);
    deleteDemandDrivenData(minEdgeLengthPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::vectorField&
Foam::fluidSolidInterface::currentSolidZonePoints() const
{
    if (!currentSolidZonePointsPtr_)
    {
        calcCurrentSolidZonePoints();
    }

    return *currentSolidZonePointsPtr_;
}


//const Foam::PrimitivePatch<Foam::face, Foam::List, const Foam::pointField&>&
const Foam::PrimitivePatch<Foam::face, Foam::SubList, const Foam::pointField&>&
Foam::fluidSolidInterface::currentSolidZonePatch() const
{
    if (!currentSolidZonePatchPtr_)
    {
        calcCurrentSolidZonePatch();
    }

    return *currentSolidZonePatchPtr_;
}

const Foam::vectorField&
Foam::fluidSolidInterface::currentFluidZonePoints() const
{
    if (!currentFluidZonePointsPtr_)
    {
        calcCurrentFluidZonePoints();
    }

    return *currentFluidZonePointsPtr_;
}

const Foam::PrimitivePatch<Foam::face, Foam::SubList, const Foam::pointField&>&
Foam::fluidSolidInterface::currentFluidZonePatch() const
{
    if (!currentFluidZonePatchPtr_)
    {
        calcCurrentFluidZonePatch();
    }

    return *currentFluidZonePatchPtr_;
}


const Foam::AMIPatchToPatchInterpolation& Foam::fluidSolidInterface::AMI() const
{
/*    if (!owner())//might be important to insert this! - example in src/meshTools/regionCoupled/patches/regionCoupledPolyPatch/regionCoupledBase.C/.H
    {
        FatalErrorIn
        (
            "const AMIPatchToPatchInterpolation& cyclicAMIPolyPatch::AMI()"
        )
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_.valid())
    {
        resetAMI();
    }
*/

    if (!AMIPtr_.valid())
    {
        calcAMIInterpolator();
    }

    return AMIPtr_;
}

////const Foam::extendedGgiZoneInterpolation&
////Foam::fluidSolidInterface::ggiInterpolator() const
////{
////    if (!ggiInterpolatorPtr_)
////    {
////        calcGgiInterpolator();
////    }

////    return *ggiInterpolatorPtr_;
////}


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

Info  << "Initialization finished..." << endl;

    // PC: do we need to re-read these?
    // nOuterCorr_ = fsiProperties_.lookupOrDefault<int>("nOuterCorr", 30);

    // outerCorrTolerance_ =
    //     fsiProperties_.lookupOrDefault<scalar>("outerCorrTolerance", 1e-6);

    //coupled_ = fsiProperties_.lookupOrDefault<Switch>("coupled", true);
}


void Foam::fluidSolidInterface::updateInterpolator()
{
    if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
	    AMIPtr_.clear();
	    AMI();
	    Info << "AMI was updated... (1)" << endl;
////            ggiInterpolator();
        }
    }
    else
    {
        if ((runTime().timeIndex() - 1) == 0)
        {
    	    AMIPtr_.clear();
	    AMI();
	    Info << "AMI was updated... (2)" << endl;
//            ggiInterpolator();
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
        pointField newPoints = fluidMesh().points();

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

        fluidMesh_.movePoints(newPoints);

        // Accumulate interface points displacement
        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;
    }
    else
    {
        // Move whole fluid mesh
        pointField newPoints = fluidMesh().points();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll (accumulatedFluidInterfaceDisplacement(), pointI)
        {
            newPoints[meshPoints[pointI]] -=
                accumulatedFluidInterfaceDisplacement()[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);

        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;

        // Check mesh motion solver type
////        bool feMotionSolver =
////            fluidMesh().objectRegistry::foundObject<tetPointVectorField>
////            (
////                "motionU"
////            );

        bool fvMotionSolver =
            fluidMesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        // bool rbfMotionSolver =
        //     fluidMesh().objectRegistry::foundObject<RBFMotionSolver>
        //     (
        //         "dynamicMeshDict"
        //     );

//         if (rbfMotionSolver)
//         {
//             // Grab RBF motion solver
//             RBFMotionSolver& ms =
//                 const_cast<RBFMotionSolver&>
//                 (
//                     fluidMesh().objectRegistry::lookupObject<RBFMotionSolver>
//                     (
//                         "dynamicMeshDict"
//                     )
//                 );

//             Info<< "RBF mesh motion" << endl;

//             const labelList& movingMeshPoints = ms.movingIDs();

//             vectorField motion(movingMeshPoints.size(), vector::zero);

//             vectorField fluidPatchDisplacement =
//                 accumulatedFluidInterfaceDisplacement();
// //                /fluid().runTime().deltaT().value();

//             const labelList& meshPoints =
//                 fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

//             forAll(meshPoints, pointI)
//             {
//                 label curMovingPoint =
//                     findIndex(movingMeshPoints, meshPoints[pointI]);

//                 if (curMovingPoint != -1)
//                 {
//                     motion[curMovingPoint] = fluidPatchDisplacement[pointI];
//                 }
//             }

//             ms.setMotion(motion);
//         }
//         else if (feMotionSolver)
////        if (feMotionSolver)
////        {
////            tetPointVectorField& motionU =
////                const_cast<tetPointVectorField&>
////                (
////                    fluidMesh().objectRegistry::
////                    lookupObject<tetPointVectorField>
////                    (
////                        "motionU"
////                    )
////                );

////            fixedValueTetPolyPatchVectorField& motionUFluidPatch =
////                refCast<fixedValueTetPolyPatchVectorField>
////                (
////                    motionU.boundaryField()[fluidPatchIndex()]
////                );
////
////            tetPolyPatchInterpolation tppi
////            (
////                refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
////            );

////            motionUFluidPatch ==
////                tppi.pointToPointInterpolate
////                (
////                    accumulatedFluidInterfaceDisplacement()
////                   /fluid().runTime().deltaT().value()
////                );
////        }
////        else if (fvMotionSolver)
	if (fvMotionSolver)
        {

	    Info << "Using fvMotionSolver for movment of fluid mesh..." << endl;

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
                accumulatedFluidInterfaceDisplacement()
               /fluid().runTime().deltaT().value();
        }
        else
        {
            FatalErrorIn("fluidSolidInterface::moveFluidMesh()")
                << "Problem with fluid mesh motion solver selection"
                    << abort(FatalError);
        }

        fluidMesh_.update();

        accumulatedFluidInterfaceDisplacement() =
            vectorField
            (
                accumulatedFluidInterfaceDisplacement().size(),
                vector::zero
            );
    }


    // Move unused fluid mesh points
    {
        vectorField newPoints = fluidMesh().points();

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

        fluidMesh_.movePoints(newPoints);
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
    const vectorField& p =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();
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

    vectorField solidZoneTraction = AMI().interpolateToTarget(-fluidZoneTraction);

////    vectorField solidZoneTraction =
////        ggiInterpolator().masterToSlave
////        (
////           -fluidZoneTraction
////        );

    vectorField solidZoneTotalTraction = AMI().interpolateToTarget(-fluidZoneTotalTraction);

////    const vectorField solidZoneTotalTraction =
////        ggiInterpolator().masterToSlave
////        (
////            -fluidZoneTotalTraction
////        );

    const scalarField solidZoneMuEff = 
			AMI().interpolateToTarget
			(
		     	    fluid().faceZoneMuEff(fluidZoneIndex(), fluidPatchIndex())
			);

////    const scalarField solidZoneMuEff =
////        ggiInterpolator().masterToSlave
////        (
////            fluid().faceZoneMuEff(fluidZoneIndex(), fluidPatchIndex())
////        );

    const tensorField solidZoneSurfaceGradientOfVelocity =
        solid().faceZoneSurfaceGradientOfVelocity
        (
            solidZoneIndex(),
            solidPatchIndex()
        );

    const vectorField solidZoneNormal =
        solid().faceZoneNormal
        (
            solidZoneIndex(),
            solidPatchIndex()
        );

    // Add part of the viscous force present only
    // at the deforming and moving walls
    solidZoneTraction +=
        solidZoneMuEff
       *(
           -2*tr(solidZoneSurfaceGradientOfVelocity)*solidZoneNormal
          + (solidZoneSurfaceGradientOfVelocity&solidZoneNormal)
        );

    const vectorField movingTraction =
        solidZoneMuEff
       *(
           -2*tr(solidZoneSurfaceGradientOfVelocity)*solidZoneNormal
          + (solidZoneSurfaceGradientOfVelocity&solidZoneNormal)
        );

    solidZonePressure_ = AMI().interpolateToTarget(fluidZonePressure);

////    solidZonePressure_ =
////        ggiInterpolator().masterToSlave
////        (
////            fluidZonePressure
////        );

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
    }

    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        const vector totalTractionForce = sum(fluidZoneTotalTraction*mag(S));

        Info<< "Total force (fluid) = "
            << totalTractionForce << endl;
    }

    // Total force at the solid side of the interface
    {
        const vectorField& p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        const vector totalTractionForce = sum(solidZoneTotalTraction*mag(S));

        Info<< "Total force (solid) = "
            << totalTractionForce << endl;
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    vectorField solidZonePointsDisplAtSolid =
        solid().faceZonePointDisplacementIncrement(solidZoneIndex());

    solidZonePointsDispl() = AMI().interpolateToSource(solidZonePointsDisplAtSolid);

////    solidZonePointsDispl() =
////        ggiInterpolator().slaveToMasterPointInterpolate
////        (
////            solidZonePointsDisplAtSolid
////        );

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

    scalar residualNorm = ::sqrt(sum(magSqr(residual())));

    if (residualNorm > maxResidualNorm_)
    {
        maxResidualNorm_ = residualNorm;
    }

    residualNorm /= maxResidualNorm_ + SMALL;

    Info<< "Current FSI relative residual norm: " << residualNorm << endl;

    return residualNorm;
}


// ************************************************************************* //

