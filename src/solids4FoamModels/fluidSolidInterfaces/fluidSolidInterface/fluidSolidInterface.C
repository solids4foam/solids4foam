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
#ifndef OPENFOAMESIORFOUNDATION
    #include "tetPointFields.H"
    #include "fixedValueTetPolyPatchFields.H"
    #include "tetPolyPatchInterpolation.H"
    #include "tetFemMatrices.H"
    #include "newSubsetMotionSolverFvMesh.H"
    #include "newSubsetMotionSolverFvMesh.H"
#endif
#include "fixedValuePointPatchFields.H"
#include "ZoneIDs.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"
#include "RBFMeshMotionSolver.H"

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


void Foam::fluidSolidInterface::calcInterfaceToInterfaceList() const
{
    if (interfaceToInterfaceList_.size())
    {
        FatalErrorIn
        (
            "void Foam::fluidSolidInterface::"
            "calcInterfaceToInterfaceList() const"
        )   << "List already set!" << abort(FatalError);
    }

    interfaceToInterfaceList_.setSize(nGlobalPatches_);

    // To maintain backwards compatibility, we will add a default dict
    {
        dictionary emptyDict;
        fsiProperties_.add("GGICoeffs", emptyDict);
        fsiProperties_.add("AMICoeffs", emptyDict);
        fsiProperties_.add("RBFCoeffs", emptyDict);
        fsiProperties_.add("directMapCoeffs", emptyDict);
    }

    // Create each interface-to-interface object
    // Note: the interpolation/mapping for each interface pair can be different
    for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Lookup the type
        const word type = fsiProperties_.lookupOrDefault<word>
#ifdef OPENFOAMESIORFOUNDATION
        (
            "interfaceTransferMethod", "AMI"
        );
#else
        (
            "interfaceTransferMethod", "GGI"
        );
#endif

        interfaceToInterfaceList_.set
        (
            interfaceI,
            interfaceToInterfaceMapping::New
            (
                type,
                fsiProperties_.subDict(type + "Coeffs"),
                fluidMesh().boundaryMesh()[fluidPatchIndices()[interfaceI]],
                solidMesh().boundaryMesh()[solidPatchIndices()[interfaceI]],
                fluid().globalPatches()[interfaceI],
                solid().globalPatches()[interfaceI]
            )
        );
    }
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

        if
        (
#ifdef OPENFOAMESIORFOUNDATION
            accumulatedFluidInterfaceDisplacementHeader.typeHeaderOk
            <
            vectorIOField
            >
            (
                true
            )
#else
            accumulatedFluidInterfaceDisplacementHeader.headerOk()
#endif
        )
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
    interfaceToInterfaceList_(),
    outerCorrTolerance_
    (
        fsiProperties_.lookupOrDefault<scalar>("outerCorrTolerance", 1e-06)
    ),
    nOuterCorr_
    (
        fsiProperties_.lookupOrDefault<int>("nOuterCorr", 30)
    ),
    additionalMeshCorrection_
    (
        fsiProperties_.lookupOrDefault<Switch>
        (
            "additionalMeshCorrection", false
        )
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
    accumulatedFluidInterfacesDisplacementsList_()
{
    Info<< "additionalMeshCorrection: " << additionalMeshCorrection_ << endl;

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

    // Set the number of global poly patches: solid or fluid
    nGlobalPatches_ = fluid().globalPatches().size();

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

    initializeFields();

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

    // Force creation of interface-to-interface objects as they may need to read
    // fields on restart
    interfaceToInterfaceList();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::~fluidSolidInterface()
{}


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


const Foam::PtrList<Foam::interfaceToInterfaceMapping>&
Foam::fluidSolidInterface::interfaceToInterfaceList() const
{
    if (!interfaceToInterfaceList_.size())
    {
        calcInterfaceToInterfaceList();
    }

    return interfaceToInterfaceList_;
}


Foam::vector Foam::fluidSolidInterface::totalForceOnInterface
(
    const standAlonePatch& zone, const vectorField& zoneTraction
) const
{
    const vectorField& localPoints = zone.localPoints();
    const faceList& localFaces = zone.localFaces();

    // Calculate face area vectors
    vectorField S(localFaces.size(), vector::zero);
    forAll(S, faceI)
    {
#ifdef OPENFOAMFOUNDATION
        S[faceI] = localFaces[faceI].area(localPoints);
#else
        S[faceI] = localFaces[faceI].normal(localPoints);
#endif
    }

    // No need for global sum as the zone is already global
    return sum(zoneTraction*mag(S));
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

        maxResidualsNorm_[interfaceI] = 0;

        maxIntsDisplsNorm_[interfaceI] = 0;

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
    if (interfaceToInterfaceList_.empty())
    {
        interfaceToInterfaceList();
    }
    else if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
            // Clear current interpolators
            interfaceToInterfaceList_.clear();

            // Re-create global patches
            fluid().clearGlobalPatches();
            solid().clearGlobalPatches();
            fluid().makeGlobalPatches(fluidPatchNames_);
            solid().makeGlobalPatches(solidPatchNames_);

            // Re-create interpolators
            interfaceToInterfaceList();
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

    if (maxDelta < interfaceDeformationLimit())
    {
        // Move only interface points
#ifdef OPENFOAMESIORFOUNDATION
        pointField newPoints = fluidMesh().points();
#else
        pointField newPoints = fluidMesh().allPoints();
#endif

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

        // PC: it is not good that this is hard-coded
        // A better way is to create an fsi point patch boundary condition that
        // knows how to lookup the motion from the fsi class or similar
        // For now, we will leave it

        // If the motionU field is in the object registry then we assume that
        // the fe motion solver is being used
#ifdef OPENFOAMESIORFOUNDATION
        const bool feMotionSolver = false;
#else
        const bool feMotionSolver =
            fluidMesh().foundObject<tetPointVectorField>("motionU");
#endif

        // If the pointMotionU field is in the object registry then we assume
        // that the fv motion solver is being used
        const bool fvMotionSolver =
            fluidMesh().foundObject<pointVectorField>("pointMotionU");

#ifndef OPENFOAMFOUNDATION
        // Check for RBF motion solver
        bool rbfMotionSolver = false;
        if (fluidMesh().foundObject<motionSolver>("dynamicMeshDict"))
        {
            rbfMotionSolver = isA<RBFMeshMotionSolver>
            (
                fluidMesh().lookupObject<motionSolver>("dynamicMeshDict")
            );
        }
#endif

        // Set motion on FSI interface
        if (feMotionSolver)
        {
#ifdef OPENFOAMESIORFOUNDATION
            notImplemented("Not implemented for this version of OpenFOAM/FOAM");
#else
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
                        fluidPatchesPointsDispls[interfaceI]
                      - fluidPatchesPointsDisplsPrev[interfaceI]
                    )/fluid().runTime().deltaT().value();
            }
#endif
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
#ifdef OPENFOAMESIORFOUNDATION
                        motionU.boundaryFieldRef()
                        [
                            fluidPatchIndices()[interfaceI]
                        ]
#else
                        motionU.boundaryField()[fluidPatchIndices()[interfaceI]]
#endif
                    );

                motionUFluidPatch ==
                    (
                        fluidPatchesPointsDispls[interfaceI]
                      - fluidPatchesPointsDisplsPrev[interfaceI]
                    )/fluid().runTime().deltaT().value();
            }
        }
#ifndef OPENFOAMESIORFOUNDATION
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
#endif
#ifndef OPENFOAMFOUNDATION
        else if (rbfMotionSolver)
        {
            // Prepare list of patch motions
            Field<vectorField> motion(fluidMesh().boundaryMesh().size());

            // Initialise all fields to zero
            forAll(fluidMesh().boundaryMesh(), patchI)
            {
                motion[patchI] = vectorField
                (
                    fluidMesh().boundaryMesh()[patchI].size(), vector::zero
                );
            }

            // Loop through all FSI interfaces
            forAll(fluid().globalPatches(), interfaceI)
            {
                // Interpolate the FSI interface point motion to the faces
                const vectorField interfacePatchMotion =
                    fluidPatchesPointsDispls[interfaceI]
                  - fluidPatchesPointsDisplsPrev[interfaceI];

                // Create interpolator
                primitivePatchInterpolation interp
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndices()[interfaceI]]
                );

                // Set motion of FSI interface
                motion[fluidPatchIndices()[interfaceI]] =
                    interp.pointToFaceInterpolate(interfacePatchMotion);
            }

            // Set motion field in RBF motion solver
            // Note: take displacement as opposed to velocity
            const_cast<RBFMeshMotionSolver&>
            (
                fluidMesh().lookupObject<RBFMeshMotionSolver>("dynamicMeshDict")
            ).setMotion(motion);
        }
#endif
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

    Info<< "Setting traction on solid interfaces" << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate total traction of fluid zone
        vectorField fluidZoneTotalTraction =
            fluid().faceZoneViscousForce(interfaceI)
          - fluid().faceZonePressureForce(interfaceI)*fluidZone.faceNormals();

        // Initialise the solid zone traction field that is to be interpolated
        // from the fluid zone
        vectorField solidZoneTotalTraction(solidZone.size(), vector::zero);

        // Transfer the field frm the fluid interface to the solid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZoneTotalTraction,    // from field
            solidZoneTotalTraction     // to field
        );

        // Flip traction sign after transferring from fluid to solid
        solidZoneTotalTraction = -solidZoneTotalTraction;

        // Set traction on solid
        if (coupled())
        {
            solid().setTraction
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZoneTotalTraction
            );
        }

        // Print total force on solid and fluid interfaces
        Info<< "Total force on fluid interface " << interfaceI << ": "
            << totalForceOnInterface(fluidZone, fluidZoneTotalTraction) << nl
            << "Total force on solid interface " << interfaceI << ": "
            << totalForceOnInterface(solidZone, solidZoneTotalTraction) << nl
            << endl;

        // Set interface pressure for elasticWallPressure boundary condition
        const label fluidPatchID = fluidPatchIndices()[interfaceI];
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchID]
            )
        )
        {
            scalarField& prevPressure =
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().p().boundaryField()[fluidPatchID]
                    )
                ).prevPressure();

            if (coupled())
            {
                prevPressure = fluid().patchPressureForce(fluidPatchID);
            }
            else
            {
                prevPressure = 0;
            }
        }
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    // Maximum residual for all interfaces
    scalar maxResidual = 0;

    for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate the point displacements of the solid interface
        const vectorField solidZonePointsDisplsAtSolid =
            solid().faceZonePointDisplacementIncrement(interfaceI);
        const vectorField solidZonePointsTotDisplsAtSolid =
            solid().faceZonePointDisplacementOld(interfaceI);

        // Initialise point displacement field at fluid interface
        vectorField solidZonePointsTotDispl
        (
            solidZonesPointsDispls()[interfaceI].size(), vector::zero
        );

        // Transfer displacement field from the solid to the fluid
        interfaceToInterfaceList()[interfaceI].transferPointsZoneToZone
        (
            solidZone,                              // from zone
            fluidZone,                              // to zone
            solidZonePointsDisplsAtSolid,           // from field
            solidZonesPointsDispls()[interfaceI]    // to field
        );

        interfaceToInterfaceList()[interfaceI].transferPointsZoneToZone
        (
            solidZone,                              // from zone
            fluidZone,                              // to zone
            solidZonePointsTotDisplsAtSolid,        // from field
            solidZonePointsTotDispl                 // to field
        );

        // Update interface residuals
        residualsPrev()[interfaceI] = residuals()[interfaceI];
        residuals()[interfaceI] =
            solidZonesPointsDispls()[interfaceI]
          - fluidZonesPointsDispls()[interfaceI];

        // We will use two definitions of residual
        scalar residualNorm1 = ::sqrt(gSum(magSqr(residuals()[interfaceI])));
        scalar residualNorm2 = residualNorm1;

        if (residualNorm1 > maxResidualsNorm_[interfaceI])
        {
            maxResidualsNorm_[interfaceI] = residualNorm1;
        }

        residualNorm1 /= maxResidualsNorm_[interfaceI] + SMALL;

        Info<< "FSI relative residual1 norm for interface " << interfaceI
            << ": " << residualNorm1 << endl;

        interfacesPointsDisplsPrev_[interfaceI] =
            interfacesPointsDispls_[interfaceI];

        interfacesPointsDispls_[interfaceI] =
            solidZonesPointsDispls()[interfaceI];

        const vectorField intTotDispl =
            interfacesPointsDispls_[interfaceI] + solidZonePointsTotDispl;

        const scalar intTotDisplNorm = Foam::sqrt(gSum(magSqr(intTotDispl)));

        if (intTotDisplNorm > maxIntsDisplsNorm_[interfaceI])
        {
            maxIntsDisplsNorm_[interfaceI] = intTotDisplNorm;
        }

        residualNorm2 /= maxIntsDisplsNorm_[interfaceI] + SMALL;

        Info<< "FSI residual2 norm for interface " << interfaceI
            << ": " << residualNorm2 << endl;

        // For this interface, the residual is defined as the minium of
        // residualNorm1 and residualNorm2
        const scalar residualInterfaceI = min(residualNorm1, residualNorm2);

        // Update the maximum residual for all interfaces
        maxResidual = max(maxResidual, residualInterfaceI);
    }

    return maxResidual;
}


void Foam::fluidSolidInterface::updateMovingWallPressureAcceleration()
{
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
            Info<< "Setting acceleration at fluid side of the interface"
                << endl;

            // Take references to zones
            const standAlonePatch& fluidZone =
               fluid().globalPatches()[interfaceI].globalPatch();
            const standAlonePatch& solidZone =
               solid().globalPatches()[interfaceI].globalPatch();

            const vectorField solidZoneAcceleration =
               solid().faceZoneAcceleration(interfaceI);

            // Initialise the fluid zone acceleration field that is to be
            // interpolated from the solid zone
            vectorField fluidZoneAcceleration(fluidZone.size(), vector::zero);

            // Transfer the field from the fluid interface to the solid
            // interface
            interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
            (
                solidZone,                // from zone
                fluidZone,                // to zone
                solidZoneAcceleration,    // from field
                fluidZoneAcceleration     // to field
            );

            const vectorField fluidPatchAcceleration =
                fluid().globalPatches()
                [
                   interfaceI
                ].globalFaceToPatch(fluidZoneAcceleration);

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
            Info<< "Setting acceleration at fluid side of the interface"
                << endl;

            // Take references to zones
            const standAlonePatch& fluidZone =
               fluid().globalPatches()[interfaceI].globalPatch();
            const standAlonePatch& solidZone =
               solid().globalPatches()[interfaceI].globalPatch();

            const vectorField solidZoneAcceleration =
               solid().faceZoneAcceleration(interfaceI);

            // Initialise the fluid zone acceleration field that is to be
            // interpolated from the solid zone
            vectorField fluidZoneAcceleration(fluidZone.size(), vector::zero);

            // Transfer the field from the fluid interface to the solid
            // interface
            interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
            (
                solidZone,                // from zone
                fluidZone,                // to zone
                solidZoneAcceleration,    // from field
                fluidZoneAcceleration     // to field
            );

            const vectorField fluidPatchAcceleration =
                fluid().globalPatches()
                [
                    interfaceI
                ].globalFaceToPatch(fluidZoneAcceleration);

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
        for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
        {
            if (!Pstream::master())
            {
                fluidZonesPointsDispls[interfaceI] = vector::zero;
            }

            // pass to all procs
            reduce(fluidZonesPointsDispls[interfaceI], sumOp<vectorField>());

            const labelList& map =
                fluid().globalPatches()
                [
                    interfaceI
                ].globalMasterToCurrentProcPointAddr();

            if (!Pstream::master())
            {
                const vectorField fluidZonePointsDisplGlobal =
                    fluidZonesPointsDispls[interfaceI];

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
    // Note: this means if the fluid defines new tmeporary fields within the
    // writeField function then they will not be created/written
    //fluid().writeFields(runTime);
    solid().writeFields(runTime);
}

// ************************************************************************* //
