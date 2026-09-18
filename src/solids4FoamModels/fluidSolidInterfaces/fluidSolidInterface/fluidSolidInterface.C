/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fluidSolidInterface.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"
#ifndef OPENFOAM_NOT_EXTEND
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
#include "addToRunTimeSelectionTable.H"
#ifndef S4F_NO_RBF
    #include "RBFMeshMotionSolver.H"
#endif
#include "FieldSumOp.H"
#ifdef OPENFOAM_COM
    #include "dynamicMotionSolverFvMesh.H"
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSolidInterface, 0);
    defineRunTimeSelectionTable(fluidSolidInterface, dictionary);
    addToRunTimeSelectionTable(physicsModel, fluidSolidInterface, physicsModel);

    namespace
    {
        word resolvedRegion
        (
            const Time& runTime,
            const dictionary& fsiProperties,
            const word& modelName,
            const word& defaultRegion
        )
        {
            const word controlDictRegion
            (
                runTime.controlDict().subOrEmptyDict(modelName)
                    .lookupOrDefault<word>("region", defaultRegion)
            );

            return fsiProperties.lookupOrDefault<word>
            (
                modelName + "Region",
                controlDictRegion
            );
        }


        void warnOnConflictingRegion
        (
            const Time& runTime,
            const dictionary& fsiProperties,
            const word& modelName,
            const word& resolvedRegionName,
            const word& defaultRegion
        )
        {
            if
            (
                fsiProperties.found(modelName + "Region")
             && runTime.controlDict().subOrEmptyDict(modelName).found("region")
            )
            {
                const word controlDictRegion
                (
                    runTime.controlDict().subOrEmptyDict(modelName)
                        .lookupOrDefault<word>("region", defaultRegion)
                );

                if (controlDictRegion != resolvedRegionName)
                {
                    WarningInFunction
                        << "Conflicting " << modelName << " region settings: "
                        << "using '" << resolvedRegionName << "' from "
                        << "fsiProperties and ignoring '" << controlDictRegion
                        << "' from controlDict/" << modelName
                        << "/region" << endl;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fluidSolidInterface::updateCoupled()
{
    if (couplingStartTime_ > SMALL && !coupled_)
    {
        if (runTime().value() > (couplingStartTime_ - SMALL))
        {
            InfoIn("fluidSolidInterface::updateCoupled()")
                << "Enabling fluid-solid coupling" << nl << endl;

            // Enable coupling
            coupled_ = true;

            return true;
        }
    }

    return false;
}


bool Foam::fluidSolidInterface::newTimeStep() const
{
    if (curTimeIndex_ != runTime().timeIndex())
    {
        curTimeIndex_ = runTime().timeIndex();
        return true;
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
        if (!fsiProperties_.found("GGICoeffs"))
        {
            fsiProperties_.add("GGICoeffs", emptyDict);
        }
        if (!fsiProperties_.found("AMICoeffs"))
        {
            fsiProperties_.add("AMICoeffs", emptyDict);
        }
        if (!fsiProperties_.found("RBFCoeffs"))
        {
            fsiProperties_.add("RBFCoeffs", emptyDict);
        }
        if (!fsiProperties_.found("directMapCoeffs"))
        {
            fsiProperties_.add("directMapCoeffs", emptyDict);
        }
    }

    // Create each interface-to-interface object
    // Note: the interpolation/mapping for each interface pair can be different
    for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Lookup the type
        const word type = fsiProperties_.lookupOrAddDefault<word>
#ifdef OPENFOAM_NOT_EXTEND
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
#ifdef OPENFOAM_NOT_EXTEND
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
    physicsModel(type, runTime, region),
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
    fluid_(),
    solid_(),
    solidPatchNames_(),
    fluidPatchNames_(),
    solidPatchIndices_(),
    fluidPatchIndices_(),
    nGlobalPatches_(-1),
    interfaceToInterfaceList_(),
    outerCorrTolerance_
    (
        fsiProperties_.lookupOrAddDefault<scalar>("outerCorrTolerance", 1e-06)
    ),
    nOuterCorr_
    (
        fsiProperties_.lookupOrAddDefault<int>("nOuterCorr", 30)
    ),
    allowUnconvergedCoupling_
    (
        fsiProperties_.lookupOrAddDefault<Switch>
        (
            "allowUnconvergedCoupling", false
        )
    ),
    additionalMeshCorrection_
    (
        fsiProperties_.lookupOrAddDefault<Switch>
        (
            "additionalMeshCorrection", false
        )
    ),
    coupled_
    (
        fsiProperties_.lookupOrAddDefault<Switch>("coupled", true)
    ),
    couplingStartTime_
    (
        fsiProperties_.lookupOrAddDefault<scalar>("couplingStartTime", -1.0)
    ),
    predictor_(fsiProperties_.lookupOrAddDefault<Switch>("predictor", false)),
    curTimeIndex_(-1),
    interfaceDeformationLimit_
    (
        fsiProperties_.lookupOrAddDefault<scalar>
        (
            "interfaceDeformationLimit", 0.0
        )
    ),
    fluidZonesPointsDispls_(),
    fluidZonesPointsDisplsRef_(),
    fluidZonesPointsDisplsPrev_(),
    solidZonesPointsDispls_(),
    solidZonesPointsDisplsRef_(),
    fluidZonesTractions_(),
    fluidZonesTractionsRef_(),
    solidZonesTractions_(),
    solidZonesTractionsRef_(),
    interfacesPointsDispls_(),
    interfacesPointsDisplsPrev_(),
    incrementalResiduals_
    (
        fsiProperties_.lookupOrAddDefault<Switch>("incrementalResiduals", true)
    ),
    requireAllResidualMeasures_
    (
        fsiProperties_.lookupOrAddDefault<Switch>
        (
            "requireAllResidualMeasures", false
        )
    ),
    robinPressureTolerance_
    (
        fsiProperties_.lookupOrAddDefault<scalar>
        (
            "robinPressureTolerance", outerCorrTolerance_
        )
    ),
    robinFluxTolerance_
    (
        fsiProperties_.lookupOrAddDefault<scalar>
        (
            "robinFluxTolerance", 10*outerCorrTolerance_
        )
    ),
    robinConvergence_
    (
        fsiProperties_.lookupOrAddDefault<word>
        (
            "robinConvergence", "residual"
        )
    ),
    robinSlowRate_
    (
        fsiProperties_.lookupOrAddDefault<scalar>("robinSlowRate", 0.9)
    ),
    robinStallRate_
    (
        fsiProperties_.lookupOrAddDefault<scalar>("robinStallRate", 0.97)
    ),
    robinStallIterations_
    (
        fsiProperties_.lookupOrAddDefault<label>("robinStallIterations", 4)
    ),
    robinStallToleranceFactor_
    (
        fsiProperties_.lookupOrAddDefault<scalar>
        (
            "robinStallToleranceFactor", 10
        )
    ),
    robinFluxReferenceVelocity_
    (
        fsiProperties_.lookupOrAddDefault<scalar>
        (
            "robinFluxReferenceVelocity", 0
        )
    ),
    robinPressureReference_
    (
        fsiProperties_.lookupOrAddDefault<scalar>
        (
            "robinPressureReference", 0
        )
    ),
    robinInterfaces_(),
    robinPressurePrev_(),
    maxRobinPressureNorm_(),
    robinPressureResidual_(0),
    robinFluxResidual_(0),
    maxRobinKinematicNorm_(0),
    robinDispHistory_(),
    robinPressureHistory_(),
    robinFluxHistory_(),
    robinConvergenceState_(0),
    robinKinematicConsistency_(true),
    robinRelaxationWarned_(false),
    robinReferenceHintIssued_(false),
    residuals_(),
    residualsPrev_(),
    maxResidualsNorm_(),
    maxIntsDisplsNorm_(),
    outerCorr_(0),
    writeResidualsToFile_
    (
        fsiProperties_.lookupOrAddDefault<Switch>("writeResidualsToFile", false)
    ),
    residualFilePtr_(),
    interpolatorUpdateFrequency_
    (
        fsiProperties_.lookupOrAddDefault<int>("interpolatorUpdateFrequency", 0)
    ),
    accumulatedFluidInterfacesDisplacementsList_()
{
    const word fluidRegion
    (
        resolvedRegion(runTime, fsiProperties_, "fluid", "fluid")
    );

    const word solidRegion
    (
        resolvedRegion(runTime, fsiProperties_, "solid", "solid")
    );

    warnOnConflictingRegion
    (
        runTime,
        fsiProperties_,
        "fluid",
        fluidRegion,
        "fluid"
    );

    warnOnConflictingRegion
    (
        runTime,
        fsiProperties_,
        "solid",
        solidRegion,
        "solid"
    );

    fluid_ = fluidModel::New(runTime, fluidRegion);
    solid_ = solidModel::New(runTime, solidRegion);

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
    fluidZonesTractions_.setSize(nGlobalPatches_);
    fluidZonesTractionsRef_.setSize(nGlobalPatches_);
    solidZonesTractions_.setSize(nGlobalPatches_);
    solidZonesTractionsRef_.setSize(nGlobalPatches_);
    interfacesPointsDispls_.setSize(nGlobalPatches_);
    interfacesPointsDisplsPrev_.setSize(nGlobalPatches_);
    residuals_.setSize(nGlobalPatches_);
    residualsPrev_.setSize(nGlobalPatches_);
    maxResidualsNorm_.setSize(nGlobalPatches_);
    maxIntsDisplsNorm_.setSize(nGlobalPatches_);
    robinInterfaces_.setSize(nGlobalPatches_, false);
    robinPressurePrev_.setSize(nGlobalPatches_);
    maxRobinPressureNorm_.setSize(nGlobalPatches_, 0);

    if (robinStallIterations_ < 0)
    {
        FatalErrorInFunction
            << "robinStallIterations must be non-negative (0 disables the"
            << " stall detection)" << abort(FatalError);
    }

    if (robinFluxReferenceVelocity_ < 0 || robinPressureReference_ < 0)
    {
        FatalErrorInFunction
            << "robinFluxReferenceVelocity and robinPressureReference must be"
            << " non-negative (0 disables them)" << abort(FatalError);
    }

    if
    (
        robinConvergence_ != "iterationError"
     && robinConvergence_ != "residual"
    )
    {
        FatalErrorInFunction
            << "Unknown robinConvergence " << robinConvergence_
            << ". Valid options are (residual iterationError)"
            << abort(FatalError);
    }

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
    // NB: dictionary must be unregistered to avoid adding to the database

    IOdictionary props
    (
        IOobject
        (
            "fsiProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // Do not register
        )
    );

    const word modelType(props.lookup("fluidSolidInterface"));

    Info<< "Selecting fluidSolidInterface " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            props,
            "fluidSolidInterface",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fluidSolidInterface::New(Time&, const word&)"
        )   << "Unknown fluidSolidInterface type " << modelType
            << endl << endl
            << "Valid fluidSolidInterface types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<fluidSolidInterface>(ctorPtr(runTime, region));
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
        fileName historyDir;
        if (Pstream::parRun())
        {
            historyDir = runTime().path()/".."/"postProcessing";
        }
        else
        {
            historyDir = runTime().path()/"postProcessing";
        }

        mkDir(historyDir);
        residualFilePtr_.set(new OFstream(historyDir/"fsiResiduals.dat"));
        residualFilePtr_() << "Time outerCorrector residual";
        if (hasRobinInterface())
        {
            residualFilePtr_()
                << " robinPressureResidual robinFluxResidual"
                << " robinConvergenceState";
        }
        residualFilePtr_() << endl;
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
#ifdef OPENFOAM_ORG
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
    maxRobinKinematicNorm_ = 0;
    robinPressureResidual_ = 0;
    robinFluxResidual_ = 0;
    robinDispHistory_.clear();
    robinPressureHistory_.clear();
    robinFluxHistory_.clear();
    robinConvergenceState_ = 0;

    //- Check is any of the patches have changed sizes
    bool needUpdate = false;
    forAll(fluid().globalPatches(), interfaceI)
    {
        const label nFaces = returnReduce
        (
            fluid().globalPatches()[interfaceI].patch().size(),
            sumOp<label>()
        );
        if (fluid().globalPatches()[interfaceI].globalPatch().size() != nFaces)
        {
            needUpdate = true;
            break;
        }
    }

    // The size of a patch has changed so clear the patches
    if (needUpdate)
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

    // Reset the point fields
    forAll(fluid().globalPatches(), interfaceI)
    {
        const label nPoints =
            fluid().globalPatches()[interfaceI].globalPatch().nPoints();
        const label nFluidFaces =
            fluid().globalPatches()[interfaceI].globalPatch().size();
        const label nSolidFaces =
            solid().globalPatches()[interfaceI].globalPatch().size();

        fluidZonesPointsDispls_[interfaceI] = vectorField(nPoints, vector::zero);

        fluidZonesPointsDisplsRef_[interfaceI] =
            vectorField(nPoints, vector::zero);

        fluidZonesPointsDisplsPrev_[interfaceI] =
            vectorField(nPoints, vector::zero);

        solidZonesPointsDispls_[interfaceI] =
            vectorField(nPoints, vector::zero);

        solidZonesPointsDisplsRef_[interfaceI] =
            vectorField(nPoints, vector::zero);

        fluidZonesTractions_[interfaceI] =
            vectorField(nFluidFaces, vector::zero);

        fluidZonesTractionsRef_[interfaceI] =
            vectorField(nFluidFaces, vector::zero);

        solidZonesTractions_[interfaceI] =
            vectorField(nSolidFaces, vector::zero);

        solidZonesTractionsRef_[interfaceI] =
            vectorField(nSolidFaces, vector::zero);

        residualsPrev_[interfaceI] = residuals_[interfaceI];

        residuals_[interfaceI] =
            vectorField(nPoints, vector::zero);

        maxResidualsNorm_[interfaceI] = 0;

        maxIntsDisplsNorm_[interfaceI] = 0;

        const label fluidPatchID = fluidPatchIndices()[interfaceI];
        robinInterfaces_[interfaceI] =
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().solutionP().boundaryField()[fluidPatchID]
            );
        maxRobinPressureNorm_[interfaceI] = 0;

        // Robin-Neumann coupling is designed for unrelaxed fixed-point
        // iterations; warn once otherwise (not during construction, when the
        // derived coupling interface is not yet available)
        if
        (
            robinInterfaces_[interfaceI]
         && runTime().timeIndex() > 0
         && !robinRelaxationWarned_
         && !fluidMeshFollowsSolid()
        )
        {
            WarningInFunction
                << "The elasticWallPressure (Robin) interface "
                << fluidMesh().boundary()[fluidPatchID].name()
                << " is used with an under-relaxed or accelerated interface"
                << " displacement (" << type() << "). The Robin-Neumann"
                << " coupling is designed for fixedRelaxation with"
                << " relaxationFactor 1: with relaxation the fluid mesh lags"
                << " the solid within the iterations, so the interface flux"
                << " cannot be made consistent with the mesh motion and the"
                << " leakage is only reported." << endl;

            robinRelaxationWarned_ = true;
        }

        if (robinInterfaces_[interfaceI])
        {
            robinPressurePrev_[interfaceI] =
                fluid().solutionP().boundaryField()[fluidPatchID];
        }
        else
        {
            robinPressurePrev_[interfaceI].clear();
        }

        interfacesPointsDispls_[interfaceI] =
            vectorField(nPoints, vector::zero);

        interfacesPointsDisplsPrev_[interfaceI] =
            vectorField(nPoints, vector::zero);

        if (accumulatedFluidInterfacesDisplacementsList_.size())
        {
            accumulatedFluidInterfacesDisplacementsList_[interfaceI] =
                vectorField
                (
                    fluidMesh().boundaryMesh()
                    [
                        fluidPatchIndices()[interfaceI]
                    ].nPoints(),
                    vector::zero
                );
        }
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
        const scalarField fluidZoneDeltaCoeffs
        (
            fluid().globalPatches()[interfaceI].patchFaceToGlobal
            (
                fluidMesh().boundary()
                [
                    fluid().globalPatches()[interfaceI].patch().index()
                ].deltaCoeffs()
            )
        );

        // Zone deltaCoeffs at points
        const scalarField fluidZonePointDeltaCoeffs
        (
            fluid().globalPatches()
            [
                interfaceI
            ].interpolator().faceToPointInterpolate(fluidZoneDeltaCoeffs)
        );

        // Patch deltaCoeffs at points
        const scalarField fluidPatchPointDeltaCoeffs
        (
            fluid().globalPatches()[interfaceI].globalPointToPatch
            (
                fluidZonePointDeltaCoeffs
            )
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
#ifdef OPENFOAM_NOT_EXTEND
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
#ifdef OPENFOAM_NOT_EXTEND
        const bool feMotionSolver = false;
#else
        const bool feMotionSolver =
            fluidMesh().foundObject<tetPointVectorField>("motionU");
#endif

        // If the pointMotionU field is in the object registry then we assume
        // that the fv motion solver is being used
        const bool fvMotionSolver =
            fluidMesh().foundObject<pointVectorField>("pointMotionU");

#if !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)
        // Check for RBF motion solver
        RBFMeshMotionSolver* rbfMotionSolverPtr = nullptr;
#ifdef OPENFOAM_COM
        if (isA<dynamicMotionSolverFvMesh>(fluidMesh()))
        {
            motionSolver& motion =
                const_cast<motionSolver&>
                (
                    refCast<const dynamicMotionSolverFvMesh>
                    (
                        fluidMesh()
                    ).motion()
                );

            if (isA<RBFMeshMotionSolver>(motion))
            {
                rbfMotionSolverPtr = &refCast<RBFMeshMotionSolver>(motion);
            }
        }
#else
        if (fluidMesh().foundObject<motionSolver>("dynamicMeshDict"))
        {
            motionSolver& motion =
                const_cast<motionSolver&>
                (
                    fluidMesh().lookupObject<motionSolver>("dynamicMeshDict")
                );

            if (isA<RBFMeshMotionSolver>(motion))
            {
                rbfMotionSolverPtr = &refCast<RBFMeshMotionSolver>(motion);
            }
        }
#endif
#endif

        // Set motion on FSI interface
        if (feMotionSolver)
        {
#ifdef OPENFOAM_NOT_EXTEND
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
                        boundaryFieldRef
                        (
                            motionU
                        )[fluidPatchIndices()[interfaceI]]
                    );

                motionUFluidPatch ==
                    (
                        fluidPatchesPointsDispls[interfaceI]
                      - fluidPatchesPointsDisplsPrev[interfaceI]
                    )/fluid().runTime().deltaT().value();
            }
        }
#ifndef OPENFOAM_NOT_EXTEND
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
#if !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)
        else if (rbfMotionSolverPtr)
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
                const vectorField interfacePatchMotion
                (
                    fluidPatchesPointsDispls[interfaceI]
                  - fluidPatchesPointsDisplsPrev[interfaceI]
                );

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
            rbfMotionSolverPtr->setMotion(motion);
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

    fluid().syncGlobalPatches();
}


void Foam::fluidSolidInterface::updateForce()
{
    Info<< "Setting traction on solid interfaces" << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate total traction of fluid zone
        vectorField fluidZoneTotalTraction
        (
            fluid().faceZoneViscousForce(interfaceI)
          - fluid().faceZonePressureForce(interfaceI)*fluidZone.faceNormals()
        );

        fluidZonesTractions()[interfaceI] = fluidZoneTotalTraction;

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

        solidZonesTractions()[interfaceI] = solidZoneTotalTraction;

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

        // ZT: move this part of the code into
        // updateElasticWallPressureAcceleration() function
        // // Set interface pressure for elasticWallPressure boundary condition
        // const label fluidPatchID = fluidPatchIndices()[interfaceI];
        // if
        // (
        //     isA<elasticWallPressureFvPatchScalarField>
        //     (
        //         fluid().p().boundaryField()[fluidPatchID]
        //     )
        // )
        // {
        //     scalarField& prevPressure =
        //         const_cast<elasticWallPressureFvPatchScalarField&>
        //         (
        //             refCast<const elasticWallPressureFvPatchScalarField>
        //             (
        //                 fluid().p().boundaryField()[fluidPatchID]
        //             )
        //         ).prevPressure();

        //     if (coupled())
        //     {
        //         prevPressure = fluid().patchPressureForce(fluidPatchID);
        //     }
        //     else
        //     {
        //         prevPressure = 0;
        //     }
        // }
    }
}

void Foam::fluidSolidInterface::updateViscousForceAndPressure()
{
    // Check if coupling switch needs to be updated
    if (!coupled_)
    {
        updateCoupled();
    }

    Info<< "Setting traction and pressure on solid interfaces" << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate total traction of fluid zone
        vectorField fluidZoneTraction
        (
            fluid().faceZoneViscousForce(interfaceI)
        );

        // Calculate pressure of fluid zone
        scalarField fluidZonePressure
        (
            fluid().faceZonePressureForce(interfaceI)
        );

        // Initialise the solid zone traction and pressure fields
        // that is to be interpolated from the fluid zone
        vectorField solidZoneTraction(solidZone.size(), vector::zero);
        scalarField solidZonePressure(solidZone.size(), 0);

        // Transfer the field frm the fluid interface to the solid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZoneTraction,         // from field
            solidZoneTraction          // to field
        );
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZonePressure,         // from field
            solidZonePressure          // to field
        );

        // Flip traction sign after transferring from fluid to solid
        solidZoneTraction = -solidZoneTraction;

        // Set traction on solid
        if (coupled())
        {
            solid().setTraction
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZoneTraction
            );

            solid().setPressure
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZonePressure
            );
        }

        // Print total viscous force on solid and fluid interfaces
        Info<< "Total viscous force on fluid interface " << interfaceI << ": "
            << totalForceOnInterface(fluidZone, fluidZoneTraction) << nl
            << "Total force on solid interface " << interfaceI << ": "
            << totalForceOnInterface(solidZone, solidZoneTraction) << nl
            << endl;
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    solid().syncGlobalPatches();

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
        const vectorField solidZonePointsDisplsAtSolid
        (
            solid().faceZonePointDisplacementIncrement(interfaceI)
        );
        const vectorField solidZonePointsTotDisplsAtSolid
        (
            solid().faceZonePointDisplacementOld(interfaceI)
        );

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

        if (incrementalResiduals_)
        {
            // Residual calculated based on the displacement increments of the
            // fluid and solid interfaces
            residuals()[interfaceI] =
                solidZonesPointsDispls()[interfaceI]
              - fluidZonesPointsDispls()[interfaceI];
        }
        else
        {
            // Residual based on the current position of the interfaces
            // Potentially this approach will avoid the accumulation of errors
            // that can occur in the incrementalResiduals approach

            // Calculate solid deformed point positions
            vectorField solidZonePointsAtSolid
            (
                solid().globalPatches()[interfaceI].patchPointToGlobal
                (
                    solidMesh().boundaryMesh()
                    [
                        solid().globalPatches()[interfaceI].patch().index()
                    ].localPoints()
                )
              + solidZonePointsDisplsAtSolid
            );

            if (!solid().movingMesh())
            {
                solidZonePointsAtSolid += solidZonePointsTotDisplsAtSolid;
            }

            // Map solid points to the fluid interface

            vectorField solidZonePoints
            (
                fluidZonesPointsDispls()[interfaceI].size(), vector::zero
            );

            interfaceToInterfaceList()[interfaceI].transferPointsZoneToZone
            (
                solidZone,                              // from zone
                fluidZone,                              // to zone
                solidZonePointsAtSolid,                 // from field
                solidZonePoints                         // to field
            );

            // Calculate fluid zone positions
            const vectorField fluidZonePoints
            (
                fluid().globalPatches()[interfaceI].patchPointToGlobal
                (
                    fluidMesh().boundaryMesh()
                    [
                        fluid().globalPatches()[interfaceI].patch().index()
                    ].localPoints()
                )
            );

            residuals()[interfaceI] = solidZonePoints - fluidZonePoints;
        }

        // We will use two definitions of residual
        scalar residualNorm1 = Foam::sqrt(gSum(magSqr(residuals()[interfaceI])));
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

        const vectorField intTotDispl
        (
            interfacesPointsDispls_[interfaceI] + solidZonePointsTotDispl
        );

        const scalar intTotDisplNorm = Foam::sqrt(gSum(magSqr(intTotDispl)));

        if (intTotDisplNorm > maxIntsDisplsNorm_[interfaceI])
        {
            maxIntsDisplsNorm_[interfaceI] = intTotDisplNorm;
        }

        residualNorm2 /= maxIntsDisplsNorm_[interfaceI] + SMALL;

        Info<< "FSI residual2 norm for interface " << interfaceI
            << ": " << residualNorm2 << endl;

        // Legacy solids4foam behavior uses the minimum of the two residual
        // measures, while stricter convergence requires both to be small.
        const scalar residualInterfaceI =
            requireAllResidualMeasures_
          ? max(residualNorm1, residualNorm2)
          : min(residualNorm1, residualNorm2);

        if (requireAllResidualMeasures_)
        {
            Info<< "FSI combined residual norm for interface " << interfaceI
                << " (max): " << residualInterfaceI << endl;
        }

        // Update the maximum residual for all interfaces
        maxResidual = max(maxResidual, residualInterfaceI);
    }

    robinPressureResidual_ = 0;
    robinFluxResidual_ = 0;
    robinKinematicConsistency_ = true;

    if (!hasRobinInterface())
    {
        return maxResidual;
    }

    // Flux scale: the throughput of the physical boundaries plus the flux
    // swept by the moving Robin interfaces
    scalar boundaryFluxNorm = 0;
    scalar interfaceMotionNorm = 0;
    scalar robinInterfaceArea = 0;
    forAll(fluid().phi().boundaryField(), patchI)
    {
        // Physical boundaries only: the flux through the Robin interfaces is
        // the leakage itself, which must not enter its own scale
        bool robinPatch = false;
        forAll(robinInterfaces_, interfaceI)
        {
            robinPatch =
                robinPatch
             || (
                    robinInterfaces_[interfaceI]
                 && fluidPatchIndices()[interfaceI] == patchI
                );
        }

        if (!fluidMesh().boundary()[patchI].coupled() && !robinPatch)
        {
            boundaryFluxNorm +=
                sum(mag(fluid().phi().boundaryField()[patchI]));
        }
    }
    forAll(robinInterfaces_, interfaceI)
    {
        if (robinInterfaces_[interfaceI])
        {
            robinInterfaceArea +=
                sum
                (
                    fluidMesh().magSf().boundaryField()
                    [
                        fluidPatchIndices()[interfaceI]
                    ]
                );
        }
    }
    if (fluidMesh().moving())
    {
        forAll(robinInterfaces_, interfaceI)
        {
            if (robinInterfaces_[interfaceI])
            {
                interfaceMotionNorm +=
                    sum
                    (
                        mag
                        (
                            fluidMesh().phi().boundaryField()
                            [
                                fluidPatchIndices()[interfaceI]
                            ]
                        )
                    );
            }
        }
    }
    reduce(boundaryFluxNorm, sumOp<scalar>());
    reduce(interfaceMotionNorm, sumOp<scalar>());
    reduce(robinInterfaceArea, sumOp<scalar>());
    maxRobinKinematicNorm_ =
        max(maxRobinKinematicNorm_, boundaryFluxNorm + interfaceMotionNorm);

    // Flow at or near rest: without throughput and interface motion the flux
    // scale vanishes, and the leakage floor of the inner solvers would give a
    // large residual, so the optional reference flux bounds the scale
    const scalar robinReferenceFlux =
        robinFluxReferenceVelocity_*robinInterfaceArea;
    const bool fluxScaleLimited = robinReferenceFlux > maxRobinKinematicNorm_;
    const scalar fluxScale = max(maxRobinKinematicNorm_, robinReferenceFlux);

    forAll(robinInterfaces_, interfaceI)
    {
        if (!robinInterfaces_[interfaceI])
        {
            continue;
        }

        const label fluidPatchID = fluidPatchIndices()[interfaceI];
        const scalarField currentPressure
        (
            fluid().solutionP().boundaryField()[fluidPatchID]
        );
        const scalar pressureNorm = Foam::sqrt(gSum(magSqr(currentPressure)));
        maxRobinPressureNorm_[interfaceI] = max
        (
            maxRobinPressureNorm_[interfaceI], pressureNorm
        );

        // Optional reference pressure, as the norm of a uniform interface
        // pressure, bounds the scale of a pressure at or near zero
        const scalar referencePressureNorm =
            robinPressureReference_
           *Foam::sqrt
            (
                scalar(returnReduce(currentPressure.size(), sumOp<label>()))
            );
        const bool pressureScaleLimited =
            referencePressureNorm > maxRobinPressureNorm_[interfaceI];
        const scalar pressureScale =
            max(maxRobinPressureNorm_[interfaceI], referencePressureNorm);

        const scalar pressureResidual =
            Foam::sqrt
            (
                gSum
                (
                    magSqr
                    (
                        currentPressure - robinPressurePrev_[interfaceI]
                    )
                )
            )
           /(pressureScale + SMALL);

        // Interface flux relative to the mesh motion, i.e. the leakage
        // through the moving wall
        const scalarField& leak = fluid().phi().boundaryField()[fluidPatchID];

        // Kinematic residual: the leakage normalised by the throughput plus
        // the interface motion flux. If the Robin condition is kinematically
        // consistent with the mesh motion, the leakage is proportional to
        // the Robin iteration error and vanishes as the iterations converge,
        // as for a Dirichlet velocity condition, and it is used as a
        // convergence criterion. Otherwise the leakage has a discretisation
        // floor (it does not vanish as the iterations converge) and is only
        // reported.
        const scalar leakage = gSum(mag(leak))/(fluxScale + SMALL);

        robinKinematicConsistency_ =
            robinKinematicConsistency_
         && refCast<const elasticWallPressureFvPatchScalarField>
            (
                fluid().solutionP().boundaryField()[fluidPatchID]
            ).kinematicConsistency();

        robinPressurePrev_[interfaceI] = currentPressure;
        robinPressureResidual_ = max
        (
            robinPressureResidual_, pressureResidual
        );
        robinFluxResidual_ = max(robinFluxResidual_, leakage);

        Info<< "Robin pressure residual for interface " << interfaceI
            << ": " << pressureResidual
            << (pressureScaleLimited ? " (reference pressure scale)" : "")
            << nl
            << "Robin leakage-flux residual for interface " << interfaceI
            << ": " << leakage
            << (fluxScaleLimited ? " (reference flux scale)" : "") << endl;
    }

    robinDispHistory_.append(maxResidual);
    robinPressureHistory_.append(robinPressureResidual_);
    robinFluxHistory_.append(robinFluxResidual_);

    evaluateRobinConvergence();

    // The FSI iterations ran out with a Robin residual that has no reference
    // scale above its tolerance: if the flow is at or near rest, the relative
    // residual has no meaningful scale, so point to the reference options
    if
    (
        !robinReferenceHintIssued_
     && coupled()
     && robinConvergenceState_ == 0
     && outerCorr_ >= nOuterCorr_
    )
    {
        const bool pressureUnscaled =
            robinPressureResidual_ > robinPressureTolerance_
         && robinPressureReference_ <= 0;
        const bool fluxUnscaled =
            robinKinematicConsistency_
         && robinFluxResidual_ > robinFluxTolerance_
         && robinFluxReferenceVelocity_ <= 0;

        if (pressureUnscaled || fluxUnscaled)
        {
            robinReferenceHintIssued_ = true;

            Info<< nl << "Robin FSI iterations reached nOuterCorr with a"
                << " Robin residual above its tolerance. If the flow is at or"
                << " near rest (little throughput or interface motion), the"
                << " relative residual has no physical scale: set";
            if (fluxUnscaled)
            {
                Info<< " robinFluxReferenceVelocity (a characteristic"
                    << " velocity of the flow)";
            }
            if (pressureUnscaled)
            {
                Info<< (fluxUnscaled ? " and" : "")
                    << " robinPressureReference (a characteristic pressure)";
            }
            Info<< " in fsiProperties. Otherwise, check the inner-solver"
                << " tolerances and nOuterCorr. This hint is written once."
                << nl << endl;
        }
    }

    return maxResidual;
}


bool Foam::fluidSolidInterface::hasRobinInterface() const
{
    forAll(robinInterfaces_, interfaceI)
    {
        if (robinInterfaces_[interfaceI])
        {
            return true;
        }
    }

    return false;
}


namespace Foam
{
    // Convergence test of one residual history of the current time step.
    //
    // For a linearly converging fixed-point iteration with contraction rate
    // rho, the remaining error after the latest change R is at most
    // R*rho/(1 - rho). The rate is the largest of the last two ratios of
    // successive residuals, using only the iterations in which the solid
    // solution changed (update[i]): when the solid solver skips its solve
    // (its own tolerance is met) the interface data do not change, and the
    // resulting plateaus carry no contraction information.
    //
    // The residual is converged if the error estimate satisfies the
    // tolerance, or if the residual itself does and the iterations are not
    // contracting slowly (rate below slowRate); slowly contracting
    // iterations must satisfy the error estimate. It is stalled if the
    // geometric-mean rate of the last nStall ratios is at or above
    // stallRate.
    static void robinResidualTest
    (
        const UList<scalar>& history,
        const UList<bool>& update,
        const scalar tolerance,
        const scalar slowRate,
        const scalar stallRate,
        const label nStall,
        bool& converged,
        bool& stalled,
        scalar& error,
        scalar& rate
    )
    {
        converged = false;
        stalled = false;
        rate = -1;
        error = GREAT;

        if (history.empty())
        {
            return;
        }

        scalar maxValue = 0;
        forAll(history, i)
        {
            maxValue = max(maxValue, history[i]);
        }
        const scalar zeroTol = 1e-12*maxValue + VSMALL;

        const scalar R = history.last() > zeroTol ? history.last() : 0;

        // Residuals of the iterations in which the solid solution changed
        DynamicList<scalar> values;
        forAll(history, i)
        {
            if (update[i] && history[i] > zeroTol)
            {
                values.append(history[i]);
            }
        }

        const label n = values.size();

        if (n < 3)
        {
            // Too few iterations to estimate the rate
            error = R;
            converged = R <= tolerance;
            return;
        }

        rate = max(values[n - 1]/values[n - 2], values[n - 2]/values[n - 3]);

        if (nStall > 0 && n > nStall)
        {
            scalar logRate = 0;
            for (label i = n - nStall; i < n; i++)
            {
                logRate += Foam::log(values[i]/values[i - 1]);
            }
            // A stall is a plateau: the mean rate is close to one from both
            // sides, so growing (diverging) residuals are not accepted
            const scalar meanRate = Foam::exp(logRate/nStall);
            stalled = meanRate >= stallRate && meanRate <= 1/stallRate;
        }

        const scalar rho = min(rate, 0.99);
        error = R*rho/(1 - rho);

        converged =
            error <= tolerance
         || (R <= tolerance && rate < slowRate);
    }
}


void Foam::fluidSolidInterface::evaluateRobinConvergence()
{
    robinConvergenceState_ = 0;

    // Before the coupling starts only the displacement residual applies
    if (!coupled())
    {
        robinConvergenceState_ =
            (
                robinDispHistory_.size()
             && robinDispHistory_.last() <= outerCorrTolerance_
            ) ? 1 : 0;

        return;
    }

    // Iterations in which the solid solution changed
    scalar maxDisp = 0;
    forAll(robinDispHistory_, i)
    {
        maxDisp = max(maxDisp, robinDispHistory_[i]);
    }
    List<bool> update(robinDispHistory_.size());
    forAll(update, i)
    {
        update[i] = robinDispHistory_[i] > 1e-12*maxDisp + VSMALL;
    }

    if (robinConvergence_ == "residual")
    {
        const bool converged =
            robinDispHistory_.size()
         && robinDispHistory_.last() <= outerCorrTolerance_
         && robinPressureResidual_ <= robinPressureTolerance_
         && (
                !robinKinematicConsistency_
             || robinFluxResidual_ <= robinFluxTolerance_
            );

        if (converged)
        {
            robinConvergenceState_ = 1;
            return;
        }

        // Stalled residuals: they no longer decrease (e.g. inner-solver
        // tolerances), so further iterations are pointless; they are
        // accepted if they are close to their tolerances. Only the stall
        // flag of robinResidualTest is used here.
        const UList<scalar>* histories[3] =
        {
            &robinDispHistory_, &robinPressureHistory_, &robinFluxHistory_
        };
        const scalar tolerances[3] =
        {
            outerCorrTolerance_, robinPressureTolerance_, robinFluxTolerance_
        };
        const label nMeasures = robinKinematicConsistency_ ? 3 : 2;

        bool allAcceptable = true;
        for (label i = 0; i < nMeasures; i++)
        {
            bool measureConverged = false;
            bool stalled = false;
            scalar error = GREAT;
            scalar rate = -1;

            robinResidualTest
            (
                *histories[i],
                update,
                tolerances[i],
                robinSlowRate_,
                robinStallRate_,
                robinStallIterations_,
                measureConverged,
                stalled,
                error,
                rate
            );

            const scalar R = histories[i]->last();

            allAcceptable =
                allAcceptable
             && (
                    R <= tolerances[i]
                 || (stalled && R <= robinStallToleranceFactor_*tolerances[i])
                );
        }

        if (allAcceptable)
        {
            robinConvergenceState_ = 2;
        }

        return;
    }

    const UList<scalar>* histories[3] =
    {
        &robinDispHistory_, &robinPressureHistory_, &robinFluxHistory_
    };
    const scalar tolerances[3] =
    {
        outerCorrTolerance_, robinPressureTolerance_, robinFluxTolerance_
    };
    const char* names[3] = {"displacement", "pressure", "kinematic"};

    bool allConverged = true;
    bool allAcceptable = true;

    // The kinematic residual is only a convergence measure if the Robin
    // condition is kinematically consistent with the mesh motion
    const label nMeasures = robinKinematicConsistency_ ? 3 : 2;

    Info<< "Robin iteration error estimates:";

    for (label i = 0; i < nMeasures; i++)
    {
        bool converged = false;
        bool stalled = false;
        scalar error = GREAT;
        scalar rate = -1;

        robinResidualTest
        (
            *histories[i],
            update,
            tolerances[i],
            robinSlowRate_,
            robinStallRate_,
            robinStallIterations_,
            converged,
            stalled,
            error,
            rate
        );

        // A stalled residual cannot be reduced further by FSI iterations
        // (e.g. inner-solver tolerances); it is accepted if it is close to
        // its tolerance
        const bool acceptable =
            converged
         || (
                stalled
             && histories[i]->last()
             <= robinStallToleranceFactor_*tolerances[i]
            );

        allConverged = allConverged && converged;
        allAcceptable = allAcceptable && acceptable;

        Info<< " " << names[i] << " " << histories[i]->last()
            << " (error " << error << ", rate " << rate
            << (stalled ? ", stalled)" : ")");
    }

    Info<< endl;

    if (allConverged)
    {
        robinConvergenceState_ = 1;
    }
    else if (allAcceptable)
    {
        robinConvergenceState_ = 2;
    }
}


bool Foam::fluidSolidInterface::couplingConverged
(
    const scalar displacementResidual
) const
{
    // Before the coupling starts (couplingStartTime) the fluid is solved
    // without solid feedback, so only the displacement residual applies
    if (!hasRobinInterface() || !coupled())
    {
        return displacementResidual <= outerCorrTolerance_;
    }

    // At least two FSI iterations for the iteration-error estimate, so that
    // the Robin condition has been updated with the solid response of the
    // current time step
    if (robinConvergence_ == "iterationError" && outerCorr_ < 2)
    {
        return false;
    }

    if (robinConvergenceState_ == 2)
    {
        Info<< "Robin FSI iterations stalled after " << outerCorr_
            << " iterations; accepted because the residuals are within "
            << robinStallToleranceFactor_ << " times the tolerances" << endl;
    }

    return robinConvergenceState_ > 0;
}


void Foam::fluidSolidInterface::writeResidualLine(const scalar residualNorm)
{
    if (!writeResidualsToFile() || !Pstream::master())
    {
        return;
    }

    residualFile()
        << runTime().value() << " "
        << outerCorr() << " "
        << residualNorm;

    if (hasRobinInterface())
    {
        // Before the coupling starts the Robin residuals are not evaluated
        // and a step is accepted on the displacement residual alone
        const label state =
            coupled()
          ? robinConvergenceState_
          : label(residualNorm <= outerCorrTolerance_ ? 1 : 0);

        residualFile()
            << " " << robinPressureResidual_
            << " " << robinFluxResidual_
            << " " << state;
    }

    residualFile() << endl;
}


void Foam::fluidSolidInterface::updateMovingWallPressureAcceleration()
{
    forAll(fluid().globalPatches(), interfaceI)
    {
        if
        (
            isA<movingWallPressureFvPatchScalarField>
            (
                fluid().solutionP().boundaryField()
                [
                    fluidPatchIndices()[interfaceI]
                ]
            )
        )
        {
            Info<< "Setting acceleration at fluid side of the interface: "
                << fluidMesh().boundary()
                   [
                       fluidPatchIndices()[interfaceI]
                   ].name()
                << endl;

            // Take references to zones
            const standAlonePatch& fluidZone =
               fluid().globalPatches()[interfaceI].globalPatch();
            const standAlonePatch& solidZone =
               solid().globalPatches()[interfaceI].globalPatch();

            const vectorField solidZoneAcceleration
            (
               solid().faceZoneAcceleration(interfaceI)
            );

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

            const vectorField fluidPatchAcceleration
            (
                fluid().globalPatches()
                [
                   interfaceI
                ].globalFaceToPatch(fluidZoneAcceleration)
            );

            vectorField& prevAcceleration =
                const_cast<movingWallPressureFvPatchScalarField&>
                (
                    refCast<const movingWallPressureFvPatchScalarField>
                    (
                        fluid().solutionP().boundaryField()
                        [
                            fluidPatchIndices()[interfaceI]
                        ]
                    )
                ).prevAcceleration();

            if (coupled())
            {
                prevAcceleration = fluidPatchAcceleration;
            }
            else
            {
                prevAcceleration = vector::zero;
            }
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
                fluid().solutionP().boundaryField()
                [
                    fluidPatchIndices()[interfaceI]
                ]
            )
        )
        {
            Info<< "Setting acceleration and previous pressure at fluid side of "
                << "the interface: "
                << fluidMesh().boundary()
                   [
                       fluidPatchIndices()[interfaceI]
                   ].name()
                << endl;

            // Take references to zones
            const standAlonePatch& fluidZone =
               fluid().globalPatches()[interfaceI].globalPatch();
            const standAlonePatch& solidZone =
               solid().globalPatches()[interfaceI].globalPatch();

            const vectorField solidZoneAcceleration
            (
               solid().faceZoneAcceleration(interfaceI)
            );

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

            const vectorField fluidPatchAcceleration
            (
                fluid().globalPatches()
                [
                    interfaceI
                ].globalFaceToPatch(fluidZoneAcceleration)
            );

            vectorField& prevAcceleration =
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().solutionP().boundaryField()
                        [
                            fluidPatchIndices()[interfaceI]
                        ]
                    )
                ).prevAcceleration();

            scalarField& prevPressure =
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().solutionP().boundaryField()
                        [
                            fluidPatchIndices()[interfaceI]
                        ]
                    )
                ).prevPressure();

            if (coupled())
            {
                prevAcceleration = fluidPatchAcceleration;
                prevPressure =
                    fluid().patchSolutionPressureForce
                    (
                        fluidPatchIndices()[interfaceI]
                    );
            }
            else
            {
                // ZT: Helps to improve stability in case of
                // uncoupled simulation where acceleration
                // shoud be exactly zero.
                prevAcceleration = vector::zero;
                // ZT: Pressure is not zero.
                // prevPressure = 0;
                prevPressure =
                    fluid().patchSolutionPressureForce
                    (
                        fluidPatchIndices()[interfaceI]
                    );
            }

            elasticWallPressureFvPatchScalarField& pRobin =
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().solutionP().boundaryField()
                        [
                            fluidPatchIndices()[interfaceI]
                        ]
                    )
                );

            // Update the Robin coefficient estimate from the new pressure and
            // solid acceleration pair
            pRobin.updateRobinCoefficient();

            // The fluid model may only make the interface flux consistent
            // with the mesh motion if the mesh follows the solid
            pRobin.setFluidMeshFollowsSolid(fluidMeshFollowsSolid());
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
            reduce(fluidZonesPointsDispls[interfaceI], FieldSumOp<vector>());

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

void Foam::fluidSolidInterface::end()
{
    this->IOobject::rename(this->IOobject::name()+".withDefaultValues");
    this->regIOobject::write();
    solid().end();
    fluid().end();
}

// ************************************************************************* //
