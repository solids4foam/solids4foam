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

#include "newtonIcoFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "findRefCell.H"
#include "fixedValueFvPatchFields.H"
#include "compatibilityFunctions.H"


namespace
{

void scaleFvScalarMatrix(Foam::fvScalarMatrix& matrix, const Foam::scalar scale)
{
    matrix.diag() *= scale;

    if (matrix.hasUpper())
    {
        matrix.upper() *= scale;
    }

    if (matrix.hasLower())
    {
        matrix.lower() *= scale;
    }

    Foam::FieldField<Foam::Field, Foam::scalar>& internalCoeffs =
        matrix.internalCoeffs();
    Foam::FieldField<Foam::Field, Foam::scalar>& boundaryCoeffs =
        matrix.boundaryCoeffs();

    forAll(internalCoeffs, patchI)
    {
        internalCoeffs[patchI] *= scale;
        boundaryCoeffs[patchI] *= scale;
    }
}


const Foam::dictionary& solverControls
(
    const Foam::fvMesh& mesh,
    const Foam::word& fieldName
)
{
#ifdef OPENFOAM_NOT_EXTEND
    return mesh.solverDict(fieldName);
#else
    return mesh.solutionDict().solver(fieldName);
#endif
}


class retryTimeStateBuilder
:
    public Foam::TimeState
{
public:

    static Foam::TimeState rollbackState(const Foam::Time& runTime)
    {
        retryTimeStateBuilder state;

        static_cast<Foam::TimeState&>(state) =
            static_cast<const Foam::TimeState&>(runTime);

        // Restore the previous successful time-step length as the "saved"
        // value so the retry does not treat the failed step as the last
        // accepted one when operator++() updates deltaT0.
        state.deltaTSave_ = runTime.deltaT0Value();
#ifdef OPENFOAM_NOT_EXTEND
        state.writeTime_ = false;
#endif

        return state;
    }
};

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newtonIcoFluid, 0);
addToRunTimeSelectionTable(fluidModel, newtonIcoFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //


void newtonIcoFluid::makeRAUf() const
{
    if (rAUfPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rAUfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rAUf",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimTime, 0.0)
        )
    );
}


const surfaceScalarField& newtonIcoFluid::rAUf() const
{
    if (rAUfPtr_.empty())
    {
        makeRAUf();
    }

    return autoPtrRef(rAUfPtr_);
}


surfaceScalarField& newtonIcoFluid::rAUf()
{
    if (rAUfPtr_.empty())
    {
        makeRAUf();
    }

    return autoPtrRef(rAUfPtr_);
}


void newtonIcoFluid::restoreOldTimeState
(
    const pointField& oldPoints,
    const bool meshMoved
)
{
    dynamicFvMesh& mesh = this->mesh();
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    surfaceScalarField& phi = this->phi();

    U = U.oldTime();
    U.correctBoundaryConditions();

    p = p.oldTime();
    p.correctBoundaryConditions();

    if (meshMoved)
    {
        mesh.movePoints(oldPoints);
    }

    phi = fvc::interpolate(U) & mesh.Sf();

    if (meshMoved)
    {
        fvc::makeRelative(phi, U);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonIcoFluid::newtonIcoFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "Up",
        fileName
        (
            fluidProperties().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solutionLocation::CELLS,
        fluidProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true
    ),
    Uf_(),
    rAUfPtr_(),
    pRefCell_(-1),
    pRefValue_(0.0),
    laminarTransport_(U(), phi()),
    turbulence_
    (
#ifdef OPENFOAM_ORG
        incompressible::momentumTransportModel::New
#else
        incompressible::turbulenceModel::New
#endif
        (
            U(), phi(), laminarTransport_
        )
    ),
    rho_(laminarTransport_.lookup("rho")),
    momentumStabilisationPtr_(),
    pressureStabilisationPtr_(),
    pressureScaleFactor_
    (
        fluidProperties().lookupOrDefault<scalar>("pressureScaleFactor", 1.0)
    ),
    scaleMixedPetScFields_
    (
        fluidProperties().lookupOrDefault<Switch>
        (
            "scaleMixedPetScFields", false
        )
    ),
    pressureUnknownScaleType_
    (
        fluidProperties().lookupOrDefault<word>
        (
            "pressureUnknownScale", "none"
        )
    ),
    pressureUnknownScale_(1.0),
    pressureEqnScale_(pressureScaleFactor_),
    blockSize_(fluidModel::twoD() ? 3 : 4),
    tsLogPtr_()
{
    setRefCell(p(), fluidProperties(), pRefCell_, pRefValue_);
    //mesh().setFluxRequired(p().name());

    dictionary defaultMomentumStabSubDict;
    defaultMomentumStabSubDict.add("type", "diffStencilLaplacian");
    defaultMomentumStabSubDict.add("scaleFactor", 0.0);

    dictionary defaultPressureStabSubDict;
    defaultPressureStabSubDict.add("type", "RhieChow");
    defaultPressureStabSubDict.add("scaleFactor", 1.0);

    if (!fluidProperties().found("stabilisation"))
    {
        dictionary stabDict;
        stabDict.add("momentum", defaultMomentumStabSubDict);
        stabDict.add("pressure", defaultPressureStabSubDict);
        fluidProperties().add("stabilisation", stabDict);
    }

    dictionary& stabDict = fluidProperties().subDict("stabilisation");

    if
    (
        (stabDict.found("type") || stabDict.found("scaleFactor"))
     && !stabDict.found("pressure")
    )
    {
        Info<< "Using legacy fluid stabilisation format as pressure "
            << "stabilisation" << endl;

        dictionary legacyPressureStabSubDict;

        if (stabDict.found("type"))
        {
            legacyPressureStabSubDict.add("type", word(stabDict.lookup("type")));
        }
        else
        {
            legacyPressureStabSubDict.add("type", word("RhieChow"));
        }

        if (stabDict.found("scaleFactor"))
        {
            legacyPressureStabSubDict.add
            (
                "scaleFactor",
                readScalar(stabDict.lookup("scaleFactor"))
            );
        }
        else
        {
            legacyPressureStabSubDict.add("scaleFactor", 1.0);
        }

        if (stabDict.found("omega"))
        {
            legacyPressureStabSubDict.add
            (
                "omega",
                dimensionedScalar(stabDict.lookup("omega"))
            );
        }

        if (stabDict.found("innerScaleFactor"))
        {
            legacyPressureStabSubDict.add
            (
                "innerScaleFactor",
                dimensionedScalar(stabDict.lookup("innerScaleFactor"))
            );
        }

        if (stabDict.found("outerScaleFactor"))
        {
            legacyPressureStabSubDict.add
            (
                "outerScaleFactor",
                dimensionedScalar(stabDict.lookup("outerScaleFactor"))
            );
        }

        stabDict.add("pressure", legacyPressureStabSubDict);
    }

    if (!stabDict.found("momentum"))
    {
        stabDict.add("momentum", defaultMomentumStabSubDict);
    }

    if (!stabDict.found("pressure"))
    {
        stabDict.add("pressure", defaultPressureStabSubDict);
    }

    momentumStabilisationPtr_ =
        stabilisationModel::New
        (
            mesh(),
            stabDict.subDict("momentum"),
            dimVelocity/dimLength
        );

    pressureStabilisationPtr_ =
        stabilisationModel::New
        (
            mesh(),
            stabDict.subDict("pressure"),
            p().dimensions()/dimLength
        );

#ifdef OPENFOAM_NOT_EXTEND
    turbulence_->validate();
#endif

    U().oldTime().oldTime();

    // if (mesh().dynamic())
    // {
    //     Info<< "Constructing face velocity Uf\n" << endl;

        // Uf_.reset
        // (
        //     new surfaceVectorField
        //     (
        //         IOobject
        //         (
        //             "Uf",
        //             runTime.timeName(),
        //             mesh(),
        //             IOobject::READ_IF_PRESENT,
        //             IOobject::AUTO_WRITE
        //         ),
        //         fvc::interpolate(U())
        //     )
        // );

        // Uf_().oldTime();

        // if
        // (
        //     word(mesh().ddtScheme("ddt(" + U().name() +')'))
        //  == fv::backwardDdtScheme<vector>::typeName
        // )
        // {
        //     Uf_().oldTime().oldTime();
        // }
    // }

    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& phi = this->phi();
    #include "CourantNo.H"

    if (mag(pressureScaleFactor_ - 1.0) > SMALL)
    {
        Info<< "pressureScaleFactor = " << pressureScaleFactor_ << endl;
    }

    // Resolve the PETSc pressure unknown scale.
    // Default ("none") leaves pressureUnknownScale_ = 1, which is a
    // no-op everywhere pack/unpackSolution is used
    pressureUnknownScale_ = 1.0;
    if (scaleMixedPetScFields_)
    {
        if (pressureUnknownScaleType_ == "none")
        {
            pressureUnknownScale_ = 1.0;
        }
        else if
        (
            pressureUnknownScaleType_ == "user"
         || pressureUnknownScaleType_ == "scalar"
        )
        {
            pressureUnknownScale_ =
                readScalar
                (
                    fluidProperties().lookup("pressureUnknownScaleValue")
                );
        }
        else if (pressureUnknownScaleType_ == "dynamicHead")
        {
            // 0.5*max(|U|)^2 evaluated over the initial U field. This is
            // the natural magnitude of kinematic pressure (p/rho) in
            // regions of significant pressure gradient for an
            // incompressible flow. It is fixed at construction so that
            // MFFD perturbations within a Newton solve see a stable
            // scale
            scalar Umax2 = 0;
            forAll(U(), cellI)
            {
                Umax2 = max(Umax2, magSqr(U()[cellI]));
            }
            reduce(Umax2, maxOp<scalar>());

            if (Umax2 <= VSMALL)
            {
                WarningInFunction
                    << "pressureUnknownScale = dynamicHead requested but "
                    << "the initial U field is zero everywhere. Falling "
                    << "back to pressureUnknownScale_ = 1. Set "
                    << "pressureUnknownScale to 'user' with an explicit "
                    << "pressureUnknownScaleValue if the case starts "
                    << "from rest."
                    << endl;
                pressureUnknownScale_ = 1.0;
            }
            else
            {
                pressureUnknownScale_ = 0.5*Umax2;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unknown pressureUnknownScale "
                << pressureUnknownScaleType_ << nl
                << "Valid options are user, scalar, dynamicHead, none"
                << abort(FatalError);
        }

        if (pressureUnknownScale_ <= VSMALL)
        {
            FatalErrorInFunction
                << "pressureUnknownScale must be positive, found "
                << pressureUnknownScale_ << abort(FatalError);
        }
    }

    if (mag(pressureUnknownScale_ - 1.0) > SMALL)
    {
        Info<< "PETSc pressure unknown scale = " << pressureUnknownScale_
            << " (scaleMixedPetScFields = " << scaleMixedPetScFields_
            << ", pressureUnknownScale = " << pressureUnknownScaleType_
            << ")" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> newtonIcoFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tmpRef(tvF) = rho_.value()
       *(
            mesh().boundary()[patchID].nf()
#ifdef OPENFOAM_ORG
          & (-turbulence_->devTau()().boundaryField()[patchID])
#else
          & (-turbulence_->devReff()().boundaryField()[patchID])
#endif
        );

    return tvF;
}


tmp<scalarField> newtonIcoFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tmpRef(tpF) = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<vectorField> newtonIcoFluid::patchViscousForce
(
    const label patchID, const solidModel& motion
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    const tensorField Fm(I + motion.gradD().T());
    const scalarField Jm(det(Fm));
    const tensorField invFm(inv(Fm));
    const scalarField nuEff(turbulence_->nuEff()().boundaryField()[patchID]);
    const vectorField& Sf = mesh().boundary()[patchID].Sf();
    const vectorField deformedSf(Jm*invFm.T() & Sf);
    const vectorField deformedNf(deformedSf/mag(deformedSf));
    const tensorField& gradU = this->gradU().boundaryField()[patchID];

    tmpRef(tvF) = rho_.value()*deformedNf & (nuEff*invFm.T() & gradU);

    // Deformed mesh
    // tmpRef(tvF) = rho_.value()
    //    *(
    //         mesh().boundary()[patchID].nf()
    //       & (-turbulence_->devReff()().boundaryField()[patchID])
    //     );

    return tvF;
}


void newtonIcoFluid::setDeltaT(Time& runTime)
{
#ifdef USE_PETSC
    if
    (
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
     && foamPetscSnesHelper::snesHasRun()
    )
    {
        const scalar maxDeltaT =
            readScalar(runTime.controlDict().lookup("maxDeltaT"));
        const scalar minDeltaT =
            readScalar(runTime.controlDict().lookup("minDeltaT"));

        const int minTargetNIter =
            runTime.controlDict().lookupOrDefault<int>("minTargetNIter", 3);
        const int maxTargetNIter =
            runTime.controlDict().lookupOrDefault<int>("maxTargetNIter", 6);

        const Switch enableTimeStepLog =
            runTime.controlDict().lookupOrDefault("logTimeStepAdjustments", true);

        PetscInt numIter;
        SNESGetIterationNumber(foamPetscSnesHelper::snes(), &numIter);

        SNESConvergedReason reason;
        SNESGetConvergedReason(foamPetscSnesHelper::snes(), &reason);

        const scalar currentDeltaT = runTime.deltaTValue();
        scalar newDeltaT = currentDeltaT;

        // if (reason == SNES_DIVERGED_FUNCTION_DOMAIN)
        if (reason < 0)
        {
            // SNES failed to converge
            newDeltaT = max(0.25*currentDeltaT, minDeltaT);
            Info<< nl << "SNES failed to converge: "
                << "reducing timestep to " << newDeltaT << endl;
        }
        else
        {
            // Guard against zero
            if (numIter <= 0)
            {
                numIter = 1;
            }

            scalar factor = 1.0;

            if (numIter > maxTargetNIter)
            {
                factor = max(0.5, 0.9*scalar(maxTargetNIter)/numIter);
            }
            else if (numIter < minTargetNIter)
            {
                factor = min(1.5, 1.1*scalar(maxTargetNIter)/numIter);
            }

            newDeltaT = min(max(factor*currentDeltaT, minDeltaT), maxDeltaT);
        }

        Info<< "Nonlinear iterations = " << numIter << nl
            << "Old time step        = " << currentDeltaT << nl
            << "New time step        = " << newDeltaT << nl << endl;

        runTime.setDeltaT(newDeltaT);

        if (enableTimeStepLog)
        {
            if (tsLogPtr_.empty())
            {
                const fileName timeStepLogFile =
                    runTime.controlDict().lookupOrDefault<fileName>
                    (
                        "timeStepLogFile", "timeStepLog.dat"
                    );

                tsLogPtr_.set(new OFstream(timeStepLogFile));

                tsLogPtr_()
                    << "Time currentDeltaT newDeltaT numIter reason" << endl;
            }

            tsLogPtr_()
                << runTime.timeName() << " "
                << currentDeltaT << " "
                << newDeltaT << " "
                << numIter << " "
                << reason << endl;
        }
    }
#else
    if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        static bool warned = false;

        if (!warned)
        {
            WarningInFunction
                << "Ignoring adjustTimeStep because PETSc support is not "
                << "enabled" << endl;

            warned = true;
        }
    }
#endif
}


bool newtonIcoFluid::evolve()
{
#ifdef USE_PETSC
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    // const Time& runTime = fluidModel::runTime();
    Time& time = physicsModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    volVectorField& U = this->U();
    surfaceScalarField& phi = this->phi();
    // autoPtr<surfaceVectorField>& Uf = Uf_;
    //scalar& cumulativeContErr = cumulativeContErr_;
    //const bool correctPhi = correctPhi_;
    //const bool checkMeshCourantNo = checkMeshCourantNo_;
    //const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;

    // Solution predictor
    const Switch predictor
    (
        fluidProperties().lookupOrDefault<Switch>("predictor", false)
    );

    const Switch adjustTimeStep
    (
        time.controlDict().lookupOrDefault("adjustTimeStep", false)
    );

    const label maxTimeStepRetries
    (
        fluidProperties().lookupOrDefault<label>("maxTimeStepRetries", 10)
    );

    label timeStepRetry = 0;
    bool retriedCurrentDeltaT = false;

    while (true)
    {
        // Update U boundary conditions
        U.correctBoundaryConditions();

        {
            const Time& runTime = mesh.time();
            #include "CourantNo.H"
        }

        const scalar failedTimeValue = time.value();
        const scalar failedDeltaT = time.deltaTValue();
        const label oldTimeIndex = time.timeIndex() - 1;
        const scalar oldTimeValue = failedTimeValue - failedDeltaT;
        const pointField oldPoints(mesh.points());
        const TimeState retryTimeState =
            retryTimeStateBuilder::rollbackState(time);

        if (predictor && time.timeIndex() > 1) // && newTimeStep())
        {
            Info<< "Applying a linear predictor to velocity" << endl;
            U = 2.0*U.oldTime() - U.oldTime().oldTime();
            Info<< "Applying a linear predictor to velocity: done" << endl;
        }

        // Seed the PETSc Vec from the current U and p fields on every
        // timestep, including the first one. Otherwise the very first
        // formResidual reads an unset (zero) PETSc Vec, overwrites the
        // user-supplied initial U and p, and SNES starts from an
        // internal-vs-boundary-inconsistent state. This is especially
        // damaging when a fixedValue pressure boundary differs from
        // zero, because grad(p) then sees a huge boundary step that
        // produces a spurious momentum residual on the first iteration
        packSolution(foamPetscSnesHelper::solution());

        // Update the mesh, unless the FSI interface already moved it.
        if (fluidModel::fsiMeshUpdate())
        {
            fluidModel::fsiMeshUpdateChanged();
        }
        else
        {
#ifdef OPENFOAM_COM
            mesh.controlledUpdate();
#else
            mesh.update();
#endif
        }

        const bool meshMoved = mesh.changing();

        // Update the flux
        phi = fvc::interpolate(U) & mesh.Sf();

        // If the mesh moved, update the flux and make it relative to the mesh
        // motion
        if (meshMoved)
        {
            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        // Keep the pre-solve state so a failed SNES solve can be retried.
        foamPetscSnesHelper::storeSolutionBackup();

        // Solve the nonlinear system and check the convergence
        Info<< "Solving the fluid for U and p" << endl;
        const int solveStatus = foamPetscSnesHelper::solve(true);

        if (solveStatus >= 0)
        {
            break;
        }

        VecCopy
        (
            foamPetscSnesHelper::solutionBackup(),
            foamPetscSnesHelper::solution()
        );

        restoreOldTimeState(oldPoints, meshMoved);

        static_cast<TimeState&>(time) = retryTimeState;
        time.setTime(oldTimeValue, oldTimeIndex);

        if (!retriedCurrentDeltaT)
        {
            retriedCurrentDeltaT = true;

            foamPetscSnesHelper::resetSnesSolverState();

            ++time;

            Info<< "Retrying the failed PETSc time step at unchanged deltaT = "
                << time.deltaTValue() << " after resetting PETSc solver state"
                << " at Time = " << time.timeName() << nl << endl;

            continue;
        }

        if (!adjustTimeStep)
        {
            FatalErrorInFunction
                << "PETSc SNES failed to converge and the previous time-step "
                << "state has been restored, and a same-deltaT PETSc reset "
                << "retry has already been attempted, but `adjustTimeStep` is "
                << "disabled."
                << nl << "Enable `adjustTimeStep` to retry the failed time "
                << "step with a reduced deltaT."
                << abort(FatalError);
        }

        ++timeStepRetry;

        if (timeStepRetry > maxTimeStepRetries)
        {
            FatalErrorInFunction
                << "Exceeded the maximum number of failed PETSc retries ("
                << maxTimeStepRetries << ") at time " << failedTimeValue
                << " with deltaT = " << failedDeltaT << nl
                << "Set a larger `maxTimeStepRetries` if you would like more "
                << "recovery attempts."
                << abort(FatalError);
        }

        setDeltaT(time);

        if (time.deltaTValue() >= failedDeltaT*(1.0 - SMALL))
        {
            FatalErrorInFunction
                << "PETSc SNES failed to converge at the minimum allowed time "
                << "step. The old-time state has been restored, but deltaT "
                << "could not be reduced below " << failedDeltaT
                << " for a retry."
                << abort(FatalError);
        }

        foamPetscSnesHelper::resetSnesSolverState();

        ++time;

        Info<< "Retrying the failed PETSc time step with deltaT = "
            << time.deltaTValue() << " at Time = "
            << time.timeName() << nl << endl;
    }

    // Retrieve the solution: map the PETSc Vec back into the U field
    // and the p field. Pressure unknown scaling (pHat -> p) and the
    // boundary-condition corrections are applied by unpackSolution
    unpackSolution(foamPetscSnesHelper::solution());

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi);

    // Update the flux
    //phi = mesh.Sf() & Uf();
    phi = mesh.Sf() & fvc::interpolate(U);

    if (mesh.changing())
    {
        // Enable: needed for inlet/outlet?
        // #include "correctPhi.esi.H"

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
    }

    // Correct transport and turbulence models once per call to evolve
    Info<< nl << "Correcting the transport model" << endl;
    laminarTransport_.correct();

    Info<< nl << "Correcting the turbulence model" << endl;
    turbulence_->correct();

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


void newtonIcoFluid::clearRAUf()
{
    rAUfPtr_.clear();
}


#ifdef USE_PETSC

void newtonIcoFluid::unpackSolution(const Vec x)
{
    volVectorField& U = const_cast<volVectorField&>(this->U());
    volScalarField& p = const_cast<volScalarField&>(this->p());

    // Copy x into the U field
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0, // Location of U
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );
    U.correctBoundaryConditions();

    // Copy the scaled pressure unknown pHat from x into the physical
    // p field via p = pressureUnknownScale_ * pHat. When the scale is
    // 1 (default) this reduces to a direct extract
    scalarField& pI = p;
    if (mag(pressureUnknownScale_ - 1.0) > SMALL)
    {
        scalarField pHat(pI.size(), 0.0);
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pHat, blockSize_ - 1
        );
        pI = pressureUnknownScale_*pHat;
    }
    else
    {
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );
    }
    p.correctBoundaryConditions();
}


void newtonIcoFluid::packSolution(Vec x)
{
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        U(),
        x,
        0, // Location of U
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Insert the scaled pressure unknown pHat = p/pressureUnknownScale_.
    // When the scale is 1 (default) this reduces to a direct insert
    if (mag(pressureUnknownScale_ - 1.0) > SMALL)
    {
        scalarField pHat(p());
        pHat /= pressureUnknownScale_;
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pHat, x, blockSize_ - 1
        );
    }
    else
    {
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            p(), x, blockSize_ - 1
        );
    }
}


label newtonIcoFluid::initialiseJacobian(Mat& jac)
{
    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label newtonIcoFluid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label newtonIcoFluid::formResidual
(
    Vec f,         // Residual
    const Vec x,   // Solution
    const bool extrapolatedFlux
)
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    // Take references
    dynamicFvMesh& mesh = this->mesh();

    // Copy x into the U and p fields (pHat -> p when scale != 1).
    // unpackSolution also corrects the boundary conditions on both
    // fields, so the rest of this routine can use them directly
    unpackSolution(x);

    volVectorField& U = const_cast<volVectorField&>(this->U());
    volScalarField& p = const_cast<volScalarField&>(this->p());
    surfaceScalarField& phi = const_cast<surfaceScalarField&>(this->phi());

    gradU() = fvc::grad(U);

    // Lookup the forceImplicitFlux flag
    const Switch forceImplicitFlux =
        fluidProperties().lookupOrDefault<Switch>("forceImplicitFlux", false);

    if (forceImplicitFlux || !extrapolatedFlux)
    {
        phi = fvc::interpolate(U) & mesh.Sf();
    }
    else
    {
        if
        (
            Switch
            (
                fluidProperties().lookup("fluidFluxExtrapolationAlgorithm1")
            )
        )
        {
            // Equation 6.10
            phi = fvc::interpolate
                  (
                      2.0*U.oldTime() - U.oldTime().oldTime()
                  ) & mesh.Sf();
        }
        else
        {
            // Equation 6.30
            phi = fvc::interpolate
                  (
                      2.25*U.oldTime()
                    - 1.5*U.oldTime().oldTime()
                    + 0.25*U.oldTime().oldTime().oldTime()
                  ) & mesh.Sf();
        }
    }

    // Absolute flux
    const surfaceScalarField phiAbs("phiAbs", phi);

    if (mesh.changing())
    {
        fvc::makeRelative(phi, U);

        forAll(U.boundaryField(), patchI)
        {
            if (mesh.boundaryMesh()[patchI].type() == "wall")
            {
                boundaryFieldRef(phi)[patchI] = 0.0;
            }
        }
    }

    // Pressure has already been unpacked from x (and BCs corrected)
    // by unpackSolution() above
    gradp() = fvc::grad(p);

    // Interpolated effective viscosity for momentum stabilisation
    const surfaceScalarField nuEfff(fvc::interpolate(turbulence_->nuEff()));

    // Update momentum stabilisation
    momentumStabilisation().updateVector(U, &gradU());

    // The residual vector
    vectorField residual
    (
        fvc::laplacian(turbulence_->nuEff(), U)
      - gradp()
      - fvc::ddt(U)
      - fvc::div(phi, U)
      + momentumStabilisation().cellVector(&nuEfff, true)
    );

    if (Switch(fluidProperties().lookupOrDefault<Switch>("addDivPhiUDamping", false)))
    {
        residual -= 0.5*fvc::div(phiAbs)*U;
    }

    residual *= mesh.V();

    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0,
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Pressure residual
    scalarField pressureResidual(- fvc::div(U));

    // Approximate momentum operator used only to derive the rAUf
    // coefficient field. Use the same first-order upwind named scheme
    // as the assembled Jacobian's div term ("jacobian-div(phi,U)" in
    // fvSchemes) so that A.A() has a predictable sign structure: the
    // upwind div diagonal is the sum of outflow fluxes (>= 0), so
    // pressureStabUEqn.A() = laplacian.A() - V/dt - upwind_div.A()
    // stays consistently negative in the ddt/advection-dominated
    // regime regardless of how non-physical a Newton or MFFD trial U
    // is. This is what keeps rAUf = -1/A from flipping sign under
    // Newton/MFFD perturbations and tripping the RhieChow check.
    // The actual momentum residual (above) still uses the user-chosen
    // div scheme, so physical accuracy is unchanged
    fvVectorMatrix pressureStabUEqn
    (
        fvm::laplacian(turbulence_->nuEff(), U)
      - fvm::ddt(U)
      - fvm::div(phi, U, "jacobian-div(phi,U)")
    );

    rAUf() = -fvc::interpolate(1.0/pressureStabUEqn.A());

    pressureStabilisation().updateScalar(p, &gradp());
    pressureResidual += pressureStabilisation().cellScalar(&rAUf(), true);

    pressureResidual *= mesh.V();

    if (pRefCell_ != -1)
    {
        pressureResidual[pRefCell_] = pRefValue_ - p[pRefCell_];
    }

    // Row-scale the pressure residual.
    // pressureEqnScale_ == pressureScaleFactor_ today; the symbol is
    // used here so that any future twoMu-style augmentation lives in
    // one place
    if (mag(pressureEqnScale_ - 1.0) > SMALL)
    {
        pressureResidual *= pressureEqnScale_;
    }

    foamPetscSnesHelper::InsertFieldComponents<scalar>
    (
        pressureResidual, f, blockSize_ - 1
    );

    return 0;
}


label newtonIcoFluid::formJacobian
(
    Mat jac,
    const Vec x,
    const bool extrapolatedFlux
)
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    const fvMesh& mesh = this->mesh();

    // Copy x into the U and p fields (pHat -> p when scale != 1).
    // unpackSolution also corrects the boundary conditions on both
    // fields, so the rest of this routine can use them directly
    unpackSolution(x);

    volVectorField& U = const_cast<volVectorField&>(this->U());

    surfaceScalarField& phi = this->phi();
    const Switch forceImplicitFlux =
        fluidProperties().lookupOrDefault<Switch>("forceImplicitFlux", false);

    if (forceImplicitFlux || !extrapolatedFlux)
    {
        phi = fvc::interpolate(U) & mesh.Sf();
    }
    else if
    (
        Switch
        (
            fluidProperties().lookup("fluidFluxExtrapolationAlgorithm1")
        )
    )
    {
        phi =
            fvc::interpolate
            (
                2.0*U.oldTime() - U.oldTime().oldTime()
            ) & mesh.Sf();
    }
    else
    {
        phi =
            fvc::interpolate
            (
                2.25*U.oldTime()
              - 1.5*U.oldTime().oldTime()
              + 0.25*U.oldTime().oldTime().oldTime()
            ) & mesh.Sf();
    }

    const surfaceScalarField phiAbs("phiAbs", phi);

    if (mesh.changing())
    {
        fvc::makeRelative(phi, U);

        forAll(U.boundaryField(), patchI)
        {
            if (mesh.boundaryMesh()[patchI].type() == "wall")
            {
                boundaryFieldRef(phi)[patchI] = 0.0;
            }
        }
    }

    // Pressure has already been unpacked from x (and BCs corrected)
    // by unpackSolution() above
    volScalarField& p = const_cast<volScalarField&>(this->p());

    laminarTransport_.correct();
    turbulence_->correct();

    fvVectorMatrix UEqn
    (
        fvm::laplacian(turbulence_->nuEff(), U)
      - fvm::ddt(U)
    );

    UEqn.relax();

    if (extrapolatedFlux && !forceImplicitFlux)
    {
        UEqn -= fvm::div(phi, U, "jacobian-div(phi,U)");

        if (Switch(fluidProperties().lookupOrDefault<Switch>("addDivPhiUDamping", false)))
        {
            UEqn -= 0.5*fvm::Sp(fvc::div(phiAbs), U);
        }
    }
    else
    {
        foamPetscSnesHelper::InsertFvmDivPhiUIntoPETScMatrix
        (
            U,
            phi,
            jac,
            0,
            0,
            fluidModel::twoD() ? 2 : 3
        );
    }

    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
    (
        UEqn, jac, 0, 0, fluidModel::twoD() ? 2 : 3
    );

    // Add advection term to UEqn to allow the correct central coefficient to
    // be calculated for pressure stabilisation
    if (!extrapolatedFlux || forceImplicitFlux)
    {
        UEqn -= fvm::div(phi, U);
    }

    // Record the reciprocal of the central coefficient
    rAUf() = -fvc::interpolate(1.0/UEqn.A());

    fvScalarMatrix pEqn
    (
        pressureStabilisation().scalarJacobian(p, &rAUf(), true)
    );

    if (pRefCell_ != -1)
    {
#ifdef OPENFOAM_COM
        pEqn.setValues(labelList(1, pRefCell_), 0.0);
#else
        pEqn.setValues(labelList(1, pRefCell_), scalarField(1, 0.0));
#endif
        pEqn.diag()[pRefCell_] = -1.0;
    }

    // J_pp block: row is the pressure equation (row-scaled by
    // pressureEqnScale_) AND column is the pressure unknown pHat
    // (column-scaled by pressureUnknownScale_). Both factors apply
    const scalar ppScale = pressureEqnScale_*pressureUnknownScale_;
    if (mag(ppScale - 1.0) > SMALL)
    {
        scaleFvScalarMatrix(pEqn, ppScale);
    }

    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
    (
        pEqn, jac, blockSize_ - 1, blockSize_ - 1, 1
    );

    const word pressureStabType(pressureStabilisation().type());
    if
    (
        pressureStabType == "RhieChow"
     || pressureStabType == "diffStencilLaplacian"
    )
    {
        // Same J_pp block: row-scale by pressureEqnScale_ and
        // column-scale by pressureUnknownScale_
        foamPetscSnesHelper::InsertFvcDivGradInterpolateIntoPETScMatrix
        (
            p,
            rAUf(),
            jac,
            blockSize_ - 1,
            blockSize_ - 1,
           -ppScale*pressureStabilisation().scaleFactorJacobian()
        );
    }

    // J_pU block: pressure row, U column. Only row-scaling applies
    foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
    (
        p,
        U,
        jac,
        blockSize_ - 1,
        0,
        fluidModel::twoD() ? 2 : 3,
        pressureEqnScale_
    );

    // J_Up block: U row, pressure (pHat) column. Only column-scaling
    // applies. Helper defaults to scale = -1.0 (inserts -V*grad(p));
    // multiply by pressureUnknownScale_ for the pHat column
    foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
    (
        p,
        jac,
        0,
        blockSize_ - 1,
        fluidModel::twoD() ? 2 : 3,
        -pressureUnknownScale_
    );

    if (debug)
    {
        Info<< "End" << endl;
    }

    return 0;
}


label newtonIcoFluid::precondition
(
    Vec y,
    const Vec x
)
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    const fvMesh& mesh = this->mesh();
    const volVectorField& U = this->U();
    const volScalarField& p = this->p();
    const surfaceScalarField& phi = this->phi();

    // Build safe BC types for the correction fields.  Using
    // U.boundaryField().types() directly would propagate codedFixedValue
    // entries, which require a codeDict file that only exists inside the
    // original field specification.  For the correction fields, Dirichlet
    // boundaries are simply fixedValue (zero) and Neumann boundaries are
    // zeroGradient.
    wordList dUBCTypes(U.boundaryField().size());
    forAll(dUBCTypes, patchI)
    {
        if (U.boundaryField()[patchI].fixesValue())
        {
            dUBCTypes[patchI] = fixedValueFvPatchVectorField::typeName;
        }
        else
        {
            dUBCTypes[patchI] = U.boundaryField().types()[patchI];
        }
    }

    wordList dpBCTypes(p.boundaryField().size());
    forAll(dpBCTypes, patchI)
    {
        if (p.boundaryField()[patchI].fixesValue())
        {
            dpBCTypes[patchI] = fixedValueFvPatchScalarField::typeName;
        }
        else
        {
            dpBCTypes[patchI] = p.boundaryField().types()[patchI];
        }
    }

    volVectorField dU
    (
        IOobject
        (
            "dU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero),
        dUBCTypes
    );

    volScalarField dp
    (
        IOobject
        (
            "dp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", p.dimensions(), 0.0),
        dpBCTypes
    );

    vectorField momentumRhs(Foam::primitiveFieldRef(dU).size(), vector::zero);
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        momentumRhs,
        0,
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    scalarField pressureRhs(Foam::primitiveFieldRef(dp).size(), 0.0);
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        x, pressureRhs, blockSize_ - 1
    );

    // Use the same approximate momentum block as the assembled Jacobian's
    // segregated part.  The extra forceImplicitFlux tensor terms are left to
    // the assembled PETSc matrix path; this preconditioner stays deliberately
    // cheap and local.
    fvVectorMatrix UEqn
    (
        fvm::laplacian(turbulence_->nuEff(), dU)
      - fvm::ddt(dU)
      - fvm::div(phi, dU, "jacobian-div(phi,U)")
    );

    UEqn.relax();

    if (Switch(fluidProperties().lookupOrDefault<Switch>("addDivPhiUDamping", false)))
    {
        surfaceScalarField phiAbs("phiAbs", phi);

        if (mesh.changing())
        {
            fvc::makeAbsolute(phiAbs, U);
        }

        UEqn -= 0.5*fvm::Sp(fvc::div(phiAbs), dU);
    }

    UEqn.source() = momentumRhs;
    Foam::primitiveFieldRef(dU) = vector::zero;
    dU.correctBoundaryConditions();

#ifdef OPENFOAM_NOT_EXTEND
    const int oldVectorSolverDebug = SolverPerformance<vector>::debug;
    const int oldScalarSolverDebug = SolverPerformance<scalar>::debug;
    SolverPerformance<vector>::debug = 0;
    SolverPerformance<scalar>::debug = 0;
#else
    const int oldBlockLduDebug = blockLduMatrix::debug();
    blockLduMatrix::debug = 0;
#endif

    UEqn.solve(solverControls(mesh, "dU"));

    volScalarField rAU("rAU", 1.0/UEqn.A());

    // Match the pressure-block sign convention in formJacobian().
    rAUf() = -fvc::interpolate(rAU);

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(dU) & mesh.Sf()
    );

    fvScalarMatrix pEqn
    (
        pressureStabilisation().scalarJacobian(dp, &rAUf(), true)
    );

#ifdef OPENFOAM_NOT_EXTEND
    const scalarField divPhiHbyAV
    (
        fvc::div(phiHbyA)().primitiveField()*mesh.V()
    );
#else
    const scalarField divPhiHbyAV
    (
        fvc::div(phiHbyA)().internalField()*mesh.V()
    );
#endif

    // Both pressureRhs and divPhiHbyAV are mapped onto the same
    // row-scaled space the assembled pEqn lives in. pressureEqnScale_
    // is used here so the rename stays symmetric with formJacobian.
    // No pressureUnknownScale_ here: we deliberately solve for the
    // physical pressure correction dp and convert dp -> dpHat at the
    // end (see the dp /= pressureUnknownScale_ step below)
    scalarField pSource
    (
        pressureRhs
      + pressureEqnScale_*divPhiHbyAV
    );

    if (pRefCell_ != -1)
    {
#ifdef OPENFOAM_COM
        pEqn.setValues(labelList(1, pRefCell_), 0.0);
#else
        pEqn.setValues(labelList(1, pRefCell_), scalarField(1, 0.0));
#endif
        pEqn.diag()[pRefCell_] = -1.0;
        pSource[pRefCell_] = pressureRhs[pRefCell_];
    }

    if (mag(pressureEqnScale_ - 1.0) > SMALL)
    {
        scaleFvScalarMatrix(pEqn, pressureEqnScale_);
    }

    pEqn.source() = pSource;
    Foam::primitiveFieldRef(dp) = 0.0;
    dp.correctBoundaryConditions();

    pEqn.solve(solverControls(mesh, "dp"));

    dU += rAU*fvc::grad(dp);
    dU.correctBoundaryConditions();

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = oldVectorSolverDebug;
    SolverPerformance<scalar>::debug = oldScalarSolverDebug;
#else
    blockLduMatrix::debug = oldBlockLduDebug;
#endif

    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        Foam::primitiveFieldRef(dU),
        y,
        0,
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Convert the physical pressure correction dp produced by the
    // SIMPLE solve into the scaled correction dpHat = dp /
    // pressureUnknownScale_ before writing it into the output vector,
    // because PETSc consumes y in (U, pHat) coordinates
    if (mag(pressureUnknownScale_ - 1.0) > SMALL)
    {
        scalarField dpHat(Foam::primitiveFieldRef(dp));
        dpHat /= pressureUnknownScale_;
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            dpHat, y, blockSize_ - 1
        );
    }
    else
    {
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            Foam::primitiveFieldRef(dp), y, blockSize_ - 1
        );
    }

    if (debug)
    {
        Info<< "End" << endl;
    }

    return 0;
}

#endif // USE_PETSC


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
