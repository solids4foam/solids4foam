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
        state.writeTime_ = false;

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
    volScalarField& p = this->p();
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

            foamPetscSnesHelper::InsertFieldComponents<vector>
            (
                U,
                foamPetscSnesHelper::solution(),
                blockSize_,
                fluidModel::twoD()
              ? makeList<label>({0,1})
              : makeList<label>({0,1,2})
            );
        }

        // Update the mesh
#ifdef OPENFOAM_COM
        mesh.controlledUpdate();
#else
        mesh.update();
#endif

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

        if (!adjustTimeStep)
        {
            FatalErrorInFunction
                << "PETSc SNES failed to converge and the previous time-step "
                << "state has been restored, but `adjustTimeStep` is disabled."
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

        static_cast<TimeState&>(time) = retryTimeState;
        time.setTime(oldTimeValue, oldTimeIndex);
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

        ++time;

        Info<< "Retrying the failed PETSc time step with deltaT = "
            << time.deltaTValue() << " at Time = "
            << time.timeName() << nl << endl;
    }

    // Access the raw solution data
    const PetscScalar *xx;
    VecGetArrayRead(foamPetscSnesHelper::solution(), &xx);

    // Retrieve the solution
    // Map the PETSc solution to the U field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        xx,
        U,
        0, // Location of U
        blockSize_,
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    U.correctBoundaryConditions();

    // Map the PETSc solution to the p field
    // p is located in the final component
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        xx,
        p,
        blockSize_ - 1, // Location of p component
        blockSize_
    );

    p.correctBoundaryConditions();

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi);

    // Restore the solution vector
    VecRestoreArrayRead(foamPetscSnesHelper::solution(), &xx);

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
    volVectorField& U = const_cast<volVectorField&>(this->U());
    volScalarField& p = const_cast<volScalarField&>(this->p());
    surfaceScalarField& phi = const_cast<surfaceScalarField&>(this->phi());

    // Copy x into the U field
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0,
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    U.correctBoundaryConditions();
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

    // Copy x into the p field
    scalarField& pI = p;
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        x, pI, blockSize_ - 1
    );

    p.correctBoundaryConditions();
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

    fvVectorMatrix pressureStabUEqn
    (
        fvm::laplacian(turbulence_->nuEff(), U)
      - fvm::ddt(U)
      - fvm::div(phi, U)
    );

    rAUf() = fvc::interpolate(1.0/pressureStabUEqn.A());

    pressureStabilisation().updateScalar(p, &gradp());
    pressureResidual -= pressureStabilisation().cellScalar(&rAUf(), true);

    pressureResidual *= mesh.V();

    if (pRefCell_ != -1)
    {
        pressureResidual[pRefCell_] = pRefValue_ - p[pRefCell_];
    }

    if (pressureScaleFactor_ != 1.0)
    {
        pressureResidual *= pressureScaleFactor_;
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

    volVectorField& U = const_cast<volVectorField&>(this->U());
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0,
        fluidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    U.correctBoundaryConditions();

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

    volScalarField& p = const_cast<volScalarField&>(this->p());
    scalarField& pI = p;
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        x, pI, blockSize_ - 1
    );

    p.correctBoundaryConditions();

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

    if (pressureScaleFactor_ != 1.0)
    {
        scaleFvScalarMatrix(pEqn, pressureScaleFactor_);
    }

    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
    (
        pEqn, jac, blockSize_ - 1, blockSize_ - 1, 1
    );

    foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
    (
        p,
        U,
        jac,
        blockSize_ - 1,
        0,
        fluidModel::twoD() ? 2 : 3,
        pressureScaleFactor_
    );

    foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
    (
        p,
        jac,
        0,
        blockSize_ - 1,
        fluidModel::twoD() ? 2 : 3
    );

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
