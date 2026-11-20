/*---------------------------------------------------------------------------* \
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

#ifdef USE_PETSC

#include "foamPetscSnesSolver.H"
#include "IFstream.H"
#include "IOdictionary.H"
#include "petscUtils.H"
#include "petscErrorHandling.H"
#include <petsc/private/pcimpl.h>
#include "petscdmshell.h"
#include <cstring>

// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

PetscErrorCode formResidualFoamPetscSnesSolver
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Vec f,        // residual
    void *ctx     // user context
)
{
    Foam::foamPetscSnesSolverContext *user = (Foam::foamPetscSnesSolverContext *)ctx;

    PetscFunctionBeginUser;

    // Compute the residual
    if (user->system_.formResidual(f, x) != 0)
    {
        if (user->solver_.stopOnPetscError())
        {
            Foam::FatalError
                << "formResidual(ff, xx) returned an error code!"
                << Foam::abort(Foam::FatalError);
        }
        else
        {
            Foam::Warning
                << "formResidual(ff, xx) returned an error code!"
                << Foam::endl;
        }

        // Let SNES know about the error
        user->solver_.diverged() = true;
        PetscCall(SNESSetFunctionDomainError(snes));

        // Exit without an error code so SNES can exit "softly"
        PetscFunctionReturn(0);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode formJacobianFoamPetscSnesSolver
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Mat jac,      // Jacobian
    Mat B,        // Preconditioner matrix (can be jac)
    void *ctx     // user context
)
{
    // Note
    // The "-snes_lag_jacobian -2" PETSc option can be used to avoid
    // re-building the matrix

    // Access the OpenFOAM data
    Foam::foamPetscSnesSolverContext *user = (Foam::foamPetscSnesSolverContext *)ctx;

    // Zero the matrix but do not reallocate the space
    // For a nested matrix, this will zero all sub-matrices
    CHKERRQ(MatZeroEntries(B));

    // Populate the Jacobian => implemented by the system
    if (user->system_.formJacobian(B, x) != 0)
    {
        Foam::FatalError
            << "formJacobian(B, xx) returned an error code!"
            << Foam::abort(Foam::FatalError);
    }

    // Complete matrix assembly
    CHKERRQ(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));

    if (jac != B)
    {
        CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
    }

    return 0;
}


PetscErrorCode convergenceCheckFoamPetscSnesSolver
(
    SNES snes,
    PetscInt it,
    PetscReal xnorm,
    PetscReal gnorm,
    PetscReal fnorm,
    SNESConvergedReason *reason,
    void *ctx
)
{
  Foam::foamPetscSnesSolverContext *user = (Foam::foamPetscSnesSolverContext *)ctx;

  // PETSc default check
  PetscCall(SNESConvergedDefault(snes, it, xnorm, gnorm, fnorm, reason, NULL));
  if (*reason)
  {
      return 0;
  }

  if (user->solver_.diverged())
  {
      *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
      user->solver_.diverged() = false;
  }

  return 0;
}


struct physicsPCData
{
    KSP ksp{nullptr};
    bool useOperator{false};
};


static PetscErrorCode PCApplyPhysics(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    physicsPCData* data = static_cast<physicsPCData*>(pc->data);

    if (data && data->useOperator)
    {
        PetscCall(KSPSolve(data->ksp, x, y));

        KSPConvergedReason reason;
        PetscCall(KSPGetConvergedReason(data->ksp, &reason));

        if (reason < 0)
        {
            SETERRQ
            (
                PetscObjectComm((PetscObject)pc),
                PETSC_ERR_CONV_FAILED,
                "Inner physics operator KSP failed with reason %d",
                (int)reason
            );
        }

        PetscFunctionReturn(0);
    }

    // Retrieve the model context from the DM
    DM dm = nullptr;
    AssertPETSc(PCGetDM(pc, &dm));
    Foam::foamPetscSnesSolverContext* user = nullptr;

    if (dm)
    {
        AssertPETSc(DMGetApplicationContext(dm, (void**)&user));
    }

    if (!user)
    {
        Foam::FatalError
            << "The system context needs to be attached to the SNES "
            << "object via a DM"
            << Foam::abort(Foam::FatalError);
    }

    // Call the system precondition function
    AssertPETSc(user->system_.precondition(y, x));

    PetscFunctionReturn(0);
}


static PetscErrorCode PCSetUpPhysics(PC pc)
{
    PetscFunctionBeginUser;

    physicsPCData* data = static_cast<physicsPCData*>(pc->data);

    char pcType[32] = "model";
    PetscCall
    (
        PetscOptionsGetString
        (
            nullptr,
            nullptr,
            "-pc_physics_type",
            pcType,
            sizeof(pcType),
            nullptr
        )
    );

    const bool useModel = (std::strcmp(pcType, "model") == 0);
    data->useOperator = (std::strcmp(pcType, "operator") == 0);

    if (!useModel && !data->useOperator)
    {
        SETERRQ
        (
            PetscObjectComm((PetscObject)pc),
            PETSC_ERR_ARG_WRONG,
            "Unknown pc_physics_type '%s'; expected 'model' or 'operator'",
            pcType
        );
    }

    if (data->useOperator)
    {
        Mat mat = nullptr;
        Mat pmat = nullptr;
        PetscCall(PCGetOperators(pc, &mat, &pmat));

        if (!pmat)
        {
            SETERRQ
            (
                PetscObjectComm((PetscObject)pc),
                PETSC_ERR_ARG_NULL,
                "No preconditioning matrix available for pc_physics_type operator"
            );
        }

        if (!data->ksp)
        {
            PetscCall
            (
                KSPCreate(PetscObjectComm((PetscObject)pc), &data->ksp)
            );
            PetscCall(KSPSetOptionsPrefix(data->ksp, "pc_physics_operator_"));
        }

        PetscCall(KSPSetOperators(data->ksp, pmat, pmat));
        PetscCall(KSPSetType(data->ksp, KSPPREONLY));

        PC innerPc = nullptr;
        PetscCall(KSPGetPC(data->ksp, &innerPc));
        PetscCall(PCSetType(innerPc, PCLU));

        PetscCall(KSPSetFromOptions(data->ksp));
        PetscCall(KSPSetUp(data->ksp));
    }
    else if (data->ksp)
    {
        PetscCall(KSPReset(data->ksp));
    }

    PetscFunctionReturn(0);
}


static PetscErrorCode PCDestroyPhysics(PC pc)
{
    PetscFunctionBeginUser;

    physicsPCData* data = static_cast<physicsPCData*>(pc->data);

    if (data)
    {
        PetscCall(KSPDestroy(&data->ksp));
        delete data;
        pc->data = nullptr;
    }

    PetscFunctionReturn(0);
}


static PetscErrorCode PCCreatePhysics(PC pc)
{
    PetscFunctionBeginUser;

    pc->data = new physicsPCData;

    pc->ops->setup   = PCSetUpPhysics;
    pc->ops->apply   = PCApplyPhysics;
    pc->ops->destroy = PCDestroyPhysics;

    PetscFunctionReturn(0);
}


extern "C" PetscErrorCode PCRegisterPhysics()
{
    PetscFunctionBeginUser;

    // Register the name "physics" so users can select "-pc_type physics"
    PetscCall(PCRegister("physics", PCCreatePhysics));

    PetscFunctionReturn(0);
}


extern "C" PetscErrorCode PCPhysicsEnsureRegistered(void)
{
    static bool done = false;
    if (!done)
    {
        // Register the "physics" preconditioner
        PetscCall(PCRegisterPhysics());
        done = true;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(foamPetscSnesSolver, 0);


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

void foamPetscSnesSolver::initialisePetsc()
{
    PetscBool initialised = PETSC_FALSE;
    PetscInitialized(&initialised);

    if (!initialised)
    {
        // Initialise PETSc without any options file or fvSolution lookup
        PetscInitialize(NULL, NULL, NULL, NULL);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

autoPtr<IOdictionary> foamPetscSnesSolver::makeSolutionDictOverride() const
{
    if (fvSolutionLocation_.empty())
    {
        return autoPtr<IOdictionary>();
    }

    autoPtr<IOdictionary> solutionDictPtr
    (
        new IOdictionary
        (
            IOobject
            (
                "fvSolution",
                fvSolutionLocation_,
                mesh_.time(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    );

    if (!solutionDictPtr().headerOk())
    {
        solutionDictPtr.clear();
    }

    return solutionDictPtr;
}


bool foamPetscSnesSolver::populateOptionsFromDict
(
    const dictionary& solutionDict,
    const fileName& solutionDictPath
)
{
    if (!solutionDict.found("solvers"))
    {
        return false;
    }

    const dictionary& solversDict = solutionDict.subDict("solvers");

    if (!solversDict.found(fieldName_))
    {
        return false;
    }

    const dictionary& solverDict = solversDict.subDict(fieldName_);

    if (solverDict.empty())
    {
        return false;
    }

    if (solverDict.lookupOrDefault<word>("solver", "none") != "petsc")
    {
        return false;
    }

    Info<< "Reading the PETSc options from " << solutionDictPath
        << " for " << fieldName_ << endl;

    const dictionary& optionsDict = solverDict.subDict("options");

    // We will load (push) the empty options database and populate it
    // with options from the the dictionary, then unload (pop) the
    // database
    AssertPETSc(PetscOptionsPush(options_));

    // Populate the database
    PetscUtils::setFlags("", optionsDict, debug);

    // Unload (pop) the database
    AssertPETSc(PetscOptionsPop());

    return true;
}


void foamPetscSnesSolver::loadOptions()
{
    if (optionsLoaded_)
    {
        return;
    }

    if (!options_)
    {
        AssertPETSc(PetscOptionsCreate(&options_));
    }

    bool useDict = false;
    const fileName fvSolutionPath =
        fvSolutionLocation_.empty()
      ? fileName("fvSolution")
      : fvSolutionLocation_/"fvSolution";

    if (fvSolutionLocation_.empty())
    {
        useDict = populateOptionsFromDict(mesh_.solutionDict(), fvSolutionPath);
    }
    else
    {
        autoPtr<IOdictionary> solutionDictPtr = makeSolutionDictOverride();

        if (solutionDictPtr.valid())
        {
            useDict = populateOptionsFromDict(solutionDictPtr(), fvSolutionPath);
        }
    }

    if (!useDict)
    {
        fileName optionsFile(optionsFile_);

        Info<< "Reading the PETSc options from the " << optionsFile
            << " file" << endl;

        // Expand the options file name
        optionsFile.expand();

        // Check the options file exists
        IFstream is(optionsFile);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot find the PETSc options file: " << optionsFile
                << ". Either provide an option file or add 'solver petsc;' "
                << "to " << fvSolutionPath << "/solvers/" << fieldName_
                << " dictionary"
                << exit(FatalError);
        }

        // Populate the options database with the options file
        AssertPETSc
        (
            PetscOptionsInsertFile
            (
                PETSC_COMM_WORLD, options_, optionsFile.c_str(), PETSC_TRUE
            )
        );
    }

    optionsLoaded_ = true;
}


label foamPetscSnesSolver::initialiseSnes()
{
    if (snes_.s)
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    // Create the PETSc SNES object
    snes_.s = SNES();
    AssertPETSc(SNESCreate(PETSC_COMM_WORLD, &snes_.s));

    // Create user data context
    contextPtr_.set(new foamPetscSnesSolverContext(system_, *this));
    foamPetscSnesSolverContext& user = contextPtr_();

    // Set the user context
    AssertPETSc(SNESSetApplicationContext(snes_.s, &user));

    // Set the residual function
    AssertPETSc
    (
        SNESSetFunction(snes_.s, NULL, formResidualFoamPetscSnesSolver, &user)
    );

    // The system initialises A
    AssertPETSc(system_.initialiseJacobian(A_.m));

    // Set the Jacobian function
    AssertPETSc
    (
        SNESSetJacobian
        (
            snes_.s, A_.m, A_.m, formJacobianFoamPetscSnesSolver, &user
        )
    );

    // Set the convergence check function
    AssertPETSc
    (
        SNESSetConvergenceTest
        (
            snes_.s, convergenceCheckFoamPetscSnesSolver, &user, NULL
        )
    );

    // Register the solids4foam physics-based solver preconditioner, which can
    // be optionally called with "-pc_type physics"
    AssertPETSc(PCPhysicsEnsureRegistered());

    // Store a pointer to the solid model context in a DM, which is attached to
    // the SNES object
    {
        // Make a DMShell and store the user ctx on it
        DM dm;
        DMShellCreate(PETSC_COMM_WORLD, &dm);
        DMSetApplicationContext(dm, &user);

        // Attach the DM to SNES
        SNESSetDM(snes_.s, dm);
        DMDestroy(&dm); // SNES holds a ref
    }

    // Set solver options
    // Uses default options, can be overridden by command line options
    AssertPETSc(SNESSetFromOptions(snes_.s));

    // Attach the solid model context if the "physics" precondition is chosen
    // This needs to happen after SNESSetFromOptions
    // {
    //     Info<< "CHECKING for physics PC" << endl;
    //     KSP ksp; PC pc;
    //     AssertPETSc(SNESGetKSP(snes_.s, &ksp));
    //     AssertPETSc(KSPGetPC(ksp, &pc));

    //     PetscBool isPhysics = PETSC_FALSE;
    //     AssertPETSc
    //     (
    //         PetscObjectTypeCompare((PetscObject)pc, "physics", &isPhysics)
    //     );

    //     // if (isPhysics)
    //     {
    //         Info<< "ATTACHING user" << endl;
    //         // pc->data = &user;
    //         Info<< "ATTACHING user: DONE" << endl;
    //     }
    // }

    // The system initialises the solution vector
    AssertPETSc(system_.initialiseSolution(x_.v));

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamPetscSnesSolver::foamPetscSnesSolver
(
    petscSnesSystem& system,
    const word& fieldName,
    const fileName& optionsFile,
    const fvMesh& mesh,
    const Switch stopOnPetscError,
    const fileName& fvSolutionLocation
)
:
    system_(system),
    fieldName_(fieldName),
    optionsFile_(optionsFile),
    mesh_(mesh),
    fvSolutionLocation_(fvSolutionLocation),
    options_(nullptr),
    optionsLoaded_(false),
    stopOnPetscError_(stopOnPetscError),
    diverged_(false),
    forceSnesSolverStateRebuild_(false),
    snes_(),
    x_(),
    xBackup_(),
    A_(),
    contextPtr_(),
    snesHasRun_(false)
{
    initialisePetsc();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

foamPetscSnesSolver::~foamPetscSnesSolver()
{
    // Destroy the PETSc objects before the context they point to
    snes_.reset();
    x_.reset();
    xBackup_.reset();
    A_.reset();
    contextPtr_.clear();

    if (options_)
    {
        PetscOptionsDestroy(&options_);
    }

    // Note: PETSc is not finalised here, as other solvers may still be
    // using it
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void foamPetscSnesSolver::resetSnes()
{
    Info<< "Resetting SNES" << endl;

    // Note: the solution backup is kept, as it has the layout of the new
    // solution vector
    snes_.reset();
    x_.reset();
    A_.reset();
    contextPtr_.clear();

    // Let the system clear any state derived from the previous SNES object
    // before the new one is built
    system_.resetCustomSolverState();

    if (initialiseSnes() != 0)
    {
        FatalErrorInFunction
            << "initialiseSnes failed" << abort(FatalError);
    }

    // Reset SNES internal state
    // SNESReset(foamPetscSnesSolver::snes());
    // Reset KSP and PC
    // KSP ksp;
    // SNESGetKSP(foamPetscSnesSolver::snes(), &ksp);
    // // KSPReset(ksp);
    // // KSPSetUp(ksp);  // Optional but safe

    // PC pc;
    // KSPGetPC(ksp, &pc);
    // PCReset(pc);
    // PCSetUp(pc);    // Optional but safe

    // // Reset line search
    // SNESLineSearch ls;
    // SNESGetLineSearch(foamPetscSnesSolver::snes(), &ls);
    // SNESLineSearchReset(ls);
}


void foamPetscSnesSolver::resetSnesSolverState()
{
    if (!snes_)
    {
        return;
    }

    Info<< "Resetting PETSc SNES/KSP retry state" << endl;

    SNESLineSearch lineSearch = nullptr;
    AssertPETSc(SNESGetLineSearch(snes_.s, &lineSearch));

    if (lineSearch)
    {
        AssertPETSc(SNESLineSearchReset(lineSearch));
    }

    KSP ksp = nullptr;
    AssertPETSc(SNESGetKSP(snes_.s, &ksp));

    if (ksp)
    {
        AssertPETSc(KSPSetReusePreconditioner(ksp, PETSC_FALSE));
    }

    forceSnesSolverStateRebuild_ = true;
    diverged_ = false;
}


int foamPetscSnesSolver::solve(const bool returnOnSnesError)
{
    // Reset the stopOnPetscError flag
    //stopOnPetscError_ = !returnOnSnesError;

    // Initialise the SNES object
    if (!snes_)
    {
        if (initialiseSnes() != 0)
        {
            FatalErrorInFunction
                << "initialiseSnes failed" << abort(FatalError);
        }
    }

    loadOptions();

    // Load the correct options database
    AssertPETSc(PetscOptionsPush(options_));

    AssertPETSc(SNESSetFromOptions(snes_.s));

    // Allow the system to adjust PC, KSP, etc.
    system_.customiseSolver(snes_.s);

    PetscInt lagJacobian = 1;
    PetscInt lagPreconditioner = 1;
    const bool forceSolverStateRebuild = forceSnesSolverStateRebuild_;

    if (forceSolverStateRebuild)
    {
        AssertPETSc(SNESGetLagJacobian(snes_.s, &lagJacobian));
        AssertPETSc(SNESGetLagPreconditioner(snes_.s, &lagPreconditioner));

        if (lagJacobian == -1)
        {
            AssertPETSc(SNESSetLagJacobian(snes_.s, -2));
        }

        if (lagPreconditioner == -1)
        {
            AssertPETSc(SNESSetLagPreconditioner(snes_.s, -2));
        }

        if (lagJacobian < 0 || lagPreconditioner < 0)
        {
            Info<< "Forcing PETSc Jacobian/preconditioner rebuild" << endl;
        }

        forceSnesSolverStateRebuild_ = false;
    }

    // Set the snesHasRun flag
    snesHasRun_ = true;

    // Clear any divergence flag from a previous failed solve before entering
    // the next SNES convergence check.
    diverged_ = false;

    // Solve the nonlinear system
    AssertPETSc(SNESSolve(snes_.s, NULL, x_.v));

    if (forceSolverStateRebuild)
    {
        if (lagJacobian == -2)
        {
            lagJacobian = -1;
        }

        if (lagPreconditioner == -2)
        {
            lagPreconditioner = -1;
        }

        AssertPETSc(SNESSetLagJacobian(snes_.s, lagJacobian));
        AssertPETSc(SNESSetLagPreconditioner(snes_.s, lagPreconditioner));
    }

    // Un-load the options file
    AssertPETSc(PetscOptionsPop());

    // Check convergence
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes_.s, &reason);
    diverged_ = (reason < 0);

    if (reason == SNES_DIVERGED_FUNCTION_DOMAIN)
    {
        Info<< nl
            << "PETSc SNES solver returned a diverged function domain: "
            << "returning" << endl;

        if (returnOnSnesError)
        {
            return reason;
        }
    }
    else if (reason < 0)
    {
        Info<< nl
            << "PETSc SNES solver return error check disabled" << endl
            << "The SNES nonlinear solver did not converge." << nl
            << " PETSc SNES convergence error code: " << reason << nl
            << " PETSc SNES convergence reason: "
            << SNESConvergedReasons[reason] << endl;

        if (returnOnSnesError)
        {
            return reason;
        }
        else if (stopOnPetscError_)
        {
            FatalErrorInFunction
                << "Stopping because of the PETSc SNES error" << nl
                << "Set `stopOnPetscError` to `false` to continue on PETSc "
                << "SNES errors"
                << abort(FatalError);
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
