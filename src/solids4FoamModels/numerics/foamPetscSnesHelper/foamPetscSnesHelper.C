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

#include "foamPetscSnesHelper.H"
#include "processorFvPatch.H"
#include "symmetryFvPatchFields.H"
#include "leastSquaresVectors.H"
#include "fvm.H"
#include "IFstream.H"
#include "petscUtils.H"
#include "petscErrorHandling.H"
#include <petsc/private/pcimpl.h>
#include "petscdmshell.h"
#ifdef OPENFOAM_NOT_EXTEND
    #include "symmetryPlaneFvPatchFields.H"
#endif

// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

PetscErrorCode formResidualFoamPetscSnesHelper
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Vec f,        // residual
    void *ctx     // user context
)
{
    const PetscScalar *xx;
    PetscScalar       *ff;
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;

    PetscFunctionBeginUser;

    // Access x and f data
    CHKERRQ(VecGetArrayRead(x, &xx));
    CHKERRQ(VecGetArray(f, &ff));

    // Compute the residual
    if (user->solMod_.formResidual(f, x) != 0)
    {
        if (user->solMod_.stopOnPetscError())
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
        user->solMod_.diverged() = true;
        PetscCall(SNESSetFunctionDomainError(snes));

        // Exit without an error code so SNES can exit "softly"
        PetscFunctionReturn(0);
    }

    // Restore the solution and residual vectors
    CHKERRQ(VecRestoreArrayRead(x, &xx));
    CHKERRQ(VecRestoreArray(f, &ff));

    PetscFunctionReturn(0);
}


PetscErrorCode formJacobianFoamPetscSnesHelper
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

    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Access the OpenFOAM data
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;

    // Zero the matrix but do not reallocate the space
    // For a nested matrix, this will zero all sub-matrices
    CHKERRQ(MatZeroEntries(B));

    // Populate the Jacobian => implemented by the solid model
    if (user->solMod_.formJacobian(B, x) != 0)
    {
        Foam::FatalError
            << "formJacobian(B, xx) returned an error code!"
            << Foam::abort(Foam::FatalError);
    }

    // Restore solution vector
    CHKERRQ(VecRestoreArrayRead(x, &xx));

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


PetscErrorCode convergenceCheckFoamPetscSnesHelper
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
  appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;

  // PETSc default check
  PetscCall(SNESConvergedDefault(snes, it, xnorm, gnorm, fnorm, reason, NULL));
  if (*reason)
  {
      return 0;
  }

  if (user->solMod_.diverged())
  {
      *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
      user->solMod_.diverged() = false;
  }

  return 0;
}


static PetscErrorCode PCApplyPhysics(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    // Retreive the solid model context from the DM
    DM dm = nullptr;
    AssertPETSc(PCGetDM(pc, &dm));
    appCtxfoamPetscSnesHelper* user = nullptr;

    if (dm)
    {
        AssertPETSc(DMGetApplicationContext(dm, (void**)&user));
    }

    if (!user)
    {
        Foam::FatalError
            << "The solid model context needs to be attached to the SNES "
            << "object via a DM"
            << Foam::abort(Foam::FatalError);
    }

    // Call the solid model precondition function
    AssertPETSc(user->solMod_.precondition(y, x));

    PetscFunctionReturn(0);
}


static PetscErrorCode PCSetUpPhysics(PC pc)
{
    PetscFunctionBeginUser;

    // // Retreive the solid model context from the DM
    // DM dm = nullptr;
    // AssertPETSc(PCGetDM(pc, &dm));
    // appCtxfoamPetscSnesHelper* user = nullptr;

    // if (dm)
    // {
    //     AssertPETSc(DMGetApplicationContext(dm, (void**)&user));
    // }

    // if (!user)
    // {
    //     Foam::FatalError
    //         << "The solid model context needs to be attached to the SNES "
    //         << "object via a DM"
    //         << Foam::abort(Foam::FatalError);
    // }

    // Read any runtime knobs and set them in my solid model
    // char variant[64] = "schur";
    // PetscOptionsGetString
    // (
    //     nullptr, nullptr, "-pc_physics_variant", variant, sizeof(variant), nullptr
    // );
    // PetscInt steps = 3;
    // PetscReal tol = 1e-1;
    // PetscOptionsGetInt(nullptr, nullptr, "-pc_physics_steps", &steps, nullptr);
    // PetscOptionsGetReal(nullptr, nullptr, "-pc_physics_tol", &tol, nullptr);

    PetscFunctionReturn(0);
}


static PetscErrorCode PCDestroyPhysics(PC pc)
{
    PetscFunctionBeginUser;

    // Nothing to do

    PetscFunctionReturn(0);
}


static PetscErrorCode PCCreatePhysics(PC pc)
{
    PetscFunctionBeginUser;

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

defineTypeNameAndDebug(foamPetscSnesHelper, 0);

// * * * * * * * * * * * * * * * Private Function  * * * * * * * * * * * * * //

void foamPetscSnesHelper::makeNeiProcFields(const fvMesh& mesh) const
{
    if
    (
        neiProcGlobalIDs_.size() != 0
     || neiProcVolumes_.size() != 0
    )
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    neiProcGlobalIDs_.setSize(mesh.boundaryMesh().size());
    neiProcVolumes_.setSize(mesh.boundaryMesh().size());

    PtrList<labelList>& neiProcGlobalIDs = neiProcGlobalIDs_;
    PtrList<scalarField>& neiProcVolumes = neiProcVolumes_;

    const scalarField& VI = mesh.V();

    // Send the data
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];

        if (fp.type() == "processor")
        {
            // Take a copy of the faceCells (local IDs) and convert them to
            // global IDs
            labelList globalFaceCells(fp.faceCells());
#ifdef OPENFOAM_COM
            foamPetscSnesHelper::globalCells().inplaceToGlobal(globalFaceCells);
#else
            forAll(globalFaceCells, cI)
            {
                const label localCellID = globalFaceCells[cI];
                globalFaceCells[cI] =
                    foamPetscSnesHelper::globalCells().toGlobal(localCellID);
            }
#endif

            // Send global IDs to the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(fp);
            procPatch.send
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );

            // Construct a list of the volumes of cells at the patch
            scalarField patchVols(fp.size());
            const labelUList& faceCells = fp.faceCells();
            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                patchVols[faceI] = VI[cellID];
            }

            // Send patch cell volumes to the neighbour proc
            procPatch.send
            (
                Pstream::commsTypes::blocking, patchVols
            );
        }
    }

    // Receive the data
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];

        if (fp.type() == "processor")
        {
            neiProcGlobalIDs.set(patchI, new labelList(fp.size()));
            labelList& globalFaceCells = neiProcGlobalIDs[patchI];

            // Receive global IDs from the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(fp);
            procPatch.receive
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );

            // Receive patch cell volmes from the neighbour proc
            neiProcVolumes.set(patchI, new scalarField(fp.size()));
            scalarField& patchNeiVols = neiProcVolumes_[patchI];
            procPatch.receive
            (
                Pstream::commsTypes::blocking, patchNeiVols
            );
        }
    }
}


const PtrList<labelList>& foamPetscSnesHelper::neiProcGlobalIDs
(
    const fvMesh& mesh
) const
{
    if (neiProcGlobalIDs_.size() == 0)
    {
        makeNeiProcFields(mesh);
    }

    return neiProcGlobalIDs_;
}


const PtrList<scalarField>& foamPetscSnesHelper::neiProcVolumes
(
    const fvMesh& mesh
) const
{
    if (neiProcVolumes_.size() == 0)
    {
        makeNeiProcFields(mesh);
    }

    return neiProcVolumes_;
}


const leastSquaresS4fVectors& foamPetscSnesHelper::lsVectors
(
    const volScalarField& p
) const
{
    // Lookup and return vectors if they exist for this field
    // Otherwise create and return them
    const word lsName("leastSquaresVectors" + p.name());
    if (p.mesh().foundObject<leastSquaresS4fVectors>(lsName))
    {
        return p.mesh().lookupObject<leastSquaresS4fVectors>(lsName);
    }

    boolList useBoundaryFaceValues(p.mesh().boundary().size(), false);
    forAll(p.mesh().boundary(), patchI)
    {
        if (p.boundaryField()[patchI].fixesValue())
        {
            useBoundaryFaceValues[patchI] = true;
        }
    }

#ifdef OPENFOAM_COM
    return leastSquaresS4fVectors::New(lsName, p.mesh(), useBoundaryFaceValues);
#else
    return leastSquaresS4fVectors::New(p.mesh(), useBoundaryFaceValues);
#endif
}


label foamPetscSnesHelper::initialiseSnes()
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
    snesUserPtr_.set(new appCtxfoamPetscSnesHelper(*this));
    appCtxfoamPetscSnesHelper& user = snesUserPtr_();

    // Set the user context
    AssertPETSc(SNESSetApplicationContext(snes_.s, &user));

    // Set the residual function
    AssertPETSc
    (
        SNESSetFunction(snes_.s, NULL, formResidualFoamPetscSnesHelper, &user)
    );

    // The derived class initialises A
    AssertPETSc(initialiseJacobian(A_.m));

    // Set the Jacobian function
    AssertPETSc
    (
        SNESSetJacobian
        (
            snes_.s, A_.m, A_.m, formJacobianFoamPetscSnesHelper, &user
        )
    );

    // Set the convergence check function
    AssertPETSc
    (
        SNESSetConvergenceTest
        (
            snes_.s, convergenceCheckFoamPetscSnesHelper, &user, NULL
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

    // The derived class initialises the solution vector
    AssertPETSc(initialiseSolution(x_.v));

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


foamPetscSnesHelper::foamPetscSnesHelper
(
    const word& fieldName,
    fileName optionsFile,
    const fvMesh& mesh,
    const solutionLocation& location,
    const Switch stopOnPetscError,
    const Switch initialise
)
:
    location_(location),
    initialised_(initialise),
    options_(nullptr),
    stopOnPetscError_(stopOnPetscError),
    diverged_(false),
    snes_(),
    x_(),
    xBackup_(),
    A_(),
    snesUserPtr_(),
    globalCellsPtr_
    (
        location == solutionLocation::CELLS
      ? new globalIndex(mesh.nCells())
      : nullptr
    ),
    globalPointsPtr_
    (
        location == solutionLocation::POINTS
      ? new globalPointIndices(mesh)
      : nullptr
    ),
    neiProcGlobalIDs_(),
    neiProcVolumes_(),
    snesHasRun_(false)
{
    if (initialise)
    {
        // Initialise PETSc without any options file
        PetscInitialize(NULL, NULL, NULL, NULL);

        // Create and store the options database from a dedicated options file
        // or from an OpenFOAM dict

        // Lookup the solver in fvSolution and check if PETSc is specified
        bool useDict = false;
#ifdef OPENFOAM_COM
        if (mesh.solversDict().found(fieldName))
        {
            const dictionary& solverDict = mesh.solverDict(fieldName);
#else
        if (mesh.solutionDict().subDict("solvers").found(fieldName))
        {
            const dictionary& solverDict =
                mesh.solutionDict().subDict("solvers").subDict(fieldName);
#endif

            if (!solverDict.empty())
            {
                if
                (
                    solverDict.lookupOrDefault<word>("solver", "none")
                 == "petsc"
                )
                {
                    useDict = true;
                }
            }
        }

        // Create an empty options database
        PetscOptionsCreate(&options_);

        // Populate the options database from the dict or options file
        if (useDict)
        {
            Info<< "Reading the PETSc options from fvSolution" << endl;

            // We will load (push) the empty options database and populate it
            // with options from the the dictionary, then unload (pop) the
            // database
#ifdef OPENFOAM_NOT_EXTEND
            const dictionary& optionsDict =
                mesh.solverDict(fieldName).subDict("options");
#else
            const dictionary& optionsDict =
                mesh.solutionDict().solverDict(fieldName).subDict("options");
#endif

            // Load (push) the database
            AssertPETSc(PetscOptionsPush(options_));

            // Populate the database
            PetscUtils::setFlags("", optionsDict, debug);

            // Unload (pop) the database
            AssertPETSc(PetscOptionsPop());
        }
        else
        {
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
                    << ". Either provide an option file or add the PETSc "
                    << "solver settings to fvSolution/solvers/" << fieldName
                    << exit(FatalError);
            }

            // Populate the options database with the options file
            PetscOptionsInsertFile
            (
                PETSC_COMM_WORLD, options_, optionsFile.c_str(), PETSC_TRUE
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

foamPetscSnesHelper::~foamPetscSnesHelper()
{
    if (initialised_)
    {
        PetscOptionsDestroy(&options_);
        snesUserPtr_.clear();

        WarningInFunction
            << "Find a neat solution to finalise PETSc only once"
            << endl;
        //PetscFinalize();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void foamPetscSnesHelper::resetSnes()
{
    Info<< "Resetting SNES" << endl;

    snes_.reset();
    x_.reset();
    A_.reset();
    snesUserPtr_.clear();

    if (initialiseSnes() != 0)
    {
        FatalErrorInFunction
            << "initialiseSnes failed" << abort(FatalError);
    }

    // Reset SNES internal state
    // SNESReset(foamPetscSnesHelper::snes());
    // Reset KSP and PC
    // KSP ksp;
    // SNESGetKSP(foamPetscSnesHelper::snes(), &ksp);
    // // KSPReset(ksp);
    // // KSPSetUp(ksp);  // Optional but safe

    // PC pc;
    // KSPGetPC(ksp, &pc);
    // PCReset(pc);
    // PCSetUp(pc);    // Optional but safe

    // // Reset line search
    // SNESLineSearch ls;
    // SNESGetLineSearch(foamPetscSnesHelper::snes(), &ls);
    // SNESLineSearchReset(ls);
}


label foamPetscSnesHelper::initialiseJacobian
(
    Mat& jac,
    const fvMesh& mesh,
    const label blockSize,
    const bool createMat
)
{
    // Number of local block equations
    label blockn = -1;

    // Number of local scalar equations, i.e. blockn*blockSize
    label n = -1;

    // Number of scalar equations across all processors
    label N = -1;

    if (location_ == solutionLocation::CELLS)
    {
        // Set the number local block equations
        blockn = mesh.nCells();

        // Set the number local scalar equations
        n = mesh.nCells()*blockSize;

        // Set the global system size
        N = returnReduce(n, sumOp<label>());
    }
    else if (location_ == solutionLocation::POINTS)
    {
        // Note: the size of x is, in general, not equal to the number of
        // points on this processor, as x only contains the points owned by
        // this processor. To access the values not-owned by the proc, we
        // need to request the values from the other processors

        // Take references for brevity and efficiency
        const labelList& localToGlobalPointMap =
            globalPoints().localToGlobalPointMap();
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();

        // Find size of global system
        // i.e. the highest global point index + 1
        const label blockN = gMax(localToGlobalPointMap) + 1;
        N = blockSize*blockN;

        // Find the start and end global point indices for this proc
        label blockStartID = N;
        label blockEndID = -1;
        forAll(ownedByThisProc, pI)
        {
            if (ownedByThisProc[pI])
            {
                blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
                blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
            }
        }

        // Find size of local system, i.e. the range of points owned by this
        // proc
        blockn = blockEndID - blockStartID + 1;
        n = blockSize*blockn;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = " << location_
            << exit(FatalError);
    }

    // Set the Jacobian matrix size
    if (createMat)
    {
        MatCreate(PETSC_COMM_WORLD, &jac);
        MatSetFromOptions(jac);
        MatSetSizes(jac, n, n, N, N);
        MatSetType(jac, MATMPIAIJ);
    }

    // Set the block size
    AssertPETSc(MatSetBlockSize(jac, blockSize));

    // Count the number of non-zeros in the matrix
    // Note: we assume a compact stencil, i.e. face only face neighbours

    // Number of on-processor non-zeros per row
    // Initialise d_nnz to one
    labelList d_nnz(blockn, 1);

    // Number of off-processor non-zeros per row
    // Initialise o_nnz to zero
    labelList o_nnz(blockn, 0);

    // Count the neighbours
    if (location_ == solutionLocation::CELLS)
    {
        // Count neighbours sharing an internal face
        const Foam::labelUList& own = mesh.owner();
        const Foam::labelUList& nei = mesh.neighbour();
        forAll(own, faceI)
        {
            const Foam::label ownCellID = own[faceI];
            const Foam::label neiCellID = nei[faceI];
            d_nnz[ownCellID]++;
            d_nnz[neiCellID]++;
        }

        // Count off-processor neighbour cells
        forAll(mesh.boundary(), patchI)
        {
            if (mesh.boundary()[patchI].type() == "processor")
            {
                const Foam::labelUList& faceCells =
                    mesh.boundary()[patchI].faceCells();

                forAll(faceCells, fcI)
                {
                    const Foam::label cellID = faceCells[fcI];
                    o_nnz[cellID]++;
                }
            }
            else if (mesh.boundary()[patchI].coupled())
            {
                // Other coupled boundaries are not implemented
                Foam::FatalError
                    << "Coupled boundary are not implemented, except for"
                    << " processor boundaries" << Foam::abort(Foam::FatalError);
            }
        }
    }
    else if (location_ == solutionLocation::POINTS)
    {
        // Take references
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();
        const labelList& stencilSizeOwned = globalPoints().stencilSizeOwned();
        const labelList& stencilSizeNotOwned =
            globalPoints().stencilSizeNotOwned();

        // Reset d_nnz to 0
        d_nnz = 0;

        // We can optionally track the max on-core non-zeros to zero
        //label d_nz = 0;

        // Count non-zeros
        label rowI = 0;
        forAll(stencilSizeOwned, blockRowI)
        {
            if (ownedByThisProc[blockRowI])
            {
                const label nCompOwned = stencilSizeOwned[blockRowI];
                const label nCompNotOwned = stencilSizeNotOwned[blockRowI];

                d_nnz[rowI] += nCompOwned;
                o_nnz[rowI++] += nCompNotOwned;

                //d_nz = max(d_nz, nCompOwned);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = " << location_
            << exit(FatalError);
    }

    // Allocate parallel matrix
    //AssertPETSc(MatMPIAIJSetPreallocation(jac, 0, d_nnz, 0, o_nnz));
    // Allocate parallel matrix with the same conservative stencil per node
    //AssertPETSc(MatMPIAIJSetPreallocation(jac, d_nz, NULL, 0, NULL));
    AssertPETSc(MatMPIBAIJSetPreallocation
    (
        jac, blockSize, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    MatSetOption(jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    return 0;
}


label foamPetscSnesHelper::initialiseSolution
(
    Vec& x,
    const fvMesh& mesh, // store this?
    const label blockSize,
    const bool createVec
)
{
    if (createVec)
    {
        // Number of local scalar equations
        label n = -1;

        // Number of scalar equations across all processors
        label N = -1;

        if (location_ == solutionLocation::CELLS)
        {
            // Set the number local scalar equations
            n = mesh.nCells()*blockSize;

            // Set the global system size
            N = returnReduce(n, sumOp<label>());
        }
        else if (location_ == solutionLocation::POINTS)
        {
            // Note: the size of x is, in general, not equal to the number of
            // points on this processor, as x only contains the points owned by
            // this processor. To access the values not-owned by the proc, we
            // need to request the values from the other processors

            // Take references for brevity and efficiency
            const labelList& localToGlobalPointMap =
                globalPoints().localToGlobalPointMap();
            const boolList& ownedByThisProc = globalPoints().ownedByThisProc();

            // Find size of global system
            // i.e. the highest global point index + 1
            const label blockN = gMax(localToGlobalPointMap) + 1;
            N = blockSize*blockN;

            // Find the start and end global point indices for this proc
            label blockStartID = N;
            label blockEndID = -1;
            forAll(ownedByThisProc, pI)
            {
                if (ownedByThisProc[pI])
                {
                    blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
                    blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
                }
            }

            //const label startID = blockSize*blockStartID;
            //const label endID = blockSize*(blockEndID + 1) - 1;

            // Find size of local system, i.e. the range of points owned by this
            // proc
            const label blockn = blockEndID - blockStartID + 1;
            n = blockSize*blockn;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown solution location = " << location_
                << exit(FatalError);
        }

        x = Vec();
        AssertPETSc(VecCreate(PETSC_COMM_WORLD, &x));
        AssertPETSc(VecSetSizes(x, n, N));
        AssertPETSc(VecSetType(x, VECMPI));
    }

    AssertPETSc(VecSetBlockSize(x, blockSize));
    AssertPETSc(PetscObjectSetName((PetscObject) x, "Solution"));
    AssertPETSc(VecSetFromOptions(x));
    AssertPETSc(VecZeroEntries(x));

    return 0;
}


label foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
(
    const volScalarField& p,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const bool flipSign
) const
{
    // Get reference to least square vectors
    const fvMesh& mesh = p.mesh();
    const leastSquaresS4fVectors& lsv = lsVectors(p);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const scalarField& VI = mesh.V();

    const scalar sign = flipSign ? -1.0 : 1.0;

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise the block coefficient
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    forAll(own, faceI)
    {
        // Local block row ID
        const label ownCellID = own[faceI];

        // Local block column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            foamPetscSnesHelper::globalCells().toGlobal(neiCellID);

        // Explicit gradient
        // lsGrad[ownCellID] += ownLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);

        // Insert coefficients for block row globalBlockRowI

        // lsGrad[ownCellID] += ownLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // mat(ownCellID, ownCellID) -= VI[ownCellID]*ownLs[faceI];
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            values[cmptI*blockSize + colOffset] =
                -sign*VI[ownCellID]*ownLs[faceI][cmptI];
        }
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                values.cdata(),
                ADD_VALUES
            )
        );

        // lsGrad[ownCellID] += ownLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // mat(ownCellID, neIFaceI) += VI[ownCellID]*ownLs[faceI];
        // Flip the sign of the values
        for (label i = 0; i < nCoeffCmpts; ++i)
        {
            values[i] = -values[i];
        }
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockRowI, 1, &globalBlockColI, values.cdata(),
                ADD_VALUES
            )
        );

        // Insert coefficients for block row globalBlockColI

        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        //mat(neIFaceI, neIFaceI) += VI[neiCellID]*neiLs[faceI];
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            values[cmptI*blockSize + colOffset] =
                -sign*VI[neiCellID]*neiLs[faceI][cmptI];
        }
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockColI, 1, &globalBlockColI, values.cdata(),
                ADD_VALUES
            )
        );

        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // mat(neIFaceI, ownCellID) -= VI[neiCellID]*neiLs[faceI];
        // Flip the sign of the values
        for (label i = 0; i < nCoeffCmpts; ++i)
        {
            values[i] = -values[i];
        }
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockColI, 1, &globalBlockRowI, values.cdata(),
                ADD_VALUES
            )
        );
    }

    // Boundary face contributions
    //const boolList& useBoundaryFaceValues = lsv.useBoundaryFaceValues();
    const boolList useBoundaryFaceValues(mesh.boundary().size(), true);
    forAll(mesh.boundary(), patchI)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchI];
        const labelUList& faceCells = mesh.boundary()[patchI].faceCells();
        const fvPatch& fp = mesh.boundary()[patchI];
        if (fp.type() == "processor")
        {
            const labelList& neiGlobalFaceCells =
                neiProcGlobalIDs(mesh)[patchI];
            const scalarField& patchNeiVols = neiProcVolumes(mesh)[patchI];
            const vectorField& patchNeiLs = neiLs.boundaryField()[patchI];
            // const vectorField patchNeiLs(... patchNeighbourField());

            forAll(fp, patchFaceI)
            {
                // Local block row ID
                const label ownCellID = faceCells[patchFaceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                // On-proc diagonal coefficient
                // mat(ownCellID, ownCellID) -= VI[ownCellID]*ownLs[faceI];
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize + colOffset] =
                        -sign*VI[ownCellID]*patchOwnLs[patchFaceI][cmptI];
                }
                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );

                // Neighbour global cell ID
                const label globalBlockColI = neiGlobalFaceCells[patchFaceI];

                // Off-proc off-diagonal coefficient
                // mat(ownCellID, neiCellID) += VI[ownCellID]*neiLs[faceI];
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize + colOffset] =
                        sign*patchNeiVols[patchFaceI]
                       *patchNeiLs[patchFaceI][cmptI];
                }
                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );
            }
        }
        else if (fp.coupled()) // coupled but not processor
        {
            FatalErrorInFunction
                << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        else if
        (
            isA<symmetryPolyPatch>(fp.patch())
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(fp.patch())
#endif
        )
        {
            // The delta in scalar p across the symmetry is zero by definition
            // so the symmetry plane does not contribution coefficients
        }
        else
        {
            if (useBoundaryFaceValues[patchI])
            {
                forAll(faceCells, patchFaceI)
                {
                    // Explicit calculation
                    // lsGrad[faceCells[patchFaceI]] +=
                    //      patchOwnLs[patchFaceI]
                    //     *(patchVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);

                    // Subtract patchOwnLs[patchFaceI] from (faceCellI, faceCellI)
                    // Nothing else to do as patch value is a known value

                    // Local block row ID
                    const label ownCellID = faceCells[patchFaceI];

                    // Global block row ID
                    const label globalBlockRowI =
                        foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                    // mat(ownCellID, ownCellID) -=
                    //     VI[ownCellID]*patchOwnLs[patchFaceI];
                    for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                    {
                        values[cmptI*blockSize + colOffset] =
                            -sign*VI[ownCellID]*patchOwnLs[patchFaceI][cmptI];
                    }
                    AssertPETSc
                    (
                        MatSetValuesBlocked
                        (
                            jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                            values.cdata(),
                            ADD_VALUES
                        )
                    );
                }
            }
        }
    }

    return 0;
}


label foamPetscSnesHelper::InsertFvmDivPhiUIntoPETScMatrix
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const bool flipSign
) const
{
    // Take references for efficiency and brevity
    const fvMesh& mesh = U.mesh();
    const scalarField& phiI = phi;
    const vectorField& UI = U;
    const vectorField& SfI = mesh.Sf();
    const scalarField& wI = mesh.weights();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalar sign = flipSign ? -1 : 1;
    tensor coeff = tensor::zero;

    // Add segregated component of convection term
    {
        fvVectorMatrix UEqn
        (
          - fvm::div(phi, U, "jacobian-div(phi,U)")
        );

        UEqn.relax();

        // Convert fvMatrix matrix to PETSc matrix
        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
        (
            UEqn, jac, 0, 0, nScalarEqns
        );
    }

    // The first-order upwind discretisation for cell P is
    //     div(phi,U)_P = \sum phi*Uf
    // where
    //     Uf = Up     if phi > 0
    //     Uf = Un     otherwise
    // and
    //     phi = Sf & [w*Up + (1 - w)*Un]
    //
    // The Newton linearisation of div(phi,U) for cell P is
    //     \sum (d/dUp)(phi*Uf) = \sum phi*(dUf/dUp) + Uf*(dphi/dUp)
    //
    // Using the following
    //     dUp/dUp = I
    //     dUn/dUn = I
    //     dphi/dUp = w*Sf
    //     dphi/dUn = (1 - w)*Sf
    // we can calculate the block coefficient for a face:
    // if phi > 0
    //     Uf = Up
    //     (d/dUp)(phi*Up) = phi*I + w*Up*Sf        => diagonal
    //     (d/dUn)(phi*Up) = (1 - w)*Up*Sf          => off-diagonal
    // else
    //     Uf = Un
    //     (d/dUp)(phi*Un) = w*Un*Sf                => diagonal
    //     (d/dUn)(phi*Un) = phi*I + (1 - w)*Un*Sf  => off-diagonal
    //
    // The w*Sf*Up term is missing from fvm::div(phi(), U) as it requires a
    // coupled solver. So we will add the w*Sf*Up term here
    // From the neighbours perspective, the coefficient is (1 - w)*Sf*Un

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Prepare coeff array
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    forAll(phiI, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            foamPetscSnesHelper::globalCells().toGlobal(neiCellID);

        // Take references for brevity, clarity and efficiency
        const scalar w = wI[faceI];
        const vector& Sf = SfI[faceI];
        const vector& Up = UI[ownCellID];
        const vector& Un = UI[neiCellID];

        if (phiI[faceI] > 0)
        {
            // Add w*Up*Sf to owner eqn
            coeff = sign*w*Up*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add (1 - w)*Up*Sf as nei to contribution to own eqn
            coeff = sign*(1.0 - w)*Up*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Up*Sf to the neighbour diagonal
            // Flip the sign as Sf is opposite
            coeff = -sign*(1.0 - w)*Up*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -w*Up*Sf as the own contribution to the nei eqn
            // Flip the sign as Sf is opposite
            coeff = -sign*w*Up*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );
        }
        else
        {
            // Add w*Un*Sf to owner diagonal
            coeff = sign*w*Un*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add (1 - w)*Un*Sf as nei to contribution to own eqn
            // coeff = sign*w*Un*Sf;
            coeff = sign*(1.0 - w)*Un*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Un*Sf to neighbour diagonal
            // Flip the sign as Sf is opposite
            coeff = -sign*(1.0 - w)*Un*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -w*Un*Sf as own to contribution to nei eqn
            // Flip the sign as Sf is opposite
            coeff = -sign*w*Un*Sf;
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );
        }
    }

    // Boundary contribution for outlets
    forAll(mesh.boundary(), patchI)
    {
        const fvPatchVectorField& pU = U.boundaryField()[patchI];
        const fvPatch& patch = pU.patch();
        const vectorField& pSf = patch.Sf();
        //const scalarField& pphi = phi().boundaryField()[patchI];
        const scalarField& pw = patch.weights();
        const labelUList& fc = patch.faceCells();

        if (patch.coupled())
        {
            notImplemented("div(phi,U) additional term for processors");
        }
        else if
        (
            patch.type() != "empty"
         && !isA<symmetryFvPatchField<vector>>(pU)
#ifdef OPENFOAM_NOT_EXTEND
         && !isA<symmetryPlaneFvPatchField<vector>>(pU)
#endif
         && !pU.fixesValue()
        )
        {
            forAll(fc, faceI)
            {
                // Local row ID
                const label ownCellID = fc[faceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                // Add w*Sf*Up to owner eqn
                coeff = sign*pw[faceI]*UI[ownCellID]*pSf[faceI];
                for (label i = 0; i < (blockSize - 1); ++i)
                {
                    for (label j = 0; j < (blockSize - 1); ++j)
                    {
                        // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                        // the 4x4 (or 3x3 in 2-D) values matrix
                        values[(i + rowOffset)*blockSize + j + colOffset] =
                            coeff[i*3 + j];
                    }
                }
                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );
            }
        }
    }

    return 0;
}


label foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
(
    const volScalarField& p,
    const volVectorField& U,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const bool flipSign
) const
{
    // Take references for efficiency and brevity
    const fvMesh& mesh = p.mesh();
    const scalar sign = flipSign ? -1 : 1;

    for (label cmptI = 0; cmptI < 3; ++cmptI)
    {
        if (nScalarEqns == 2 && cmptI == 2)
        {
            break;
        }

        fvScalarMatrix divUCoeffs(p, dimArea*dimPressure);

        const vectorField& Sf = mesh.Sf();
        const surfaceScalarField& weights = mesh.weights();
        const scalarField& w = weights;

        scalarField& upper = divUCoeffs.upper();
        scalarField& lower = divUCoeffs.lower();

        lower = w*Sf.component(cmptI);
        upper = lower - Sf.component(cmptI);

        divUCoeffs.negSumDiag();

        // Boundary contributions
        scalarField& diag = divUCoeffs.diag();
        forAll(mesh.boundary(), patchI)
        {
            const fvPatchVectorField& pU = U.boundaryField()[patchI];
            const fvPatch& patch = pU.patch();
            const vectorField& Sf = patch.Sf();
            const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
            const labelUList& fc = patch.faceCells();

            const vectorField internalCoeffs(pU.valueInternalCoeffs(pw));

            // Diag contribution
            forAll(pU, faceI)
            {
                diag[fc[faceI]] -=
                    sign*internalCoeffs[faceI][cmptI]*Sf[faceI][cmptI];
            }

            if (patch.coupled())
            {
                // Todo: add off-core coeffs
                notImplemented("patch.coupled(): D-in-p");

                // CoeffField<vector>::linearTypeField& pcoupleUpper =
                //     bs.coupleUpper()[patchI].asLinear();
                // CoeffField<vector>::linearTypeField& pcoupleLower =
                //     bs.coupleLower()[patchI].asLinear();

                // const vectorField pcl = -pw*Sf;
                // const vectorField pcu = pcl + Sf;

                // // Coupling  contributions
                // pcoupleLower -= pcl;
                // pcoupleUpper -= pcu;
            }
        }

        // Insert component coeffs
        // cmptI is the 1st column for Ux, 2nd for Uy, 3rd for Uz
        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
        (
            divUCoeffs, jac, rowOffset, colOffset + cmptI, 1
        );
    }

    return 0;
}


int foamPetscSnesHelper::solve(const bool returnOnSnesError)
{
    if (!initialised_)
    {
        FatalErrorIn("void foamPetscSnesHelper::solve()")
            << "This function cannot be called because the foamPetscSnesHelper "
            << "object was not initialised during construction"
            << abort(FatalError);
    }

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

    // Load the correct options database
    PetscOptionsPush(options_);
    SNESSetFromOptions(snes_.s);

    // Set the snesHasRun flag
    snesHasRun_ = true;

    // Solve the nonlinear system
    AssertPETSc(SNESSolve(snes_.s, NULL, x_.v));

    // Un-load the options file
    PetscOptionsPop();

    // Check convergence
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes_.s, &reason);

    if (reason == SNES_DIVERGED_FUNCTION_DOMAIN)
    {
        Info<< nl
            << "PETSc SNES solver returned a diverged function domain: "
            << "returning" << endl;

        return 0;
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
