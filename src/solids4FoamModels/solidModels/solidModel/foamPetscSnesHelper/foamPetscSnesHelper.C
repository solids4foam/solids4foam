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

#ifdef USE_PETSC

#include "foamPetscSnesHelper.H"
#include "SparseMatrixTemplate.H"
#include "processorFvPatch.H"

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

    // Access x and f data
    CHKERRQ(VecGetArrayRead(x, &xx));
    CHKERRQ(VecGetArray(f, &ff));

    // Compute the residual
    if (user->solMod_.formResidual(ff, xx) != 0)
    {
        Foam::FatalError
            << "formResidual(ff, xx) returned an error code!"
            << Foam::abort(Foam::FatalError);
    }

    // Restore the solution and residual vectors
    CHKERRQ(VecRestoreArrayRead(x, &xx));
    CHKERRQ(VecRestoreArray(f, &ff));

    return 0;
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
    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Map the solution to an OpenFOAM field
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;
    const bool twoD = user->solMod_.twoD();
    const int blockSize = twoD ? 2 : 3;

    // Initialise the matrix if it has yet to be allocated; otherwise zero all
    // entries
    MatInfo info;
    MatGetInfo(B, MAT_LOCAL, &info);
    // PetscInt rstart, rend;
    // MatGetOwnershipRange(B, &rstart, &rend);
    if (info.nz_used)
    {
        // Zero the matrix but do not reallocate the space
        // The "-snes_lag_jacobian -2" PETSc option can be used to avoid
        // re-building the matrix
        CHKERRQ(MatZeroEntries(B));
    }
    else
    {
        Foam::Info<< "    Initialising the matrix..." << Foam::flush;

        // Set the block size
        CHKERRQ(MatSetBlockSize(B, blockSize));

        // Number of vector unknowns on this processor
        const int blockn = user->solMod_.globalCells().localSize();

        // Count the number of non-zeros in the matrix
        // Note: we assume a compact stencil, i.e. face only face neighbours

        // Number of on-processor non-zeros per row
        int* d_nnz = (int*)malloc(blockn*sizeof(int));

        // Number of off-processor non-zeros per row
        int* o_nnz = (int*)malloc(blockn*sizeof(int));

        // Initialise d_nnz and o_nnz to zero
        for (int i = 0; i < blockn; ++i)
        {
            d_nnz[i] = 1; // count diagonal cell
            o_nnz[i] = 0;
        }

        // Take a reference to the mesh
        const Foam::fvMesh& mesh = user->solMod_.solution().mesh();

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

        // Allocate parallel matrix
        //CHKERRQ(MatMPIAIJSetPreallocation(B, 0, d_nnz, 0, o_nnz));
        // Allocate parallel matrix with the same conservative stencil per node
        //CHKERRQ(MatMPIAIJSetPreallocation(B, d_nz, NULL, 0, NULL));
        CHKERRQ(MatMPIBAIJSetPreallocation(B, blockSize, 0, d_nnz, 0, o_nnz));

        // Raise an error if mallocs are required during matrix assembly
        MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

        Foam::Info<< "done" << Foam::endl;
    }

    // Populate the Jacobian => implemented by the solid model
    if (user->solMod_.formJacobian(B, xx) != 0)
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(foamPetscSnesHelper, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamPetscSnesHelper::foamPetscSnesHelper
(
    fileName optionsFile,
    volVectorField& solution,
    const Switch twoD,
    const Switch stopOnPetscError,
    const Switch initialise
)
:
    initialised_(initialise),
    solution_(solution),
    globalCells_(solution.mesh().nCells()),
    blockSize_(twoD ? 2 : 3),
    twoD_(twoD),
    stopOnPetscError_(stopOnPetscError),
    snesPtr_(),
    xPtr_(),
    APtr_(),
    snesUserPtr_()
{
    if (initialise)
    {
        // Expand the options file name
        optionsFile.expand();

        // Initialise PETSc with an options file
        PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

        // Create the PETSc SNES object
        snesPtr_.set(new SNES());
        SNES& snes = snesPtr_();
        SNESCreate(PETSC_COMM_WORLD, &snes);

        // Create user data context
        snesUserPtr_.set(new appCtxfoamPetscSnesHelper(*this));
        appCtxfoamPetscSnesHelper& user = snesUserPtr_();

        // Set the user context
        SNESSetApplicationContext(snes, &user);

        // Set the residual function
        SNESSetFunction(snes, NULL, formResidualFoamPetscSnesHelper, &user);

        // Coefficient block size, e.g. 2x2 in 2-D, 3x3 in 3-D
        const int blockSize = twoD_ ? 2 : 3;

        // Global system size
        const int blockN = globalCells_.size();
        const int N = blockSize*blockN;

        // Local (this processor) system size
        const int blockn = globalCells_.localSize();
        const int n = blockSize*blockn;

        // Create the Jacobian matrix
        APtr_.set(new Mat());
        Mat& A = APtr_();
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, n, n, N, N);
        MatSetFromOptions(A);
        MatSetType(A, MATMPIAIJ);
        //MatSetType(A, MATMPIBAIJ);

        // Set the Jacobian function
        SNESSetJacobian(snes, A, A, formJacobianFoamPetscSnesHelper, &user);

        // Set solver options
        // Uses default options, can be overridden by command line options
        SNESSetFromOptions(snes);

        // Create the solution vector
        xPtr_.set(new Vec());
        Vec& x = xPtr_();
        VecCreate(PETSC_COMM_WORLD, &x);
        VecSetSizes(x, n, N);
        VecSetBlockSize(x, blockSize);
        VecSetType(x, VECMPI);
        PetscObjectSetName((PetscObject) x, "Solution");
        VecSetFromOptions(x);
        VecZeroEntries(x);
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

foamPetscSnesHelper::~foamPetscSnesHelper()
{
    if (initialised_)
    {
        SNESDestroy(&snesPtr_());
        snesPtr_.clear();

        VecDestroy(&xPtr_());
        xPtr_.clear();

        MatDestroy(&APtr_());
        APtr_.clear();

        snesUserPtr_.clear();

        PetscFinalize();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void foamPetscSnesHelper::mapSolutionFoamToPetsc()
{
    if (!initialised_)
    {
        FatalErrorIn("void foamPetscSnesHelper::mapSolutionFoamToPetsc()")
            << "This function cannot be called because the foamPetscSnesHelper "
            << "object was not initialised during construction"
            << abort(FatalError);
    }

    // Take a reference to the PETSc solution vector
    Vec& x = xPtr_();

    // Take a reference to the Foam solution field
    const volVectorField& f = solution_;

    PetscScalar *xx;
    VecGetArray(x, &xx);

    const vectorField& fI = f;

    int index = 0;

    forAll(fI, i)
    {
        xx[index++] = fI[i][vector::X];
        xx[index++] = fI[i][vector::Y];

        if (!twoD_)
        {
            xx[index++] = fI[i][vector::Z];
        }
    }

    VecRestoreArray(x, &xx);
}


void foamPetscSnesHelper::mapSolutionPetscToFoam()
{
    if (!initialised_)
    {
        FatalErrorIn("void foamPetscSnesHelper::mapSolutionPetscToFoam()")
            << "This function cannot be called because the foamPetscSnesHelper "
            << "object was not initialised during construction"
            << abort(FatalError);
    }

    // Take a reference to the PETSc solution vector
    const Vec& x = xPtr_();

    // Take a reference to the Foam solution field
    volVectorField& f = solution_;

    const PetscScalar *xx;
    VecGetArrayRead(x, &xx);

    vectorField& fI = f;

    int index = 0;

    forAll(fI, i)
    {
        fI[i].x() = xx[index++];
        fI[i].y() = xx[index++];

        if (!twoD_)
        {
            fI[i].z() = xx[index++];
        }
    }

    f.correctBoundaryConditions();

    VecRestoreArrayRead(x, &xx);
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

    // Take a reference to the SNES solver object
    SNES& snes = snesPtr_();

    // Take a reference to the solution vector
    Vec& x = xPtr_();

    // Solve the nonlinear system
    SNESSolve(snes, NULL, x);

    // Check convergence
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes, &reason);

    if (reason < 0)
    {
        WarningIn
        (
            "void foamPetscSnesHelper::checkConvergence(SNES snes) const"
        )   << "PETSc SNES solver return error check disabled" << endl
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
            FatalErrorIn
            (
                "void foamPetscSnesHelper::checkConvergence(SNES snes) const"
            )   << "Stopping because of the PETSc SNES error" << nl
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
