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

#include "colouredFiniteDifferenceJacobian.H"
#include "petscErrorHandling.H"
#include "clockTime.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

PetscErrorCode Foam::colouredFiniteDifferenceJacobian::residual
(
    void*,
    Vec x,
    Vec f,
    void* ctx
)
{
    foamPetscSnesHelper& system = *static_cast<foamPetscSnesHelper*>(ctx);

    if (system.formResidual(f, x) != 0)
    {
        return PETSC_ERR_LIB;
    }

    return PETSC_SUCCESS;
}


void Foam::colouredFiniteDifferenceJacobian::setUp(Mat compact)
{
    // The residual of a cell depends on the cells within two layers: its
    // face values use the cell gradients of the two cells sharing the face.
    // The pattern of the square of the compact (one-layer) matrix is exactly
    // that two-layer stencil, in serial and in parallel.
    AssertPETSc
    (
        MatMatMult(compact, compact, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &J_)
    );
    AssertPETSc(MatZeroEntries(J_));

    // Distance-two colouring: columns of one colour share no row
    MatColoring mc;
    ISColoring isColouring;
    AssertPETSc(MatColoringCreate(J_, &mc));
    AssertPETSc
    (
        PetscObjectSetOptionsPrefix(PetscObject(mc), "solid_fd_")
    );
    AssertPETSc(MatColoringSetDistance(mc, 2));
    AssertPETSc(MatColoringSetType(mc, MATCOLORINGGREEDY));
    AssertPETSc(MatColoringSetFromOptions(mc));
    AssertPETSc(MatColoringApply(mc, &isColouring));
    AssertPETSc(MatColoringDestroy(&mc));

    AssertPETSc(MatFDColoringCreate(J_, isColouring, &fdc_));
    AssertPETSc(MatFDColoringSetFunction(fdc_, residual, &system_));
    AssertPETSc
    (
        PetscObjectSetOptionsPrefix(PetscObject(fdc_), "solid_fd_")
    );
    AssertPETSc(MatFDColoringSetFromOptions(fdc_));
    AssertPETSc(MatFDColoringSetUp(J_, isColouring, fdc_));

    PetscInt nColours = 0;
    AssertPETSc
    (
        ISColoringGetIS(isColouring, PETSC_USE_POINTER, &nColours, nullptr)
    );
    AssertPETSc(ISColoringRestoreIS(isColouring, PETSC_USE_POINTER, nullptr));
    AssertPETSc(ISColoringDestroy(&isColouring));
    nColours_ = nColours;

    MatInfo info;
    AssertPETSc(MatGetInfo(J_, MAT_GLOBAL_SUM, &info));

    Info<< "    " << name_ << ": coloured finite-difference Jacobian with "
        << nColours_ << " colours and " << label(info.nz_used)
        << " non-zeros" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::colouredFiniteDifferenceJacobian::colouredFiniteDifferenceJacobian
(
    const word& name,
    foamPetscSnesHelper& system
)
:
    name_(name),
    system_(system),
    J_(nullptr),
    fdc_(nullptr),
    nColours_(0),
    nBuilds_(0),
    buildTime_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::colouredFiniteDifferenceJacobian::~colouredFiniteDifferenceJacobian()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::colouredFiniteDifferenceJacobian::clear()
{
    if (fdc_)
    {
        MatFDColoringDestroy(&fdc_);
        fdc_ = nullptr;
    }
    if (J_)
    {
        MatDestroy(&J_);
        J_ = nullptr;
    }
    nColours_ = 0;
    nBuilds_ = 0;
}


void Foam::colouredFiniteDifferenceJacobian::apply
(
    Mat A,
    const Vec x,
    const bool rebuild
)
{
    // Complete the compact matrix (its entries may have been added)
    AssertPETSc(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    AssertPETSc(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    if (!J_)
    {
        setUp(A);
    }

    if (rebuild || nBuilds_ == 0)
    {
        // Wall-clock time, as the solve time it is compared with by the
        // automatic selection of the solid block
        clockTime timer;

        AssertPETSc(MatFDColoringApply(J_, fdc_, x, nullptr));

        // Leave the fields of the system at the unperturbed state
        Vec f;
        AssertPETSc(VecDuplicate(x, &f));
        if (system_.formResidual(f, x) != 0)
        {
            FatalErrorInFunction
                << name_ << ": residual evaluation failed" << abort(FatalError);
        }
        AssertPETSc(VecDestroy(&f));

        ++nBuilds_;
        const scalar t = timer.elapsedTime();
        buildTime_ += t;

        Info<< "    " << name_ << ": finite-difference Jacobian built ("
            << nColours_ << " colours, " << t << " s, build " << nBuilds_
            << ", total " << buildTime_ << " s)" << endl;
    }

    // Replace the compact entries; the pattern of A grows on the first call
    // and is then fixed
    AssertPETSc(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    AssertPETSc(MatCopy(J_, A, DIFFERENT_NONZERO_PATTERN));
}


Foam::scalar Foam::colouredFiniteDifferenceJacobian::check
(
    const Vec x,
    const scalar relStep
)
{
    if (!J_ || nBuilds_ == 0)
    {
        return -1;
    }

    Vec v, Jv, f0, f1, xp;
    AssertPETSc(VecDuplicate(x, &v));
    AssertPETSc(VecDuplicate(x, &Jv));
    AssertPETSc(VecDuplicate(x, &f0));
    AssertPETSc(VecDuplicate(x, &f1));
    AssertPETSc(VecDuplicate(x, &xp));

    // Random direction of unit norm, reproducible
    PetscRandom rnd;
    AssertPETSc(PetscRandomCreate(PetscObjectComm(PetscObject(x)), &rnd));
    AssertPETSc(PetscRandomSetSeed(rnd, 1234));
    AssertPETSc(PetscRandomSeed(rnd));
    AssertPETSc(VecSetRandom(v, rnd));
    AssertPETSc(PetscRandomDestroy(&rnd));
    AssertPETSc(VecShift(v, -0.5));
    PetscReal nv = 0, nx = 0;
    AssertPETSc(VecNorm(v, NORM_2, &nv));
    AssertPETSc(VecScale(v, 1.0/nv));
    AssertPETSc(VecNorm(x, NORM_2, &nx));

    // Forward difference of the residual along v
    const PetscReal h = relStep*(1.0 + nx);
    if (system_.formResidual(f0, x) != 0)
    {
        FatalErrorInFunction
            << name_ << ": residual evaluation failed" << abort(FatalError);
    }
    AssertPETSc(VecWAXPY(xp, h, v, x));
    if (system_.formResidual(f1, xp) != 0)
    {
        FatalErrorInFunction
            << name_ << ": residual evaluation failed" << abort(FatalError);
    }
    AssertPETSc(VecAXPY(f1, -1.0, f0));
    AssertPETSc(VecScale(f1, 1.0/h));

    // J v - FD
    AssertPETSc(MatMult(J_, v, Jv));
    PetscReal nRef = 0, nErr = 0;
    AssertPETSc(VecNorm(f1, NORM_2, &nRef));
    AssertPETSc(VecAXPY(Jv, -1.0, f1));
    AssertPETSc(VecNorm(Jv, NORM_2, &nErr));

    // Restore the state
    if (system_.formResidual(f0, x) != 0)
    {
        FatalErrorInFunction
            << name_ << ": residual evaluation failed" << abort(FatalError);
    }

    AssertPETSc(VecDestroy(&v));
    AssertPETSc(VecDestroy(&Jv));
    AssertPETSc(VecDestroy(&f0));
    AssertPETSc(VecDestroy(&f1));
    AssertPETSc(VecDestroy(&xp));

    const scalar relErr = (nRef > 0) ? nErr/nRef : 0;

    Info<< "    " << name_ << ": directional check |J v - D_h R v|/|D_h R v| = "
        << relErr << endl;

    return relErr;
}


#endif // USE_PETSC

// ************************************************************************* //
