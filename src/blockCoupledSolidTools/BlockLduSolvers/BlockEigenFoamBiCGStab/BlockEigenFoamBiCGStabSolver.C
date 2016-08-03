/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blockLduSolvers.H"
#include "BlockEigenFoamBiCGStabSolver.H"
#include "BlockSolverPerformance.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "OFstream.H"
#include "polyMesh.H"

// Eigen files
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/IterativeLinearSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenFoamBiCGStabSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenFoamBiCGStabSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenFoamBiCGStabSolver, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenFoamBiCGStabSolver::BlockEigenFoamBiCGStabSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector>(fieldName, matrix, dict),
    preconPtr_
    (
        BlockLduPrecon<vector>::New
        (
            matrix,
            this->dict()
        )
    ),
    tol_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    maxIters_(dict.lookupOrDefault<scalar>("maxIter", 1000))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockEigenFoamBiCGStabSolver::solve
(
    Field<Foam::vector>& x,
    const Field<Foam::vector>& rhs
)
{
    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<vector>& matrix = this->matrix_;

    // Prepare solver performance
    BlockSolverPerformance<vector> solverPerf
    (
        typeName,
        this->fieldName()
    );

    // Grab matrix diagonal and off-diagonals
    const Field<tensor>& d = matrix.diag().asSquare();

    // Solve the system

    Info<< this->typeName << ": iterative solution of sparse system" << endl;

    // We will set an initial guess of x = rhs, like in Eigen BiCGStab
    // This works but why does this work...
    // Using a null vector diverges in some cases whereas x = rhs will converge.
    //x = rhs;
    x = gAverage(rhs);
    Info<< "gAverage(rhs) is " << gAverage(rhs) << endl;
    //scalar magAv = mag(gAverage(rhs));
    //x = magAv*vector(1, 1, 0);


    // Taken and modified from the Eigen header file BiCGStab.h

    scalar tol = tol_;
    int maxIters = maxIters_;

    int n = d.size();
    //x = precond.solve(x);
    vectorField xTemp = x;
    preconPtr_->precondition(x, xTemp);

    //VectorType r  = rhs - mat*x;
    vectorField r(n);
    matrix.Amul(r, x);
    r = rhs - r;
    vectorField r0 = r;

    // scalar r0_sqnorm = r0.squaredNorm();
    // scalar rhs_sqnorm = rhs.squaredNorm();
    scalar r0_sqnorm = gSum(magSqr(r0));
    scalar rhs_sqnorm = gSum(magSqr(rhs));

    solverPerf.initialResidual() = r0_sqnorm*vector::one;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // if (rhs_sqnorm == 0)
    if (rhs_sqnorm < SMALL)
    //do while (!this->stop(solverPerf))
    {
        //x.setZero();
        x = vector::zero;
        //return true;
        // TODO: return solverPerf
        FatalErrorIn("solve")
            << "rhs_sqnorm is SMALL: " << rhs_sqnorm << abort(FatalError);
    }

    scalar rho = 1.0;
    scalar alpha = 1.0;
    scalar w = 1.0;

    //VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
    vectorField v(n, vector::zero);
    vectorField p(n, vector::zero);
    vectorField y(n),  z(n);
    vectorField kt(n), ks(n);

    // VectorType s(n), t(n);
    vectorField s(n);
    vectorField t(n);

    scalar tol2 = tol*tol;
//  RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
//     scalar eps2 = VSMALL*VSMALL;
    scalar eps2 = SMALL*SMALL;
    int i = 0;
    int restarts = 0;

    // while (r.squaredNorm()/rhs_sqnorm > tol2 && i < maxIters)
    while (gSum(magSqr(r))/rhs_sqnorm > tol2 && i < maxIters)
    {
        scalar rho_old = rho;

        //rho = r0.dot(r);
        rho = sum(r0 & r);

        if (mag(rho) < eps2*r0_sqnorm)
        {
            // The new residual vector became too orthogonal to the arbitrarily
            // choosen direction r0
            // Let's restart with a new r0:
            r0 = r;
            // rho = r0_sqnorm = r.squaredNorm();
            r0_sqnorm = sum(magSqr(r));
            rho = r0_sqnorm;
            if (restarts++ == 0)
            {
                i = 0;
            }
        }
        scalar beta = (rho/rho_old)*(alpha/w);

        p = r + beta*(p - w*v);

        //y = precond.solve(p);
        preconPtr_->precondition(y, p);

        //v.noalias() = mat*y;
        matrix.Amul(v, y);

        // alpha = rho/r0.dot(v);
        alpha = rho/sum(r0 & v);
        s = r - alpha*v;

        //z = precond.solve(s);
        preconPtr_->precondition(z, s);
        //t.noalias() = mat*z;
        matrix.Amul(t, z);

        // RealScalar tmp = t.squaredNorm();
        scalar tmp = sum(magSqr(t));
        //if (tmp > RealScalar(0))
        if (tmp > VSMALL)
        {
            // w = t.dot(s)/tmp;
            w = sum(t & s)/tmp;
        }
        else
        {
            // w = Scalar(0);
            w = 0.0;
        }

        x += alpha*y + w*z;
        r = s - w*t;
        ++i;

        Info<< "    iteration = " << i
            << ", residual = " << ::sqrt(gSum(magSqr(r))/rhs_sqnorm) << endl;
    }

    // scalar residual = ::sqrt(sum(magSqr(r))/rhs_sqnorm);
    // int iters = i;

    solverPerf.finalResidual() = ::sqrt(gSum(magSqr(r))/rhs_sqnorm)*vector::one;
    solverPerf.nIterations() = i;

    // Info<< "    solver time = "
    //     << (matrix.mesh().thisDb().time().elapsedCpuTime() - startT) << endl;

    return solverPerf;
}


// ************************************************************************* //
