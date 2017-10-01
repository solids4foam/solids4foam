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
#include "BlockEigenGMResSolver.H"
#include "addToRunTimeSelectionTable.H"
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenGMResSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenGMResSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenGMResSolver, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenGMResSolver::BlockEigenGMResSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockEigenSolver(fieldName, matrix, dict),
    tol_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    restart_(dict.lookupOrDefault<int>("restart", 5))
{
    Info<< type() << " Eigen linear solver" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockEigenGMResSolver::solve
(
    Field<Foam::vector>& U,
    const Field<Foam::vector>& blockB
)
{
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector>"
            "Foam::BlockEigenSparseLUSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "SparesLU direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate te number of degrees of freedom
    const int m = calcDegreesOfFreedom(matrix_, twoD);

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A(m, m);

    // Convert foam matrix to the Eigen matrix format
    //convertFoamMatrixToEigenMatrix(d, u, l, upperAddr, lowerAddr, twoD, A);
    convertFoamMatrixToEigenMatrix(matrix_, A);

    // Copy source vector into Eigen vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(m);
    label index = 0;
    forAll(blockB, rowI)
    {
        b(index++) = blockB[rowI].x();
        b(index++) = blockB[rowI].y();

        if (!twoD)
        {
            b(index++) = blockB[rowI].z();
        }
    }

    // Optionally export system to Matlab
    if (writeMatlabFiles())
    {
        writeLinearSystemToMatlabFiles(A, b);
    }

    // Compressing matrix is meant to help performance
    A.makeCompressed();

    // Solve the system

    Info<< this->typeName << ": iterative solution of sparse system" << endl;

    // Generalised mimimal residual method iterative solver
    Eigen::GMRES
    <
        //Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>
        Eigen::SparseMatrix<double>
    > solver(A);

    solver.setTolerance(tol_);
    if (restart_ > 0)
    {
        solver.set_restart(restart_);
    }

    Info<< "    Using a restart frequency of " << solver.get_restart() << endl;

    // solver.preconditioner().setDroptol(0.0001);
    // solver.preconditioner().setFillfactor(100);

    // Solve the system
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x = solver.solve(b);

    // We copy the results from the std::vector into the geometric field
    // This can be costly
    Info<< this->typeName << ": copying results into foam format" << nl << endl;
    index = 0;
    forAll(U, cellI)
    {
        U[cellI].x() = x(index++);
        U[cellI].y() = x(index++);

        if (!twoD)
        {
            U[cellI].z() = x(index++);
        }
    }

    return BlockSolverPerformance<Foam::vector>
    (
        this->typeName,
        this->fieldName(),
        vector::one,
        solver.error()*vector::one,
        solver.iterations(),
        bool(solver.error() < tol_),
        false
    );
}


// ************************************************************************* //
