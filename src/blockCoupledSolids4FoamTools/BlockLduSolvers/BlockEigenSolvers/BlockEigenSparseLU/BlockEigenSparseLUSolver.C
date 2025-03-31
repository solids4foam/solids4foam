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

#ifndef S4F_NO_USE_EIGEN

#include "blockLduSolvers.H"
#include "BlockEigenSparseLUSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSparseLUSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenSparseLUSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenSparseLUSolver, asymMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenSparseLUSolver::BlockEigenSparseLUSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockEigenSolver(fieldName, matrix, dict)
{
    Info<< type() << " Eigen linear solver" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockEigenSparseLUSolver::solve
(
    Field<Foam::vector>& U,
    const Field<Foam::vector>& blockB
) const
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

    // Create scaling
    //Eigen::IterScaling<Eigen::SparseMatrix<double> >* scalPtr = NULL;

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": direct solution of sparse system" << endl;
    }

    // Now, solve the equilibrated linear system with any available solver
    // In this case we use a Direct sparse LU solver
    // Trial an error suggests COLAMD ordering is by far the best
    //Eigen::SimplicialLDLT
    Eigen::SparseLU
    <
        Eigen::SparseMatrix<scalar>, Eigen::COLAMDOrdering<int>
        //Eigen::SparseMatrix<scalar>, Eigen::AMDOrdering<int>
        //Eigen::SparseMatrix<scalar>, Eigen::MetisOrdering<int>
    > solver(A);

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x = solver.solve(b);


    // We copy the results from the std::vector into the geometric field
    // This can be costly
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": copying results into foam format" << nl
            << endl;
    }
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


    /*
    // Testing: eigenproblem
    {
        Info<< nl << "Computing the eigenvalues and eigenvectors" << endl;

        // Only for symmetric matrices; only uses the lower
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(A);

        Info<< "Copying one of the eigenvalues and eigenvectors" << endl;
        index = 0;
        forAll(U, cellI)
        {
            U[cellI].x() = es.eigenvectors().col(m - 3)(index++);
            U[cellI].y() = es.eigenvectors().col(m - 3)(index++);

            if (!twoD)
            {
                U[cellI].z() =
                    es.eigenvectors().col(m - 3)(index++);
            }
        }

        // These methods are much slower and I am not sure if I trust their
        // answers
        //Eigen::JacobiSVD<Eigen::MatrixXf> svd
        // Eigen::BDCSVD<Eigen::MatrixXf> svd
        // (
        //     A, Eigen::ComputeThinU | Eigen::ComputeThinV
        // );

        // Info<< "Copying the 2nd last eigenvalues and eigenvectors" << endl;
        // index = 0;
        // forAll(U, cellI)
        // {
        //     U[cellI].x() = svd.matrixU()(index++, m - 3);
        //     U[cellI].y() = svd.matrixU()(index++, m - 3);

        //     if (!twoD)
        //     {
        //         U[cellI].z() =
        //             svd.matrixU()(index++, m - 3);
        //     }
        // }
    }*/


    return BlockSolverPerformance<Foam::vector>
    (
        this->typeName,
        this->fieldName(),
        pTraits<Foam::vector>::one,
        pTraits<Foam::vector>::zero,
        0,
        true,
        false
    );
}

#endif // end of #ifndef S4F_NO_USE_EIGEN

// ************************************************************************* //
