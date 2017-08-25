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
#include "BlockEigenProblemSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenProblemSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenProblemSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenProblemSolver, asymMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenProblemSolver::BlockEigenProblemSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockEigenSolver(fieldName, matrix, dict),
    eigenvaluesPtr_(),
    eigenvectorsPtr_()
{
    Info<< type() << " : Eigen linear solver" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockEigenProblemSolver::solve
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
            "Foam::BlockEigenProblemSolver::solve"
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

    // Solve eigen problem only once and then store the eigenvalues and
    // eigenvectors
    if (!eigenvaluesPtr_.valid())
    {
        Info<< "Converting the matrix from foam format to Eigen format" << endl;

        // Create Eigen sparse matrix and set coeffs
        Eigen::SparseMatrix<scalar> A(m, m);

        // Convert foam matrix to the Eigen matrix format
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

        // Testing
        {
            // Flip sign so diagonal is positive
            Warning
                << "flipping the sign of the matrix!" << endl;
            A *= -1;
            b *= -1;
        }

        // Optionally export system to Matlab
        if (writeMatlabFiles())
        {
            writeLinearSystemToMatlabFiles(A, b);
        }

        // Compressing matrix is meant to help performance
        A.makeCompressed();

        // Solve the eigen problem
        Info<< "Computing the eigenvalues and eigenvectors" << endl;

        // Only for symmetric matrices; only uses the lower
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(A);

        /*
        // Create geometric/stress/inertia/loading stiffness matrix
        Eigen::SparseMatrix<scalar> B(m, m);
        // Insert coefficients
        {
            std::vector< Eigen::Triplet<scalar> > coefficients;
            coefficients.reserve(m);

            // Insert the density on the diagonal
            // TODO: we need rho and nCells from the solidModels to create the
            // loading matrix
            //const scalar rho = 8000;
            // TEST there are 500 cells in the 2-D test case
            // We do not add density to the boundary condition equations
            //for (int ID = 0; ID < 1000; ID++)
            for (int ID = 0; ID < A.rows(); ID++)
            {
                // if (twoD)
                // {
                //     // 2-D
                //     ID = 2*dI;
                // }
                // else
                // {
                //     // 3-D
                //     ID = 3*dI;
                // }
                // if (ID < 1000) // (A.rows() - 1))
                // {
                //     coefficients.push_back
                //     (
                //         Eigen::Triplet<scalar>(ID, ID, 1e-6*rho)
                //     );
                // }
                // else
                //{
                    //coefficients.push_back
                    //(
                    //    Eigen::Triplet<scalar>(ID, ID, rho)
                    //);
                coefficients.push_back(Eigen::Triplet<scalar>(ID, ID, 7800.0));
                //}
            }

            // Insert triplets into the matrix
            B.setFromTriplets(coefficients.begin(), coefficients.end());
        }
        // Create the eigen problem solver
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> es(A, B);
            */

        // Store the eigenvalues
        eigenvaluesPtr_.set
        (
            new Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf>::RealVectorType
            (
                es.eigenvalues()
            )
        );

        // Store the eigenvectors
        eigenvectorsPtr_.set
        (
            new Eigen::MatrixXf
            (
                es.eigenvectors()
            )
        );
    }

    // Return the lowest mode
    // Further modes can be retrieved through the mode function

    const label modeI = 1;

    // Copy the results into the foam format
    Info<< "Copying eigenvalue and eigenvectors for mode: " << modeI << endl;

    Info<< "Eigenvalue for mode " << modeI << ": "
        << eigenvaluesPtr_()(m - modeI) << endl;

    label index = 0;

    forAll(U, cellI)
    {
        U[cellI].x() = eigenvectorsPtr_().col(m - modeI)(index++);
        U[cellI].y() = eigenvectorsPtr_().col(m - modeI)(index++);

        if (!twoD)
        {
            U[cellI].z() = eigenvectorsPtr_().col(m - modeI)(index++);
        }
    }

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


void Foam::BlockEigenProblemSolver::mode
(
    const int i,
    Field<vector>& U
) const
{
    if (!eigenvaluesPtr_.valid())
    {
        FatalErrorIn
        (
            "void Foam::BlockEigenProblemSolver::mode\n"
            "(\n"
            "    const int modeI,\n"
            "    Field<vector>& x\n"
            ") const"
        )   << "The solve(...) function must be called to calculate the modes "
            << "before they can be retrieved!"
            << abort(FatalError);
    }

    const int modeI = i + 1;

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate te number of degrees of freedom
    const int m = calcDegreesOfFreedom(matrix_, twoD);

    Info<< "Mode " << modeI << nl
        << "    eigenvalue: "
        //<< eigenvaluesPtr_()(m - modeI) << nl
        << eigenvaluesPtr_()(modeI - 1)
        << "    frequency: "
        //<< sqrt(eigenvaluesPtr_()(m - modeI))/mathematicalConstant::twoPi
        //<< sqrt(eigenvaluesPtr_()(modeI - 1))/mathematicalConstant::twoPi
        << " Hz" << endl;

    // Testing
    for (int i = 0; i < 10; i++)
    {
        Info<< "eigval[" << i << "] " << eigenvaluesPtr_()(i) << endl;
    }

    label index = 0;

    forAll(U, cellI)
    {
        // U[cellI].x() = eigenvectorsPtr_().col(m - modeI)(index++);
        // U[cellI].y() = eigenvectorsPtr_().col(m - modeI)(index++);
        U[cellI].x() = eigenvectorsPtr_().col(modeI - 1)(index++);
        U[cellI].y() = eigenvectorsPtr_().col(modeI - 1)(index++);

        if (!twoD)
        {
            //U[cellI].z() = eigenvectorsPtr_().col(m - modeI)(index++);
            U[cellI].z() = eigenvectorsPtr_().col(modeI - 1)(index++);
        }
    }
}


// ************************************************************************* //
