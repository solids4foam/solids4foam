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
#include "BlockEigenSolver.H"
#include "OFstream.H"
#include "polyMesh.H"
#include <unsupported/Eigen/SparseExtra>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolver, 0);
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

bool Foam::BlockEigenSolver::checkTwoD() const
{
    // Check if mesh is 2-D or 3-D
    // This should be OK for multi-region cases as all regions should be 2-D in
    // same direction
    const Vector<label>& geomD =
        matrix_.mesh().thisDb().parent().lookupObject<polyMesh>
        (
            "region0" // polyMesh::defaultRegion,
        ).geometricD();

    label nDir = 0;
    forAll(geomD, dirI)
    {
        if (geomD[dirI] > SMALL)
        {
            nDir++;
        }
    }

    bool twoD = false;
    if (nDir == 2)
    {
        if
        (
            matrix_.mesh().thisDb().parent().lookupObject<polyMesh>
            (
                "region0" // polyMesh::defaultRegion,
            ).solutionD()[vector::Z]
          > -1
        )
        {
            FatalErrorIn
            (
                "Foam::BlockSolverPerformance<Foam::vector>"
                "Foam::BlockEigenSolver::solve"
                "("
                "    Field<Foam::vector>& U,"
                "    const Field<Foam::vector>& blockB"
                ")"
            )   << this->typeName << ": for 2-D models, the empty direction "
                << "must be z!" << abort(FatalError);
        }

        twoD = true;
    }
    else if (nDir != 3)
    {
        FatalErrorIn(this->typeName + " solve")
            << "solver only implemented for 2-D and 3-D models"
            << abort(FatalError);
    }

    return twoD;
}


int Foam::BlockEigenSolver::calcDegreesOfFreedom
(
    const BlockLduMatrix<vector>& matrix,
    const bool twoD
) const
{
    // Number of vectors on the diagonal
    const int diagSize = matrix.diag().asSquare().size();

    if (twoD)
    {
        // Two components in each vector
        return 2*diagSize;
    }

    // Three components in each vector
    return 3*diagSize;
}


void Foam::BlockEigenSolver::convertFoamMatrixToEigenMatrix
(
    const BlockLduMatrix<vector>& matrix,
    Eigen::SparseMatrix<scalar>& A
)
{
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into Eigen format"
            << endl;
    }

    // Matrix addressing
    const unallocLabelList& lowerAddr = matrix.mesh().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix.mesh().lduAddr().upperAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor>& d = matrix.diag().asSquare();
    const Field<tensor>& l = matrix.lower().asSquare();
    const Field<tensor>& u = matrix.upper().asSquare();

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate te number of degrees of freedom
    const int m = calcDegreesOfFreedom(matrix, twoD);

    // Create coefficient matrix: we must copy coeffs from ldu storage
    // This can be costly
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(m);

    // Insert coefficients

    // Insert diagonal
    label ID = -1;
    forAll(d, dI)
    {
        if (twoD)
        {
            // 2-D
            ID = 2*dI;
        }
        else
        {
            // 3-D
            ID = 3*dI;
        }

        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID, d[dI].xx()));
        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID+1, d[dI].xy()));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID, d[dI].yx()));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID+1, d[dI].yy()));

        if (!twoD)
        {
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(ID, ID+2, d[dI].xz())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(ID+1, ID+2, d[dI].yz())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(ID+2, ID, d[dI].zx())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(ID+2, ID+1, d[dI].zy())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(ID+2, ID+2, d[dI].zz())
            );
        }
    }

    // Insert off-diagonal
    label rowI = -1;
    label columnI = -1;
    forAll(u, uI)
    {
        const label own = lowerAddr[uI];
        const label nei = upperAddr[uI];

        const tensor& upper = u[uI];
        const tensor& lower = l[uI];

        if (twoD)
        {
            rowI = 2*own;
            columnI = 2*nei;
        }
        else
        {
            rowI = 3*own;
            columnI = 3*nei;
        }

        // Upper
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI, upper.xx())
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI+1, upper.xy())
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI, upper.yx())
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI+1, upper.yy())
        );

        // Lower
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI, lower.xx())
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI+1, lower.xy())
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI, lower.yx())
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI+1, lower.yy())
        );

        if (!twoD)
        {
            // Upper
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(rowI, columnI+2, upper.xz())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(rowI+1, columnI+2, upper.yz())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(rowI+2, columnI, upper.zx())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(rowI+2, columnI+1, upper.zy())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(rowI+2, columnI+2, upper.zz())
            );

            // Lower
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(columnI, rowI+2, lower.xz())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(columnI+1, rowI+2, lower.yz())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(columnI+2, rowI, lower.zx())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(columnI+2, rowI+1, lower.zy())
            );
            coefficients.push_back
            (
                Eigen::Triplet<scalar>(columnI+2, rowI+2, lower.zz())
            );
        }
    }

    // Insert triplets into the matrix
    A.setFromTriplets(coefficients.begin(), coefficients.end());
}

void Foam::BlockEigenSolver::writeLinearSystemToMatlabFiles
(
    const Eigen::SparseMatrix<scalar>& A,
    const Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b
) const
{
    if (writeMatlabFiles_)
    {
        Info<< type() << ": writing sparse matrix to "
            << "matlabSparseMatrix.txt and source to matlabSource.txt"
            << nl
            << "The system can be read and solved in Matlab or Octave with "
            << "commands:" << nl
            << nl
            << "    load matlabSparseMatrix.txt;" << nl
            << "    A = spconvert(matlabSparseMatrix);" << nl
            << "    B = dlmread('matlabSource.txt', ' ');" << nl
            << "    x = A\\B;" << nl
            << endl;

        // Write matrix
        Eigen::saveMarket(A, "matlabSparseMatrix.txt");

        // Write source
        OFstream sourceFile("matlabSource.txt");
        for (int rowI = 0; rowI < A.rows(); rowI++)
        {
            sourceFile
                << b(rowI) << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenSolver::BlockEigenSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector>(fieldName, matrix, dict),
    writeMatlabFiles_
    (
        dict.lookupOrDefault<Switch>("writeMatlabFiles", false)
    )
{
    Info<< type() << " : writeMatlabFiles is " << writeMatlabFiles_ << endl;
}


// ************************************************************************* //
