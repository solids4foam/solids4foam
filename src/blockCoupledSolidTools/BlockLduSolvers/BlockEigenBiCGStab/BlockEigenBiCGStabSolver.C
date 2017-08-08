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
#include "BlockEigenBiCGStabSolver.H"
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
    defineTypeNameAndDebug(BlockEigenBiCGStabSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenBiCGStabSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockEigenBiCGStabSolver, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenBiCGStabSolver::BlockEigenBiCGStabSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector>(fieldName, matrix, dict),
    exportSystemMatlab_
    (
        dict.lookupOrDefault<Switch>("exportSystemMatlab", false)
    ),
    tol_(dict.lookupOrDefault<scalar>("tolerance", 1e-6)),
    maxIter_(dict.lookupOrDefault<scalar>("maxIter", 1000))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockEigenBiCGStabSolver::solve
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
            "Foam::BlockEigenBiCGStabSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "This linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Tolerance for the inclusion of coefficients
    const scalar tol = VSMALL;

    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<vector>& matrix = this->matrix_;

    // Matrix addressing
    const unallocLabelList& lowerAddr = matrix.mesh().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix.mesh().lduAddr().upperAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor>& d = matrix.diag().asSquare();
    const Field<tensor>& l = matrix.lower().asSquare();
    const Field<tensor>& u = matrix.upper().asSquare();

    // Check if mesh is 2-D or 3-D
    // We lookup the mesh by name: we need a better way for multiple regions
    // models; I don't think it is possible from just examining the matrix and
    // source as the discretisation is performed in 3-D
    // We could set diag empty components to zero in the solver and then check
    // here... maybe
    const Vector<label>& geomD =
        matrix.mesh().thisDb().parent().lookupObject<polyMesh>
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
        Warning
            << this->typeName << ": for 2-D models, it is assumed that z is the"
            << " empty direction" << endl;
        twoD = true;
    }
    else if (nDir != 3)
    {
        FatalErrorIn(this->typeName + " solve")
            << "solver only implemented for 2-D and 3-D models"
            << abort(FatalError);
    }

    Info<< this->typeName << ": copying matrix coefficients into Eigen format"
        << endl;

    // Degrees of freedom
    int m = 0;
    if (twoD)
    {
        // 2-D
        m = 2*d.size();
    }
    else
    {
        // 3-D
        m = 3*d.size();
    }

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

        if (mag(d[dI].xx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID, ID, d[dI].xx())
                );
        }

        if (mag(d[dI].xy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID, ID+1, d[dI].xy())
                );
        }

        if (mag(d[dI].yx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID+1, ID, d[dI].yx())
                );
        }

        if (mag(d[dI].yy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID+1, ID+1, d[dI].yy())
                );
        }

        if (!twoD)
        {
            if (mag(d[dI].xz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID, ID+2, d[dI].xz())
                    );
            }

            if (mag(d[dI].yz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+1, ID+2, d[dI].yz())
                    );
            }

            if (mag(d[dI].zx()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+2, ID, d[dI].zx())
                    );
            }

            if (mag(d[dI].zy()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+2, ID+1, d[dI].zy())
                    );
            }

            if (mag(d[dI].zz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+2, ID+2, d[dI].zz())
                    );
            }
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

        if (mag(upper.xx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI, columnI, upper.xx())
                );
        }

        if (mag(upper.xy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI, columnI+1, upper.xy())
                );
        }

        if (mag(upper.yx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI+1, columnI, upper.yx())
                );
        }

        if (mag(upper.yy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI+1, columnI+1, upper.yy())
                );
        }

        // Lower

        if (mag(lower.xx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI, rowI, lower.xx())
                );
        }

        if (mag(lower.xy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI, rowI+1, lower.xy())
                );
        }

        if (mag(lower.yx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI+1, rowI, lower.yx())
                );
        }

        if (mag(lower.yy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI+1, rowI+1, lower.yy())
                );
        }

        if (!twoD)
        {
            // Upper

            if (mag(upper.xz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI, columnI+2, upper.xz())
                    );
            }

            if (mag(upper.yz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+1, columnI+2, upper.yz())
                    );
            }

            if (mag(upper.zx()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+2, columnI, upper.zx())
                    );
            }

            if (mag(upper.zy()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+2, columnI+1, upper.zy())
                    );
            }

            if (mag(upper.zz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+2, columnI+2, upper.zz())
                    );
            }

            // Lower

            if (mag(lower.xz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI, rowI+2, lower.xz())
                    );
            }

            if (mag(lower.yz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+1, rowI+2, lower.yz())
                    );
            }

            if (mag(lower.zx()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+2, rowI, lower.zx())
                    );
            }

            if (mag(lower.zy()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+2, rowI+1, lower.zy())
                    );
            }

            if (mag(lower.zz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+2, rowI+2, lower.zz())
                    );
            }
        }
    }

    // Copy source vector into std::vector

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


    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A(m, m);
    //Eigen::SparseMatrix<scalar, Eigen::ColMajor, long int> A(m, m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());


    // Optionally export system to Matlab
    if (exportSystemMatlab_)
    {
        Info<< this->typeName << ": writing sparse matrix to "
            << "matlabSparseMatrix.txt and source to matlabSource.txt"
            << endl;

        // Write matrix
        Eigen::saveMarket(A, "matlabSparseMatrix.txt");

        // Write source
        OFstream sourceFile("matlabSource.txt");
        forAll(b, rowI)
        {
            sourceFile
                << b(rowI) << endl;
        }
    }


    // Solve the system

    Info<< this->typeName << ": iterative solution of sparse system" << endl;

    const scalar startT = matrix.mesh().thisDb().time().elapsedCpuTime();

    // Compressing matrix is meant to help performance
    A.makeCompressed();
    //Info<< this->typeName << ": makeCompressed done" << endl;

    // Bi-conjugate gradient stabilised iterative solver with ILUT
    // preconditioning
    Eigen::BiCGSTAB
    <
        Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>
    > solver(A);

    solver.setTolerance(tol_);
    solver.setMaxIterations(maxIter_);
    // solver.preconditioner().setDroptol(0.0001);
    // solver.preconditioner().setFillfactor(100);


    // Solve the system
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x = solver.solve(b);
    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> x2 = b;
    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> x2 = 0.0*b;
    // Initial guess of b (used by default in Eigen) works well where the NULL
    // vector sometimes causes divergence... hmmnn
    // Maybe it is to do with causing the initial residual to be large
    //Eigen::Matrix<scalar, Eigen::Dynamic, 1> x = solver.solveWithGuess(b, x2);

    Info<< "    solver time = "
        << (matrix.mesh().thisDb().time().elapsedCpuTime() - startT) << endl;

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


    // Other Eigen sparse solvers

    // Iterative solution
    // Eigen::GMRES
    //     //< Eigen::SparseMatrix<double>,
    //     //Eigen::DiagonalPreconditioner<double> >
    //     <
    //         Eigen::SparseMatrix<double>,
    //         Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> >
    //     >
    //     solver(A);

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
