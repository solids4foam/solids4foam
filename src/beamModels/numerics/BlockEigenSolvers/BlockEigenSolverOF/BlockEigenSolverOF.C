/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlockEigenSolverOF.H"

// #include "OFstream.H"
// #include "polyMesh.H"
// #include "addToRunTimeSelectionTable.H"
// #include "fvMesh.H"
// #include "beamModel.H"
#include "denseMatrixHelperFunctions.H"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolverOF, 0);

    // addToRunTimeSelectionTable
    // (
    //     blockVector6Solver, BlockEigenSolverOF, symMatrix
    // );

    // addToRunTimeSelectionTable
    // (
    //     blockVector6Solver, BlockEigenSolverOF, asymMatrix
    // );
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

void Foam::BlockEigenSolverOF::convertFoamMatrixToEigenMatrix
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei,
    Eigen::SparseMatrix<scalar>& A
)
{
    // if (BlockLduSolver::debug)
    // {
    //     Info<< this->typeName
    //         << ": copying matrix coefficients into Eigen format"
    //         << endl;
    // }

    const label nRows = 6*d.size();

    // Block CSR matrix storage

    //const labelList& rowPointers = mbMatrix.blockRowPointers();
    //const labelList& columnIndices = mbMatrix.blockColumnIndices();
    //const scalarField& coeffs = mbMatrix.blockCoeffs();

    //label blockSize = ::sqrt(coeffs.size()/columnIndices.size());

    // Create coefficient matrix: we must copy coeffs from CSR storage
    // Maybe it is possible to use CSR storage directly.
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(36*(d.size() + l.size() + u.size()));
    //-----------------------------------------------------------------------------
    //                  diagonal
    //-----------------------------------------------------------------------------
    label globalRowI = 0;
    forAll(d, cellI)
    {
        const scalarSquareMatrix& curD = d[cellI];

        for (label localRowI = 0; localRowI < 6; localRowI++)
        {
            for (label localColI = 0; localColI < 6; localColI++)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRowI + localRowI,
                        globalRowI + localColI,
                        curD(localRowI, localColI)
                    )
                );
            }
        }

        globalRowI += 6;
    }
//-----------------------------------------------------------------------------
//                  off-diagonal
//-----------------------------------------------------------------------------

    forAll(u, faceI)
    {
        const scalarSquareMatrix& curU = u[faceI];
        const scalarSquareMatrix& curL = l[faceI];

        const label globalRowI = 6*own[faceI];
        const label globalColI = 6*nei[faceI];

        //upper
        for (label localRowI = 0; localRowI < 6; localRowI++)
        {
            for (label localColI = 0; localColI < 6; localColI++)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRowI + localRowI,
                        globalColI + localColI,
                        curU(localRowI, localColI)
                    )
                );
            }
        }
        //lower
        for (label localRowI = 0; localRowI < 6; localRowI++)
        {
            for (label localColI = 0; localColI < 6; localColI++)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalColI + localRowI,
                        globalRowI + localColI,
                        curL(localRowI, localColI)
                    )
                );
            }
        }


    }
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------

        // Insert triplets into the matrix
        //label bnRows = rowPointers.size()-1;
        //nRows = blockSize*bnRows;
    A.resize(nRows, nRows);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Compressing matrix is meant to help performance
    A.makeCompressed();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenSolverOF::BlockEigenSolverOF
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei
)
:
    d_(d),
    l_(l),
    u_(u),
    own_(own),
    nei_(nei)
{}

// ************************************************************************* //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::BlockEigenSolverOF::solve
(
    Foam::Field<Foam::scalarRectangularMatrix>& foamX,
    const Foam::Field<Foam::scalarRectangularMatrix>& foamB
)
{
    // Allow to run in parallel, where each core solves its own independent
    // problem: this will not be correct if the beam is split across cores!
    // if (Pstream::parRun())
    // {
    //     FatalErrorIn
    //     (
    //         "bool Foam::BlockEigenSolverOF::solve"
    //         "(...)"
    //     )   << "Eigen direct linear solver may not be run in parallel"
    //         << abort(FatalError);
    // }

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A; // initialized in convertFoamMatrixToEigenMatrix funtion
    convertFoamMatrixToEigenMatrix(d_, l_, u_, own_, nei_, A);
    if (A.rows() == 0)
    {
        //    Pout << "No cells on this rank (" << Pstream::myProcNo() << "), skipping solver." << endl;
        return 0.0;
    }
    // Create Eigen source and solution vector from foam vectors
    //Eigen::Matrix<scalar, Eigen::Dynamic, 1> b;
    //Eigen::Matrix<scalar, Eigen::Dynamic, 1> x;

    const label nRows = A.rows();

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(nRows);
    label index = 0;
    forAll(foamB, cellI)
    {
        b(index++) = foamB[cellI](0,0);
        b(index++) = foamB[cellI](1,0);
        b(index++) = foamB[cellI](2,0);
        b(index++) = foamB[cellI](3,0);
        b(index++) = foamB[cellI](4,0);
        b(index++) = foamB[cellI](5,0);
    }

    // Copy solution vector into Eigen vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(nRows);
    index = 0;
    forAll(foamX, cellI)
    {
        x(index++) = foamX[cellI](0,0);
        x(index++) = foamX[cellI](1,0);
        x(index++) = foamX[cellI](2,0);
        x(index++) = foamX[cellI](3,0);
        x(index++) = foamX[cellI](4,0);
        x(index++) = foamX[cellI](5,0);
    }

    // Calculate initial residual
    const label nCells = d_.size();
    scalar initialResidual = 0;
    {

        Eigen::Matrix<scalar, Eigen::Dynamic, 1> Ax(nRows);
        Ax = A*x;

        Field<scalarRectangularMatrix> foamAx
        (
            nCells, scalarRectangularMatrix(6, 1, 0.0)
        );

        // Convert poroduct
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<6; j++)
            {
                foamAx[i](j,0) = Ax[k++];
            }
        }

        Field<scalarRectangularMatrix> blockR(foamB - foamAx);

        //scalarRectangularMatrix norm(6, 1, 1.0); //this->normFactor(U, blockB);
        //initialResidual = cmptDivide(gSum(cmptMag(blockR)), norm);
        // initialResidual = sqrt(gSum(magSqr(blockR)));
        // Note: use sum instead of gSum
        initialResidual = sqrt(sum(magSqr(blockR)));
    }


    typedef enum {SparseLU, BiCGSTAB, GMRES, DGMRES, MINRES} solvers;
    solvers sol = SparseLU;

    switch(sol)
    {
        case SparseLU:
        {
            Eigen::SparseLU
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::COLAMDOrdering<int>
            > solver(A);
            x = solver.solve(b);
            //solverPerf.nIterations()++;
            break;
        }
        case BiCGSTAB:
        {
            Eigen::BiCGSTAB
            <
                Eigen::SparseMatrix<scalar>,
                // Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        case GMRES:
        {
            Eigen::GMRES
            <
                Eigen::SparseMatrix<scalar>,
                // Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        case DGMRES:
        {
            Eigen::DGMRES
            <
                Eigen::SparseMatrix<scalar>,
                // Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        case MINRES:
        {
            Eigen::MINRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::Lower|Eigen::Upper,
                Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                // Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        default:
            FatalErrorIn
            (
                "Foam::BlockSolverPerformance<Foam::vector6>"
                "Foam::BlockEigenSolverOF::solve"
                "("
                "    Field<Foam::vector6>& U,"
                "    const Field<Foam::vector6>& blockB"
                ")"
            )   << "Undefined Eugen solver."
                << abort(FatalError);
    }

    // We copy the results from the std::vector into the geometric field
    index = 0;
    forAll(foamX, cellI)
    {
        foamX[cellI](0,0) = x(index++);
        foamX[cellI](1,0) = x(index++);
        foamX[cellI](2,0) = x(index++);
        foamX[cellI](3,0) = x(index++);
        foamX[cellI](4,0) = x(index++);
        foamX[cellI](5,0) = x(index++);
    }

    //

    // Calculate final residual
    // {
    //     Eigen::Matrix<scalar, Eigen::Dynamic, 1> p(nRows);
    //     p = A*x;

    //     Field<vector6> blockP(U.size());

    //     // Convert poroduct
    //     label k = 0;
    //     for (label i=0; i<nCells; i++)
    //     {
    //         for (label j=0; j<blockSize; j++)
    //         {
    //             blockP[i](j) = p[k++];
    //         }
    //     }

    //     Field<vector6> blockR(blockB - blockP);

    //     vector6 norm = vector6::one; //this->normFactor(U, blockB);

    //     solverPerf.finalResidual() =
    //         cmptDivide(gSum(cmptMag(blockR)), norm);
    // }


    return initialResidual;
}


