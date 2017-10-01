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
#include "BlockTrilinosSolver.H"
#include "BlockSolverPerformance.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "OFstream.H"
#include "polyMesh.H"

// Trilinos files
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_Redistributor.hpp>
#include <Isorropia_Partitioner.hpp>


//#include "ml_include.h"
//#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_TRIUTILS)
//#if 1
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "Ifpack.h"

// includes required by ML
//#include "ml_include.h"
#include "Epetra_LinearProblem.h"
//#include "ml_MultiLevelOperator.h"
//#include "ml_epetra_utils.h"
//using namespace ML_Epetra;
#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

// Function prototypes - Hisano
// void get_NumNz
// (
//     char Name[],
//     int MinMyElements,
//     int MaxMyElements,
//     int &NumMyNonZero,
//     int* &NumNz
// );
// void get_CRS_Matrix_info
// (
//     char NAME[],
//     int MinMyElements,
//     int MaxMyElements,
//     double *Values,
//     int *Indices
// );
// void read_b
// (
//     char NAME2[], int MinMyElements, int MaxMyElements, Epetra_Vector &b
// );
// void get_NumNz2
// (
//     char Name[],
//     int NumMyElements,
//     int *MyGLobalElements,
//     int &NumMyNonZero,
//     int* &NumNz,
//     int* &NumNz_total
// );
// void get_CRS_Matrix_info2
// (
//     char NAME[],
//     int NumMyElements,
//     int* MyGlobalElements,
//     int* NumNz_total,
//     double *Values,
//     int *Indices
// );
// void read_b2
// (
//     char NAME2[], int NumMyElements, int* MyGLobalElements, Epetra_Vector &b
// );
// void get_CRS_Matrix_info_for_OpenFORM
// (
//     char NAME[],
//     int MinMyElements,
//     int MaxMyElements,
//     double *&Values,
//     int *&Indices,
//     int *&Num_NZ
// );
// void get_CRS_Matrix_info_for_OpenFOAM2
// (
//     char NAME[],
//     char NAME2[],
//     char NAME3[],
//     int MinMyElements,
//     int MaxMyElements,
//     double *&Values,
//     int *&Indices,
//     int *&Num_NZ
// );

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockTrilinosSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockTrilinosSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockTrilinosSolver, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockTrilinosSolver::BlockTrilinosSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector>(fieldName, matrix, dict),
    tol_(dict.lookupOrDefault<scalar>("tolerance", 1000)),
    nDir_(readInt(dict.lookup("nDirections")))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockTrilinosSolver::solve
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
            "Foam::BlockTrilinosSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "This linear solver may not be run in parallel"
            << abort(FatalError);
    }

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
    label nRows = 0;
    if (twoD)
    {
        // 2-D
        nRows = 2*d.size();
    }
    else
    {
        // 3-D
        nRows = 3*d.size();
    }


    // Convert system to compressed row format: not directly needed
    // Info<< "Converting matrix to compressed row format" << endl;

    // int n_rows = d.size();
    // int nnz = l.size() + u.size() + d.size();

    // std::vector<scalar> vals(nnz);
    // std::vector<int> c_idx(nnz);
    // std::vector<int> r_idx(n_rows + 2, 1);

    // //ldu2csr(matrix_, c_idx, r_idx, vals) ;
    // //matrix.addBoundaryDiag(diag, 0); // not needed for block
    // ldu2csr(d, u, l, upperAddr, lowerAddr, c_idx, r_idx, vals) ;


    // Convert matrix format into Epetra format

    // Number of non-zeros in each row
    // Note: the diagonal coefficient contributes 'nDir' non-zeros
    labelList nNonZeros(nRow, nDir);

    // Count off-diagonal non-zeros in each row
    forAll(u, bondI)
    {
        // Note: lowerAddr gives the owner of a bond and
        // upperAddr gives the neighbour
        label lI = lowerAddr[bondI];
        label uI = upperAddr[bondI];

        // Each coefficient is a 2x2 or 3x3 tensor, so it contributes nDir
        // non-zeros to each row
        nNonZeros[lI] += nDir;
        nNonZeros[uI] += nDir;
    }

    // Non-zero coefficients for each row
    scalarListList rowValues(nRows);

    // List to count how many entried we have added to each row
    labelList rowColI(nRows, 0);

    // Column IDs for non-zero cofficients in each row
    labelListList rowColumnIDs(nRows);

    // Set the size of rowValues and rowColumnIDs for each row
    forAll(rowValues, rowI)
    {
        label n = nNonZeros[rowI];
        rowValues.setSize(n, 0.0);
        rowColumnIDs.setSize(n, -1);
    }

    // Fill in the rowValues and rowColumnIDs for each row

    // Insert the diagonal

    label rowID = -1;
    forAll(d, dI)
    {
        if (twoD)
        {
            rowID = 2*dI;
        }
        else
        {
            rowID = 3*dI;
        }

        // Insert xx, xy, and xz coefficients

        label colID = rowColI[rowID];
        rowValues[rowID][colID] = d[rowID].xx();
        rowColumn[rowID][colID] = rowI;

        rowValues[rowID][colID + 1] = d[rowID].xy();
        rowColumn[rowID][colID + 1] = rowI;

        rowValues[rowID + 1][colID] = d[rowID].yx();
        rowColumn[rowID + 1][colID] = rowI;

        rowValues[rowID + 1][colID + 1] = d[rowID].yy();
        rowColumn[rowID + 1][colID + 1] = rowI;

        if (!twoD)
        {
            rowValues[rowID][colID + 2] = d[rowID].xz();
            rowColumn[rowID][colID + 2] = rowI;

            rowValues[rowID + 1][colID + 2] = d[rowID].yz();
            rowColumn[rowID + 1][colID + 2] = rowI;

            rowValues[rowID + 2][colID] = d[rowID].zx();
            rowColumn[rowID + 2][colID] = rowI;

            rowValues[rowID + 2][colID + 1] = d[rowID].zy();
            rowColumn[rowID + 2][colID + 1] = rowI;

            rowValues[rowID + 2][colID + 2] = d[rowID].zz();
            rowColumn[rowID + 2][colID + 2] = rowI;
        }

        // Increment the number of columns set in the row
        rowColI[rowID] += nDir;
    }


    // Insert the upper and lower coefficients
    todo(); thisFar();
    forAll(u, bondI)
    {
        label lI = lowerAddr[bondI];
        label uI = upperAddr[bondI];

        rowValues[lI][rowColI[lI]] = u[bondI];
        rowColumn[lI][rowColI[lI]] = uI;
        rowColI[lI]++;

        rowValues[uI][rowColI[uI]] = l[bondI];
        rowColumn[uI][rowColI[uI]] = lI;
        rowColI[uI]++;
    }

    // Debug: Check everything was set
    forAll(rowColumnID, rowI)
    {
        if (min(rowColumnID[rowI]))
        {
            FatalErrorIn("solve()")
                << "Problem converting matrix" << abort(FatalError);
        }
    }



    // Create Trilinos Epetra system

    // Set Epetra_Map
    Epetra_Map eMap(-1, d.size(), 0, Comm);

    // Set Epetra_CrsMatrix
    Epetra_CrsMatrix A(Copy, eMap, NumNz); // Number of Non-Zeros

    // Insert Values into Epetra_CRS matrix
    label current_pos = 0;
    forAll(d, rowI)
    {
        A.InsertGlobalValues
        (
            rowI,                       // row number (can be in be parallel)
            nNonZeros[rowI],            // number of non-zeros in the row
            rowValues[rowI],            // value
            rowColumnIDs[rowI]          // column number
            // &Values[current_pos],    // value
            // &Indices[current_pos]    // column number
        );
        current_pos += NumNz[i];
    }
    A.FillComplete();

    // Create solution vector
    Epetra_Vector x(Map);

    // Copy initial guess into solution vector
    todo();

    // Create source vector
    Epetra_Vector b(Map);
    todo();

    // Create the linear system
    // this is the linear problem to be solved: set the linear operator,
    // the solution and the right-hand side
    Epetra_LinearProblem A_Problem(&A, &x, &b);

    // Create the linear solver
    AztecOO A_Solver(A_Problem);

    A_Solver.SetAztecOption(AZ_scaling, 4);
    A_Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    A_Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
    A_Solver.SetAztecOption(AZ_graph_fill, 0);
    A_Solver.SetAztecOption(AZ_overlap, 1);

    A_Solver.SetAztecOption(AZ_kspace, nDir_);
    A_Solver.SetAztecOption(AZ_solver, AZ_gmres);


    // Solve the system

    Info<< this->typeName << ": iterative solution of sparse system" << endl;

    const scalar startT = matrix.mesh().thisDb().time().elapsedCpuTime();

    A_Solver.Iterate(Niters, tol_);

    Info<< "    solver time = "
        << (matrix.mesh().thisDb().time().elapsedCpuTime() - startT) << endl;


    // Copy solution vector back into foam format

    Info<< this->typeName << ": copying results into foam format" << nl << endl;

    todo();

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
