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

Class
    sparseMatrixExtendedTools

Description
    Helper functions for sparseMatrix

Author
    Philip Cardiff, UCD.

SourceFiles
    sparseMatrixExtendedTools.C

\*---------------------------------------------------------------------------*/

#ifndef sparseMatrixExtendedTools_H
#define sparseMatrixExtendedTools_H

#include "sparseMatrixExtended.H"
#include "RectangularMatrix.H"
#include "vectorField.H"
#include "polyMesh.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "SolverPerformance.H"
#else
    #include "BlockSolverPerformance.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace sparseMatrixExtendedTools
{

    //- Check if the case is 2-D
    bool checkTwoD(const polyMesh& mesh);
    
    //- Solve the linear system using Eigen's SparseLU direct solver
    void solveLinearSystemEigen
    (
        const sparseMatrixExtended& matrix,
        const Field<RectangularMatrix<scalar>>& source,
        Field<RectangularMatrix<scalar>>& solution,
        const bool twoD,
        const bool exportToMatlab = false,
        const bool debug = false
    );

#ifdef USE_PETSC

    //- Solve the linear system using PETSc
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solveLinearSystemPETSc
#else
    BlockSolverPerformance<vector> solveLinearSystemPETSc
#endif
    (
        const sparseMatrixExtended& matrix,
        const Field<RectangularMatrix<scalar>>& source,
        Field<RectangularMatrix<scalar>>& solution,
        const bool twoD,
        fileName& optionsFile,
        const pointField& points,
        const boolList& ownedByThisProc,
        const labelList& localToGlobalPointMap,
        const labelList& stencilSizeOwned,
        const labelList& stencilSizeNotOwned,
        const bool debug = false
    );

#endif

    //- Set the number of non-zeros in per row in the given array
    void setNonZerosPerRow
    (
        label d_nnz[],
        label o_nnz[],
        label& d_nz,
        const int nRows,
        const int blockSize,
        const boolList& ownedByThisProc,
        const labelList& stencilSizeOwned,
        const labelList& stencilSizeNotOwned
    );

    //- Check the PETSc error
    //  Exit with a FatalError if it is greater than 0
    void checkErr(const int ierr);

    //- Enforce fixed DOF contributions on the linear system
    void enforceFixedDof
    (
                sparseMatrixExtended& matrix,
                Field<RectangularMatrix<scalar>>& source,
                const bool twoD,
        const boolList& fixedDofs,
        const symmTensorField& fixedDofDirections,
        const pointField& fixedDofValues,
        const scalar fixedDofScale
    );

    //- Add fixed DOF contributions to the source
    // void addFixedDofToSource
    // (
    //     vectorField& source,
    //     const boolList& fixedDofs,
    //     const symmTensorField& fixedDofDirections,
    //     const scalar fixedDofScale
    // );

    //- Add fixed DOF contributions to the matrix
    // void addFixedDofToMatrix
    // (
    //     sparseMatrix& matrix,
    //     const boolList& fixedDofs,
    //     const symmTensorField& fixedDofDirections,
    //     const scalar fixedDofScale
    // );

} // End namespace sparseMatrixExtendedTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
