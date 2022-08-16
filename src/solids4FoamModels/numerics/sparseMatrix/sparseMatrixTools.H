/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    sparseMatrixTools

Description
    Helper functions for sparseMatrix

Author
    Philip Cardiff, UCD.

SourceFiles
    sparseMatrixTools.C

\*---------------------------------------------------------------------------*/

#ifndef sparseMatrixTools_H
#define sparseMatrixTools_H

#include "sparseMatrix.H"
#include "vectorField.H"
#include "polyMesh.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "SolverPerformance.H"
#else
    #include "BlockSolverPerformance.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace sparseMatrixTools
{

    //- Check if the case is 2-D
    bool checkTwoD(const polyMesh& mesh);

    //- Solve the linear system using Eigen's SparseLU direct solver
    void solveLinearSystemEigen
    (
        const sparseMatrix& matrix,
        const vectorField& source,
        vectorField& solution,
        const bool twoD,
        const bool exportToMatlab = false,
        const bool debug = false
    );

#ifdef USE_PETSC

    //- Solve the linear system using PETSc
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector> solveLinearSystemPETSc
#else
    BlockSolverPerformance<vector> solveLinearSystemPETSc
#endif
    (
        const sparseMatrix& matrix,
        const vectorField& source,
        vectorField& solution,
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
        sparseMatrix& matrix,
        vectorField& source,
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

} // End namespace sparseMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
