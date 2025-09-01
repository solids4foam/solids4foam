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

Description
    Implicit functions for cell-centre high order finite volume calculations:
    High Order Finite Volume Matrix (hofvm).

SourceFiles
    hofvm.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef hofvm_H
#define hofvm_H

#ifdef USE_PETSC

#include <petscsnes.h>
#include <petscksp.h>

#include "volFields.H"
#include "pointFields.H"
#include "sparseMatrix.H"
#include "scalarMatrices.H"
#include "RectangularMatrix.H"
#include "cellPointLeastSquaresVectors.H"
#include "SparseMatrixTemplate.H"
#include "LRE.H"
#include "foamPetscSnesHelper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace hofvm
{
    //- Laplacian coefficient
    tensor laplacianCoeff
    (
        const scalar& gammaMagSf,
        const scalar& quadWeight,
        const vector& gradInterpCoeff,
        const vector& faceNormal
    );

    //- Laplacian transpose coefficient
    tensor laplacianTransposeCoeff
    (
        const scalar& gammaMagSf,
        const scalar& quadWeight,
        const vector& gradInterpCoeff,
        const vector& faceNormal
     );

    //- Laplacian trace coefficient
    tensor laplacianTraceCoeff
    (
        const scalar& gammaMagSf,
        const scalar& quadWeight,
        const vector& gradInterpCoeff,
        const vector& faceNormal
     );

    // Add coefficients to the PETSc matrix for the Laplacian
    void laplacianIntoPETScMatrix
    (
        Mat matrix,
        const foamPetscSnesHelper& petscSnesHelper,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre,
        const label nScalarEqns,
        const label rowOffset = 0,
        const label colOffset = 0
     );

    // Add coefficients to the sparse matrix for the Laplacian
    void laplacianIntoSparseMatrix
    (
        sparseMatrix& matrix,
        vectorField& source,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre
    );

    // Add coefficients to the PETSc matrix for the Laplacian transpose
    void laplacianTransposeIntoPETScMatrix
    (
        Mat matrix,
        const foamPetscSnesHelper& petscSnesHelper,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre,
        const label nScalarEqns,
        const label rowOffset = 0,
        const label colOffset = 0
    );

    // Add coefficients to the sparse matrix for the Laplacian transpose
    void laplacianTransposeIntoSparseMatrix
    (
        sparseMatrix& matrix,
        vectorField& source,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre
    );

    // Add coefficients to the matrix for the Laplacian trace
    void laplacianTraceIntoPETScMatrix
    (
        Mat matrix,
        const foamPetscSnesHelper& petscSnesHelper,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre,
        const label nScalarEqns,
        const label rowOffset = 0,
        const label colOffset = 0
    );

    // Add coefficients to the matrix for the Laplacian trace
    void laplacianTraceIntoSparseMatrix
    (
        sparseMatrix& matrix,
        vectorField& source,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre
    );

    // Function that is doing all work, matrix coefficients are calculated
    // using static functions for coefficient calculations
    void hofvmLaplacianSparseMatrix
    (
        sparseMatrix& matrix,
        vectorField& source,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre,
        tensor (*calcCoeff)
        (
             const scalar& gammaMagSf,
             const scalar& quadWeight,
             const vector& gradInterpCoeff,
             const vector& faceNormal
        )
    );

    // Function that is doing all work, matrix coefficients are calculated
    // using static functions for coefficient calculations
    label hofvmLaplacianPETSc
    (
        Mat matrix,
        const foamPetscSnesHelper& petscSnesHelper,
        const fvMesh& mesh,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const LRE& lre,
        tensor (*calcCoeff)
        (
             const scalar& gammaMagSf,
             const scalar& quadWeight,
             const vector& gradInterpCoeff,
             const vector& faceNormal
        ),
        const label rowOffset,
        const label colOffset,
        const label nScalarEqns
    );

    // At traction boundary we include traction directly to source vector.
    // Traction boundaries are skipped when calculating laplace, laplaceTrace
    // and laplaceTranspose
    void addTractionBoundaries
    (
        vectorField& source,
        const fvMesh& mesh,
        const volVectorField& D,
        const LRE& lre
     );

} // End namespace hofvm

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // #ifdef USE_PETSC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
