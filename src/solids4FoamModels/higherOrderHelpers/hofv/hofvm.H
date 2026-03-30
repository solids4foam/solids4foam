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
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

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
#include "tensor.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "movingLeastSquares.H"
#include "foamPetscSnesHelper.H"

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

    //- High-order version of the d2dt2 term
    //  Currently I use fvm::d2dt2 here
    tmp<fvVectorMatrix> d2dt2(const volVectorField& D);

    //- Initialise the Jacobian matrix
    //  This procedure should set the size and non-zero structure of
    //  the matrix
    label initialiseJacobian
    (
        Mat& jac,
        const foamPetscSnesHelper& petscSnesHelper,
        const movingLeastSquares& mls,
        const volVectorField& D,
        const label blockSize,
        const bool createMat = true
    );

    // Add coefficients to the PETSc matrix for the Laplacian
    void laplacianIntoPETScMatrix
    (
        Mat jac,
        const foamPetscSnesHelper& petscSnesHelper,
        const movingLeastSquares& mls,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const label rowOffset = 0,
        const label colOffset = 0
     );

    // Add coefficients to the PETSc matrix for the Laplacian transpose
    void laplacianTransposeIntoPETScMatrix
    (
        Mat jac,
        const foamPetscSnesHelper& petscSnesHelper,
        const movingLeastSquares& mls,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const label rowOffset = 0,
        const label colOffset = 0
    );

    // Add coefficients to the matrix for the Laplacian trace
    void laplacianTraceIntoPETScMatrix
    (
        Mat jac,
        const foamPetscSnesHelper& petscSnesHelper,
        const movingLeastSquares& mls,
        const volVectorField& D,
        const volScalarField& diffusivity,
        const label rowOffset = 0,
        const label colOffset = 0
    );

} // End namespace hofvm

} // End namespace Foam

#endif

#endif

// ************************************************************************* //
