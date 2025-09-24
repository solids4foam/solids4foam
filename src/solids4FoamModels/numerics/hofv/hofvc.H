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
    Explicit functions for cell-centre high order finite volume calculations:
    High Order Finite Volume Calculus (hofvc).

SourceFiles
    hofvc.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef hofvc_H
#define hofvc_H

#ifdef USE_PETSC

#include <petscsnes.h>
#include <petscksp.h>

#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "sparseMatrix.H"
#include "scalarMatrices.H"
#include "symmTensor.H"
#include "RectangularMatrix.H"
#include "cellPointLeastSquaresVectors.H"
#include "SparseMatrixTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace hofvc
{
    //- Surface integration of sigma using quadrature points and their weights.
    //  Integration of stress at quadrature points results with surface traction
    tmp<surfaceVectorField> surfaceIntegrate
    (
        const List<List<symmTensor>>& quadSigma,
        const CompactListList<scalar>& quadW,
        const fvMesh& mesh
    );

    //- Surface integration of sigma using quadrature points and their weights.
    //  Sigma is mapped to reference configuration, surface traction refers
    //  to original configuration
    tmp<surfaceVectorField> surfaceIntegrate
    (
        const List<List<symmTensor>>& quadSigma,
        const List<List<tensor>>& gradDQuad,
        const CompactListList<scalar>& quadW,
        const fvMesh& mesh
    );

} // End namespace hofvc

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // #ifdef USE_PETSC

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
