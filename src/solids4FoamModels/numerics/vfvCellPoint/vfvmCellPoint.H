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
    Implicit functions for vertex-centred finite volume calculations: Vertex
    Finite Volume Matrix (vfvm).

    Populates the matrix for divergence of stress on a dual mesh, where the
    dual mesh cells corresponds to the primary mesh points.

    The d2dt2 operator is also defined.

    Note that these functions only set the matrix coefficients and they do not
    modify the right hand side.

SourceFiles
    vfvmCellPoint.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef vfvmCellPoint_H
#define vfvmCellPoint_H

#include "volFields.H"
#include "pointFields.H"
#include "sparseMatrix.H"
#include "scalarMatrices.H"
#include "RectangularMatrix.H"
#include "cellPointLeastSquaresVectors.H"
#include "SparseMatrixTemplate.H"
#include "meshDual.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vfvm
{
    // Add coefficients to the matrix for the divergence of stress
    // Note: this function does not calculate contributions to the right-hand
    // side
    void divSigma
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const Field<scalarSquareMatrix>& materialTangentField,
        const scalar zeta,
        const bool debug = false
    );

    // Add coefficients to the matrix for the divergence of stress
    // for a finite strain displacement formulation. Each coefficient in the
    // matrix is a 3x3 matrix.
    // Note: this function does not calculate contributions to the right-hand
    // side
    void divSigma
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const Field<scalarSquareMatrix>& materialTangentField,
        const Field<RectangularMatrix<scalar>>& geometricStiffnessField,
        const symmTensorField& sigma,
        const tensorField& dualGradDField,
        const boolList& fixedDofs,
        const symmTensorField& fixedDofDirections,
        const scalar fixedDofScale,
        const scalar zeta, // fraction of compact edge direction component used
        const bool debug = false
    );

    // Add coefficients to the matrix for the Laplacian of a scalar field with
    // a scalar diffusivity field
    void laplacian
    (
        sparseScalarMatrix& matrix,
        const Switch compactStencil,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalarField& diffusivity, // diffusivity in the primary cells
        const bool debug = false
    );

    // Create and return a sparse matrix with the (Newton) linearised
    // discretisation of the Laplacian term div(grad(pointP)) for a
    // pointVectorField pointP
    tmp<sparseMatrix> laplacian
    (
        const pointVectorField& pointP,
        const meshDual& dualMesh,
        const scalar zeta, // compact stencil coefficient
        const bool debug = false
    );

    // Create and return a sparse matrix with the (Newton) linearised
    // discretisation of the advection term div(rho*v*v)
    tmp<sparseMatrix> div
    (
        // const surfaceScalarField& dualRho, // dual surface density field
        // const scalar rho, // density
        const pointVectorField& pointV, // primary mesh velocity point field
        const meshDual& dualMesh, // dual mesh
        const bool debug = false
    );

    // Add coefficients to the matrix for the second time derivative
    // Note: this function does not calculate contributions to the right-hand
    // side
    void d2dt2
    (
        ITstream& d2dt2Scheme,
        const scalar& deltaT,           // time-step
        const word& pointDname,
        sparseMatrix& matrix,
        const scalarField& pointRhoI,
        const scalarField& pointVolI,
        const int debug  // debug switch
    );

    // Add coefficients to the matrix for the first time derivative
    // Note: this function does not calculate contributions to the right-hand
    // side
    tmp<sparseMatrix> ddt
    (
        ITstream& ddtScheme,
        const pointVectorField& pointP,
        const int debug = false
    );

} // End namespace vfvm

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
