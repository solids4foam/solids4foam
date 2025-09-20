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

#ifdef USE_PETSC

#ifndef vfvmCellPoint_H
#define vfvmCellPoint_H

#include "volFields.H"
#include "pointFields.H"
#include "scalarMatrices.H"
#include "cellPointLeastSquaresVectors.H"
#include "meshDual.H"
#include "matNN.H"
#include <petscksp.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vfvm
{
    //- Insert vfvm::divSigma coefficients into the given PETSc matrix
    void divSigma
    (
        Mat jac,
        const pointVectorField& pointD,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const label nScalarEqns,
        const labelList& localToGlobalPointMap,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const List<mat66>& materialTangent,
        const scalar zetaImplicit,
        const bool flipSign = false
    );

    // Add coefficients to the matrix for the divergence of stress
    // for a finite strain displacement formulation. Each coefficient in the
    // matrix is a 3x3 matrix.
    // Note: this function does not calculate contributions to the right-hand
    // side
    void divSigma
    (
        Mat jac,
        const pointVectorField& pointD,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const label nScalarEqns,
        const labelList& localToGlobalPointMap,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const List<mat66>& materialTangent,
        const List<mat39>& geometricStiffnessField,
        const symmTensorField& sigma,
        const tensorField& dualGradDField,
        const scalar zetaImplicit,
        const bool flipSign = false
    );

    //- Insert vfvm::d2dt2 coefficients into the given PETSc matrix
    void d2dt2
    (
        Mat jac,
        const pointVectorField& pointD,
        const scalarField& pointRhoI,
        const scalarField& pointVolI,
        ITstream& d2dt2Scheme,
        const label nScalarEqns,
        const labelList& localToGlobalPointMap,
        const bool flipSign = false
    );

    // Insert laplacian coefficients into the given PETSc matrix
    void laplacian
    (
        Mat jac,
        const Switch compactStencil,
        const scalar zeta, // compact stencil coefficient
        const meshDual& dualMesh,
        const label nScalarEqns,
        const labelList& localToGlobalPointMap,
        const scalarField& diffusivity, // diffusivity in the primary cells
        const bool flipSign = false
    );

    // // Create and return a sparse matrix with the (Newton) linearised
    // // discretisation of the advection term div(rho*v*v)
    // tmp<sparseMatrix> div
    // (
    //     // const surfaceScalarField& dualRho, // dual surface density field
    //     // const scalar rho, // density
    //     const pointVectorField& pointV, // primary mesh velocity point field
    //     const meshDual& dualMesh, // dual mesh
    //     const bool debug = false
    // );

} // End namespace vfvm

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // ifdef USE_PETSC

// ************************************************************************* //
