/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
        const Field<symmTensor4thOrder>& materialTangentField,
        const boolList& fixedDofs,
        const symmTensorField& fixedDofDirections,
        const scalar fixedDofScale,
        const scalar zeta, // fraction of compact edge direction component used
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

} // End namespace vfvc

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "vfvmCellPoint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
