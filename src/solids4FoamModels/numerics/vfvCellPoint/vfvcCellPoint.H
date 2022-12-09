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
    Explicit functions for vertex-centred finite volume calculations: Vertex
    Finite Volume Calculus (vfvc).

    Calculate the grad volTensorField given a dual mesh cell field, where the
    dual mesh cells corresponds to the primary mesh points.

    Calculate the divergence a dual mesh using a volSymmTensorField defined on
    the primary mesh.

SourceFiles
    vfvcCellPoint.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef vfvcCellPoint_H
#define vfvcCellPoint_H

#include "volFields.H"
#include "pointFields.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vfvc
{
    // Gradient field
    // Returns volTensorField on primary mesh
    tmp<volTensorField> grad
    (
        const pointVectorField& pointD, // primary mesh point displacement
        const fvMesh& mesh              // primary mesh
    );

    // Gradient field
    // Returns surfaceTensorField on dual mesh
    // Only calculates the gradient for internal faces
    tmp<surfaceTensorField> fGrad
    (
        const pointVectorField& pointD, // primary mesh point displacement
        const fvMesh& mesh,             // primary mesh
        const fvMesh& dualMesh,         // dual mesh
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalar zeta, // fraction of compact edge direction component used
        const bool debug = false
    );

    // Second time derivative pointVectorField pointD
    tmp<vectorField> d2dt2
    (
        ITstream& d2dt2Scheme,
        const pointVectorField& pointD, // displacement
        const pointVectorField& pointU, // velocity
        const pointVectorField& pointA, // acceleration
        const scalarField& pointRho,    // density
        const scalarField& pointVol,    // volumes
        const int debug // debug switch
    );

    // First time derivative pointVectorField pointP
    tmp<vectorField> ddt
    (
        ITstream& ddtScheme,
        ITstream& d2dt2Scheme,
        const pointVectorField& pointP // displacement
    );

    // <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
    // Returns laplacian of (gamma & pf) field on the dual mesh
    // tmp<pointVectorField> laplacian
    // (
    //     const tensor& gamma,
    //     const pointVectorField& pf,
    //     const fvMesh& dualMesh,
    //     const labelList& dualCellToPoint
    // );

} // End namespace vfvc

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "vfvcCellPoint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
