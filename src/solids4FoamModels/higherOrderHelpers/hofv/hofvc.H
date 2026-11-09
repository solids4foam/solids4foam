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
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef hofvc_H
#define hofvc_H

#include "CompactListList.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "pointFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class leastSquaresScheme;

namespace hofvc
{
    //- Update the gradient of a vector field for face quadrature points,
    //  using the given least squares scheme
    void fGrad
    (
        const leastSquaresScheme& ls,
        const volVectorField& vf,
        CompactListList<tensor>& gradQuad
    );

    //- Update the gradient of displacement for face quadrature points,
    //  using the displacement least squares scheme of the solid model
    void fGrad
    (
        const volVectorField& D,
        CompactListList<tensor>& gradDQuad
    );

    //- Surface integration of sigma using the quadrature points and weights
    //  of the given least squares scheme
    //  Integration of stress at quadrature points results with surface traction
    //  with dimensions dims
    tmp<surfaceVectorField> surfaceIntegrate
    (
        const leastSquaresScheme& ls,
        const CompactListList<symmTensor>& quadSigma,
        const fvMesh& mesh,
        const dimensionSet& dims = dimPressure
    );

    //- Surface integration of sigma using quadrature points and their weights
    //  Integration of stress at quadrature points results with surface traction
    tmp<surfaceVectorField> surfaceIntegrate
    (
        const CompactListList<symmTensor>& quadSigma,
        const fvMesh& mesh
    );

    //- Surface integration of the first Piola-Kirchhoff stress using the
    //  quadrature points and weights of the given least squares scheme
    //  The Cauchy stress and the displacement gradient at the quadrature points
    //  are given, and the integration results with the nominal traction, i.e.
    //  the force per unit undeformed area
    //  Unlike the Cauchy version, this has no dims argument, because the
    //  first Piola-Kirchhoff stress is used only by solids, so the traction is
    //  always in dimPressure
    tmp<surfaceVectorField> surfaceIntegrate
    (
        const leastSquaresScheme& ls,
        const CompactListList<symmTensor>& quadSigma,
        const CompactListList<tensor>& quadGradD,
        const fvMesh& mesh
    );

    //- Surface integration of the first Piola-Kirchhoff stress using quadrature
    //  points and their weights
    //  The Cauchy stress and the displacement gradient at the quadrature points
    //  are given, and the integration results with the nominal traction, i.e.
    //  the force per unit undeformed area
    tmp<surfaceVectorField> surfaceIntegrate
    (
        const CompactListList<symmTensor>& quadSigma,
        const CompactListList<tensor>& quadGradD,
        const fvMesh& mesh
    );

    tmp<volVectorField> d2dt2
    (
        const volVectorField& D
    );

} // End namespace hofvc

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
