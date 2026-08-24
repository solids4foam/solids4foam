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
    Simple central-difference snGrad scheme without non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "simpleUncorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
simpleUncorrectedSnGrad<Type>::~simpleUncorrectedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
tmp<surfaceScalarField> simpleUncorrectedSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField& delta =
        mesh.lookupObject<surfaceVectorField>("delta");

    // const surfaceTensorField& gradDf =
    //     mesh.lookupObject<surfaceTensorField>("grad(D)f");

    // surfaceTensorField Ff = I + gradDf.T();

    // surfaceVectorField deltaN = (Ff & delta);

    return 1.0/mag(delta);

    // return this->mesh().deltaCoeffs();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
simpleUncorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    notImplemented
    (
        "simpleUncorrectedSnGrad<Type>::correction"
        "(const GeometricField<Type, fvPatchField, volMesh>&)"
    );
    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >(nullptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
