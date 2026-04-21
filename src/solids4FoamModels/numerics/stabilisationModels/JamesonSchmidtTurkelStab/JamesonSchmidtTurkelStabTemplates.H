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

\*---------------------------------------------------------------------------*/

#include "JamesonSchmidtTurkelStab.H"
#include "fvc.H"
#include "orthogonalSnGrad.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class GeoField, class SurfaceField>
void Foam::JamesonSchmidtTurkelStab::computeDiffStencil
(
    const GeoField& field,
    SurfaceField& faceFlux,
    scalar scaleFactor
) const
{
    // Take a reference to the mesh
    const fvMesh& mesh = field.mesh();

    // Create the orthogonal normal gradient operator
    using Type = typename GeoField::value_type;
    fv::orthogonalSnGrad<Type> op(mesh);

    // Inner pass: orthogonal diffusion operator
    const GeoField L1(fvc::div(op.snGrad(field)*mesh.magSf()));

    // Face-based h^2 scale via deltaCoeffs
    const surfaceScalarField h2f(1.0/sqr(mesh.deltaCoeffs()));

    // Outer pass: orthogonal diffusion of L1, scaled by h^2
    faceFlux = -scaleFactor*h2f*op.snGrad(L1);
}


// ************************************************************************* //
