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

#include "generalisedEvenOrderLaplacianStab.H"
#include "fvc.H"
#include "orthogonalSnGrad.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class GeoField, class SurfaceField>
void Foam::generalisedEvenOrderLaplacianStab::computeDiffStencil
(
    const GeoField& field,
    SurfaceField& faceFlux,
    scalar scaleFactor,
    const label laplacianPower   // m ≥ 1
) const
{
    // Take a reference to the mesh
    const fvMesh& mesh = field.mesh();

    // Create the orthogonal normal gradient operator
    using Type = typename GeoField::value_type;
    fv::orthogonalSnGrad<Type> op(mesh);

    // Working field
    GeoField L("L", field);

    // Apply orthogonal Laplacian m times
    for (label i = 0; i < laplacianPower; ++i)
    {
        L = fvc::div(op.snGrad(L)*mesh.magSf());
    }

    // h^(2(m-1)) scaling
    const surfaceScalarField h2fPow
    (
        pow(1.0/sqr(mesh.deltaCoeffs()), laplacianPower - 1)
    );

    // Ensure damping sign
    const scalar sign = (laplacianPower % 2 == 0) ? 1.0 : -1.0;

    // Calculate the flux
    faceFlux = sign*scaleFactor*h2fPow*op.snGrad(L)*mesh.magSf();
}


// ************************************************************************* //
