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

#include "diffStencilLaplacianStab.H"
#include "fvc.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class GeoField, class GradField, class SurfaceField>
void Foam::diffStencilLaplacianStab::computeDiffStencil
(
    const GeoField& field,
    const GradField& gradField,
    SurfaceField& faceFlux,
    scalar scaleFactor
) const
{
    // Take a reference to the mesh
    const fvMesh& mesh = field.mesh();

    // Unit normals
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Update the stabilisation
    faceFlux =
        scaleFactor*(fvc::snGrad(field) - (n & fvc::interpolate(gradField)));
}


// ************************************************************************* //
