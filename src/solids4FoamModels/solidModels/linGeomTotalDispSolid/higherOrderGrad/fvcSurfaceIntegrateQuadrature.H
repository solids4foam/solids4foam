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

InNamespace
    Foam::fvc

Description
    Integrate a surface vector field given in terms Gaussian quadrature values
    over the cells and return a volVectorField.

    Example usage:

        volVectorField divSigma = fvc::surfaceIntegrateQuadrature(traction)

    where traction is a Field<vector> giving the traction at Gaussian
    quadrature points on the faces.

    All non-triangular polygonal faces are implicitly decomposed into triangles
    using the face centroid, and quadrature points are defined on each triangle
    sub-face.

SourceFiles
    fvcSurfaceIntegrateQuadrature.C

\*---------------------------------------------------------------------------*/


#ifndef fvcSurfaceIntegrateQuadrature_H
#define fvcSurfaceIntegrateQuadrature_H

#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    tmp<volVectorField> surfaceIntegrateQuadrature
    (
        const List<vectorField>& traction,
        const List<scalarField>& quadratureWeights,
        const fvMesh& mesh
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
