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

#include "deltaVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceVectorField> Foam::deltaVectors
(
    const fvMesh& mesh
)
{
#ifdef FOAMEXTEND
    tmp<surfaceVectorField> tdelta
    (
        new surfaceVectorField
        (
            IOobject
            (
               "deltaCoeffs",
               mesh.pointsInstance(),
               mesh
            ),
            mesh,
            dimLength
        )
    );
    surfaceVectorField& delta = tdelta();

    const volVectorField& C = mesh.C();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        delta[faceI] = C[neighbour[faceI]] - C[owner[faceI]];
    }

    auto& deltabf =  delta.boundaryField();

    forAll(deltabf, patchI)
    {
        const labelList faceCells(mesh.boundary()[patchI].faceCells());

        if (deltabf[patchI].coupled())
        {
            deltabf[patchI] =
                C.boundaryField()[patchI].patchNeighbourField()
              - C.boundaryField()[patchI].patchInternalField();
        }
        else
        {
            deltabf[patchI] =
                C.boundaryField()[patchI]
              - C.boundaryField()[patchI].patchInternalField();
        }
    }

    return tdelta;
#else
    return mesh.delta();
#endif
}

// ************************************************************************* //
