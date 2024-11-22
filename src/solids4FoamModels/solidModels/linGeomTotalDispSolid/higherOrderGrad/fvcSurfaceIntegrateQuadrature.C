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

#include "fvcSurfaceIntegrateQuadrature.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> surfaceIntegrateQuadrature
(
    const List<vectorField>& traction,
    const List<scalarField>& quadratureWeights,
    const fvMesh& mesh
)
{
    tmp<volVectorField> tresult
    (
        new volVectorField
        (
            IOobject
            (
                "surfaceIntegrateQuadrature(traction)",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "0", dimPressure*dimArea/dimVolume, vector::zero
            ),
            "zeroGradient"
        )
    );
    volVectorField& result = tresult.ref();
    vectorField& resultI = result;

    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    forAll(faceOwner, faceI)
    {
        // Traction at the Gauss quadrature points on the face
        const vectorField& gaussQuadTracs = traction[faceI];

        // Gauss quadrature weights on the face
        const scalarField& gaussQuadW = quadratureWeights[faceI];

        // Add forces to the owner and neighbour cells
        forAll(gaussQuadTracs, pI)
        {
            // Force contribution for this quadrature point
            const vector contrib = gaussQuadTracs[pI]*gaussQuadW[pI];

            // Add to the owner cell
            const label ownCellID = faceOwner[faceI];
            resultI[ownCellID] += contrib;

            // Add to the neighbour cell, if there is one
            if (mesh.isInternalFace(faceI))
            {
                const label neiCellID = faceNeighbour[faceI];
                resultI[neiCellID] -= contrib;
            }
            // else
            // {
            //     // Do nothing as processor boundaries already take care
            //     // of their own cells
            // }
        }
    }

    // Divide by the volume
    resultI /= mesh.V();

    result.correctBoundaryConditions();

    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
