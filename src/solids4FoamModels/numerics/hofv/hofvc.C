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

#ifdef USE_PETSC

#include <functional>
#include "hofvc.H"
#include "fvc.H"
#include "multiplyCoeff.H"
#include "multiplyCoeffExtended.H"
#include "sparseMatrixTools.H"
#include "cellPointLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"
#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

tmp<surfaceVectorField> hofvc::surfaceIntegrate
(
    const List<List<symmTensor>>& quadSigma,
    const List<List<scalar>>& quadW,
    const fvMesh& mesh
)
{
    tmp<surfaceVectorField> tsf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "surfaceIntegrate(sigma)",
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("0", dimPressure, vector::zero)
        )
    );

    surfaceVectorField& tf = tsf.ref();

    forAll (tf, faceI)
    {
        // Sigma at the quadrature points on the face
        const List<symmTensor>& faceQuadStress = quadSigma[faceI];

        // Weights for the quadrature points on the face
        const List<scalar>& faceQuadW = quadW[faceI];

        const vector& faceNormal =
            mesh.faceAreas()[faceI]/mag(mesh.faceAreas()[faceI]);

        forAll(faceQuadStress, pI)
        {
            // Add traction contribution of this quadrature point
            tf[faceI] += faceNormal & (faceQuadStress[pI] * faceQuadW[pI]);
        }
    }

    forAll (tf.boundaryField(), patchI)
    {
        vectorField& tfPatch = tf.boundaryFieldRef()[patchI];
        forAll(tfPatch, faceI)
        {
            const label globalFaceID =
		mesh.boundaryMesh()[patchI].start() + faceI;

            // Sigma at the Gauss quadrature points on the face
            const List<symmTensor>& faceQuadStress = quadSigma[globalFaceID];

            // Gauss quadrature weights on the face
            const List<scalar>& faceQuadW = quadW[globalFaceID];

            const vector& faceNormal =
                mesh.faceAreas()[globalFaceID]
	      / mag(mesh.faceAreas()[globalFaceID]);

            forAll(faceQuadStress, pI)
            {
                // Add traction contribution of this quadrature point
                tfPatch[faceI] +=
		    faceNormal & (faceQuadStress[pI] * faceQuadW[pI]);
            }
        }
    }

    return tsf;
};


// ************************************************************************* //

#endif // #ifdef USE_PETSC
