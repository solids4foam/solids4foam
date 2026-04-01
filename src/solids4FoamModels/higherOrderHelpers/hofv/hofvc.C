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

#include "hofvc.H"
#include "lookupSolidModel.H"
#include "fvcD2dt2.H"
#ifdef OPENFOAM_ORG
    #include "fvSchemes.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

void hofvc::fGrad
(
    const volVectorField& D,
    CompactListList<tensor>& gradDQuad
 )
{
#ifndef FOAMEXTEND
    const solidModel& solMod = lookupSolidModel(D.mesh());

    solMod.displacementMLS().fGrad(D, gradDQuad);
#else
    notImplemented("Not implemented for foam extend");
#endif
}


tmp<surfaceVectorField> hofvc::surfaceIntegrate
(
    const CompactListList<symmTensor>& quadSigma,
    const fvMesh& mesh
)
{
#ifndef FOAMEXTEND
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

    const vectorField normal(mesh.faceAreas()/mag(mesh.faceAreas()));
    const surfaceScalarField& magSf = mesh.magSf();
    // Reference to integration weights
    const solidModel& solMod = lookupSolidModel(mesh);
    const CompactListList<scalar>& quadW =
        solMod.displacementMLS().quadrature().faceQuadWeights();

    forAll(tf, faceI)
    {
        // Sigma at the quadrature points on the face
        const UList<symmTensor>& faceQuadStress = quadSigma[faceI];

        const vector& faceNormal = normal[faceI];

        forAll(faceQuadStress, pI)
        {
            // Add traction contribution of this quadrature point
            tf[faceI] += faceNormal & (faceQuadStress[pI] * quadW[faceI][pI]);
        }
        // We use physical weights so we need to divide with area to get
        // traction
        tf[faceI] *= (1.0/magSf[faceI]);
    }

    forAll(tf.boundaryField(), patchI)
    {
        vectorField& tfPatch = tf.boundaryFieldRef()[patchI];
        forAll(tfPatch, faceI)
        {
            const label globalFaceID =
                mesh.boundaryMesh()[patchI].start() + faceI;

            // Sigma at the Gauss quadrature points on the face
            const UList<symmTensor>& faceQuadStress = quadSigma[globalFaceID];

            const vector& faceNormal =  normal[globalFaceID];

            forAll(faceQuadStress, pI)
            {
                // Add traction contribution of this quadrature point
                tfPatch[faceI] +=
                    faceNormal & (faceQuadStress[pI] * quadW[globalFaceID][pI]);
            }

            tfPatch[faceI] *= (1.0/magSf.boundaryField()[patchI][faceI]);
        }
    }

    return tsf;
#else
    notImplemented("Not implemented for foam extend");
    return tmp<surfaceVectorField>();
#endif
};


tmp<volVectorField> hofvc::d2dt2
(
    const volVectorField& D
)
{
    // Default to second-order method for now
    // In the future, a consistent, higher-order scheme will be added
    return fvc::d2dt2(D);
}

// ************************************************************************* //
