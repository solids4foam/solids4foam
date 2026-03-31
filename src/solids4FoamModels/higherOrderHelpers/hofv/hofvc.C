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
#ifdef OPENFOAM_ORG
    #include "fvSchemes.H"
#endif

// #include <functional>
//#ifdef USE_PETSC
// #include "fvc.H"
// #include "multiplyCoeff.H"
// #include "multiplyCoeffExtended.H"
// #include "sparseMatrixTools.H"
// #include "cellPointLeastSquaresVectors.H"
// #include "surfaceFields.H"
// #include "processorPolyPatch.H"
// #include "emptyPolyPatch.H"
// #include "symmetryPolyPatch.H"
// #include "fixedDisplacementFvPatchVectorField.H"
// #include "solidTractionFvPatchVectorField.H"

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

    forAll (tf, faceI)
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

    forAll (tf.boundaryField(), patchI)
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
    return tmp<surfaceVectorField>();
#endif
};

tmp<volVectorField> hofvc::d2dt2
(
    const volVectorField& D
)
{
    // Hard-coded BDF2 scheme
    const fvMesh& mesh = D.mesh();

    tmp<volVectorField> tvf
    (
        new volVectorField
        (
            IOobject
            (
                 "d2dt2(" + D.name() + ")",
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zero",
                D.dimensions()/(dimTime*dimTime),
                vector::zero
            )
         )
    );

    //Check what scheme is prescribed, continue only in the case of backward
    word schemeName;
#ifdef OPENFOAM_ORG
    static_cast<const fvSchemes&>(mesh).d2dt2Scheme(D.name()) >> schemeName;
#else
    (mesh.d2dt2Schemes().found(D.name())
      ? mesh.d2dt2Schemes().lookup(D.name())
      : mesh.d2dt2Schemes().lookup("default")) >> schemeName;
#endif

    if (schemeName != "backward")
    {
        return tvf;
    }

    return tvf;
}

// ************************************************************************* //
