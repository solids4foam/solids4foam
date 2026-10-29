/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pseudoVector.H"
#include "emptyFvPatchFields.H"
#include "spinTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<volVectorField> pseudoVector(const volTensorField& rotationMatrix)
{
    tmp<volVectorField> tresult
    (
        new volVectorField
        (
            IOobject
            (
                "pseudoVector("+rotationMatrix.name()+")",
                rotationMatrix.time().timeName(),
                rotationMatrix.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotationMatrix.mesh(),
            dimensionedVector("0", dimless, vector::zero)
        )
    );

    vectorField& resultI = tresult.ref().primitiveFieldRef();

    const tensorField& rotationMatrixI = rotationMatrix.internalField();

    forAll(resultI, cellI)
    {
        resultI[cellI] = pseudoVector(rotationMatrixI[cellI]);
    }

    forAll (rotationMatrix.boundaryField(), patchI)
    {
        vectorField& pResultI = tresult.ref().boundaryFieldRef()[patchI];

        const tensorField& pRotationMatrixI =
            rotationMatrix.boundaryField()[patchI];

        forAll(pResultI, faceI)
        {
            pResultI[faceI] = pseudoVector(pRotationMatrixI[faceI]);
        }
    }

    return tresult;
}

tmp<vectorField> pseudoVector(const tensorField& rotationMatrix)
{
    tmp<vectorField> tresult
    (
        new vectorField(rotationMatrix.size(), vector::zero)
    );

    vectorField& result = tresult.ref();

    forAll(result, cellI)
    {
        result[cellI] = pseudoVector(rotationMatrix[cellI]);
    }

    return tresult;
}

vector pseudoVector(const tensor& R)
{
    vector result = vector::zero;

    scalar q0 = 0;
    vector q = vector::zero;

    scalarField D(4, 0);
    D[0] = R.xx();
    D[1] = R.yy();
    D[2] = R.zz();
    D[3] = tr(R);

    const scalar a = max(D);

    {
        label i = -1;
        for (label j = 0; j < 4; j++)
        {
            if ( mag(a - D[j]) < SMALL )
            {
                i = j;
                break;
            }
        }

        switch(i)
        {
            case 0:
                q.x() = sqrt(a/2 + (1 - tr(R))/4);
                q0 = (R.zy() - R.yz())/(4*q.x());
                q.y() = (R.yx() + R.xy())/(4*q.x());
                q.z() = (R.zx() + R.xz())/(4*q.x());
                break;
            case 1:
                q.y() = sqrt(a/2 + (1 - tr(R))/4);
                q0 = (R.xz() - R.zx())/(4*q.y());
                q.z() = (R.zy() + R.yz())/(4*q.y());
                q.x() = (R.xy() + R.yx())/(4*q.y());
                break;
            case 2:
                q.z() = sqrt(a/2 + (1 - tr(R))/4);
                q0 = (R.yz() - R.zy())/(4*q.z());
                q.x() = (R.xz() + R.zx())/(4*q.z());
                q.y() = (R.yz() + R.zy())/(4*q.z());
                break;
            case 3:
                q0 = 0.5*::sqrt(1 + a);
                q.x() = (R.zy() - R.yz())/(4*q0);
                q.y() = (R.xz() - R.zx())/(4*q0);
                q.z() = (R.yx() - R.xy())/(4*q0);
                break;
            default:
                FatalErrorIn
                (
                    "vector pseudoVector(const tensor& R)"
                )
                  << "Error in calculation of quaternions"
                  << "for rotation matrix " << R
                  << abort(FatalError);
        }
    }

    const scalar theta = atan(mag(q)/q0)*2;

    if (mag(theta) > SMALL)
    {
        result = theta*(q/mag(q));
    }

    return result;
}


tmp<tensorField> rotationMatrix(const vectorField& rotationAngle)
{
    tmp<tensorField> tRotMat
    (
        new tensorField(rotationAngle.size(), tensor::I)
    );

    const scalarField magAngle(mag(rotationAngle) + SMALL);
    const tensorField angleHat(spinTensor(rotationAngle));

    tRotMat.ref() =
        I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
            *(angleHat & angleHat);

    return tRotMat;
}


tmp<volTensorField> rotationMatrix(const volVectorField& rotationAngle)
{
    tmp<volTensorField> tRotMat
    (
        new volTensorField
        (
            IOobject
            (
                "RotMat("+rotationAngle.name()+")",
                rotationAngle.time().timeName(),
                rotationAngle.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotationAngle.mesh(),
            dimensionedTensor("0", dimless, tensor::I)
        )
    );

   const volScalarField magAngle(mag(rotationAngle) + SMALL);
   const volTensorField angleHat(spinTensor(rotationAngle));

    tRotMat.ref() =
        I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
            *(angleHat & angleHat);

    return tRotMat;
}


tmp<surfaceTensorField> rotationMatrix(const surfaceVectorField& rotationAngle)
{
    tmp<surfaceTensorField> tRotMat
    (
        new surfaceTensorField
        (
            IOobject
            (
                "RotMat("+rotationAngle.name()+")",
                rotationAngle.time().timeName(),
                rotationAngle.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotationAngle.mesh(),
            dimensionedTensor("0", dimless, tensor::I)
        )
    );

    const surfaceScalarField magAngle (mag(rotationAngle) + SMALL);
    const surfaceTensorField angleHat (spinTensor(rotationAngle));

    tRotMat.ref() =
        I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
           *(angleHat & angleHat);

    return tRotMat;
}

tensor rotationMatrix(const vector& rotationAngle)
{
    tensor R = tensor::zero;

    const scalar magAngle = mag(rotationAngle) + SMALL;
    const tensor angleHat = spinTensor(rotationAngle);

    // Rodrigues formula
    R =
        I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
           *(angleHat & angleHat);

    return R;
}

/*
    NOTE: currently unused - kept for reference.

    Objective cell-to-face interpolation: the face rotation is the owner
    rotation followed by half of the relative rotation between the owner
    and neighbour cells. Boundary faces take the boundary values of R.

    objectiveInterpolation (beamProperties) updated the rotations after each
    Newton correction as Lambda = rotationMatrix(DTheta) & Lambda, then
    recalculated K = refLambdaf.T() & meanLineCurvature(Lambda) and
    Lambdaf = interpolateRotationMatrix(Lambda). The option was removed from
    the solver and the tutorials in the OpenFOAM.com port (1e2a795); the
    original code marked both calls "this does not work". The last full
    implementation is on the openfoam-v2306 branch
    (coupledTotalLagNewtonRaphsonBeamEvolve.C) and in commit e54de7d.
*/
tmp<surfaceTensorField> interpolateRotationMatrix
(
    const volTensorField& R
)
{
    tmp<surfaceTensorField> tRf
    (
        new surfaceTensorField
        (
            IOobject
            (
                "Rf",
                R.time().timeName(),
                R.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R.mesh(),
            dimensionedTensor("0", dimless, tensor::I)
        )
    );

    surfaceTensorField& Rf = tRf.ref();

    forAll(Rf.boundaryField(), patchI)
    {
    	// For boundary faces, why are we substituting the cell centre rotation values?
        Rf.boundaryFieldRef()[patchI] =
            R.boundaryField()[patchI];
    }

    const tensorField& RI = R.internalField();

    tensorField& RfI = Rf.primitiveFieldRef();

    const fvMesh& mesh = R.mesh();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    forAll(RfI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        // scalar s = 0.5/deltaCoeffsI[faceI];
        RfI[faceI] = RI[own[faceI]] & rotationMatrix(0.5*DRotAngle);

    }

    return tRf;
}


/*
    NOTE: currently unused - kept for reference.

    Face curvature from the relative rotation vector between the owner and
    neighbour cells divided by the cell-centre distance (local frame of the
    owner cell). Used together with interpolateRotationMatrix(R) above.

    objectiveInterpolation (beamProperties) updated the rotations after each
    Newton correction as Lambda = rotationMatrix(DTheta) & Lambda, then
    recalculated K = refLambdaf.T() & meanLineCurvature(Lambda) and
    Lambdaf = interpolateRotationMatrix(Lambda). The option was removed from
    the solver and the tutorials in the OpenFOAM.com port (1e2a795); the
    original code marked both calls "this does not work". The last full
    implementation is on the openfoam-v2306 branch
    (coupledTotalLagNewtonRaphsonBeamEvolve.C) and in commit e54de7d.
*/
tmp<surfaceVectorField> meanLineCurvature
(
    const volTensorField& R
)
{
    tmp<surfaceVectorField> tKf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Kf",
                R.time().timeName(),
                R.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R.mesh(),
            dimensionedVector("0", dimless/dimLength, vector::zero)
        )
    );

    surfaceVectorField& Kf = tKf.ref();

    const tensorField& RI = R.internalField();

    vectorField& KfI = Kf.primitiveFieldRef();

    const fvMesh& mesh = R.mesh();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs();

    forAll(KfI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        KfI[faceI] = DRotAngle*deltaCoeffsI[faceI];
    }

    forAll(Kf.boundaryField(), patchI)
    {
        vectorField& pKf = Kf.boundaryFieldRef()[patchI];
        const tensorField& pR = R.boundaryField()[patchI];

        const labelList faceCells = mesh.boundary()[patchI].faceCells();
        const scalarField& pDeltaCoeffs =
            mesh.deltaCoeffs().boundaryField()[patchI];

        forAll(pKf, faceI)
        {
            tensor DR = (RI[faceCells[faceI]].T() & pR[faceI]);
            vector DRotAngle = pseudoVector(DR);

            tensor RN = (RI[faceCells[faceI]] & rotationMatrix(2*DRotAngle));

            DR = (RI[faceCells[faceI]].T() & RN);
            DRotAngle = pseudoVector(DR);

            pKf[faceI] = DRotAngle*pDeltaCoeffs[faceI]/2;
        }
    }

    return tKf;
}

/*
    NOTE: currently unused - kept for reference.

    Combines interpolateRotationMatrix(R) and rotationAngleDerivative(R):
    fills the face rotations Rf and the face rotation angle gradient in one
    loop. No caller has been found in the history of beamFoam.
*/
void interpolateRotationMatrix
(
    const volTensorField& R,
    surfaceTensorField& Rf,
    surfaceVectorField& snGradRotAngle
)
{
    forAll(Rf.boundaryField(), patchI)
    {
        Rf.boundaryFieldRef()[patchI] = R.boundaryField()[patchI];//here
    }

    const tensorField& RI = R.internalField();

    tensorField& RfI = Rf.primitiveFieldRef();
    vectorField& snGradRotAngleI = snGradRotAngle;

    const fvMesh& mesh = R.mesh();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs();

    forAll(snGradRotAngleI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];

        // scalar s = 0.5/deltaCoeffsI[faceI];
        RfI[faceI] = RI[own[faceI]] & rotationMatrix(0.5*DRotAngle);
    }

    forAll(snGradRotAngle.boundaryField(), patchI)
    {
        const tensorField& RI = R.boundaryField()[patchI];

        // tensorField& RfI = Rf.boundaryField()[patchI];
        vectorField& snGradRotAngleI =
            snGradRotAngle.boundaryFieldRef()[patchI];

        const labelList& own =
            mesh.boundary()[patchI].faceCells();
        const scalarField& deltaCoeffsI =
            mesh.deltaCoeffs().boundaryField()[patchI];

        forAll(snGradRotAngleI, faceI)
        {
            tensor DR = (R[own[faceI]].T() & RI[faceI]);

            vector DRotAngle = pseudoVector(DR);

            snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];
        }
    }
}


void interpolateRotationMatrix
(
    const beamModel& bm,
    const surfaceTensorField& Rf,
    volTensorField& R
)
{
    forAll(R.boundaryField(), patchI)
    {
        if (!R.boundaryField()[patchI].coupled())
        {
            R.boundaryFieldRef()[patchI] = Rf.boundaryField()[patchI];
        }
    }

    const tensorField& RfI = Rf.internalField();

    const labelList& upperNeiFaces = bm.upperNeiCellFaces();
    const labelList& lowerNeiFaces = bm.lowerNeiCellFaces();

    // Pout << "y1" << endl;
    // sleep(5);

    tensorField& RI = R.primitiveFieldRef();

    const fvMesh& mesh = R.mesh();
    // const unallocLabelList& own = mesh.owner();
    // const unallocLabelList& nei = mesh.neighbour();

    forAll(RI, cellI)
    {
        if (lowerNeiFaces[cellI] < 0)
        {
            label startPatchIndex = -lowerNeiFaces[cellI] - 1;

            // if (Rf.boundaryField()[startPatchIndex].size() > 1)
            // {
            //     FatalErrorIn
            //     (
            //         "interpolateRotationMatrix(...)"
            //     )
            //       << "Face to cell rotation matrix interpolation "
            //       << "is not implemented for multi-beam problem"
            //       << abort(FatalError);
            // }

            label bFaceI =
                findIndex
                (
                    mesh.boundary()[startPatchIndex].faceCells(),
                    cellI
                );

            tensor DR =
            (
                Rf.boundaryField()[startPatchIndex][bFaceI].T()
              & RfI[upperNeiFaces[cellI]]
            );

            vector DRotAngle = pseudoVector(DR);

            RI[cellI] =
            (
                Rf.boundaryField()[startPatchIndex][0].T()
              & rotationMatrix(0.5*DRotAngle)
            );
        }
        else if (upperNeiFaces[cellI] < 0)
        {
            label endPatchIndex = -upperNeiFaces[cellI] - 1;

            // if (Rf.boundaryField()[endPatchIndex].size() > 1)
            // {
            //     FatalErrorIn
            //     (
            //         "interpolateRotationMatrix(...)"
            //     )
            //       << "Face to cell rotation matrix interpolation "
            //       << "is not implemented for mutiple-beam problem"
            //       << abort(FatalError);
            // }

            label bFaceI =
                findIndex
                (
                    mesh.boundary()[endPatchIndex].faceCells(),
                    cellI
                );

            tensor DR =
            (
                RfI[lowerNeiFaces[cellI]].T()
              & Rf.boundaryField()[endPatchIndex][bFaceI]
            );

            vector DRotAngle = pseudoVector(DR);

            RI[cellI] =
            (
                RfI[lowerNeiFaces[cellI]]
              & rotationMatrix(0.5*DRotAngle)
            );
        }
        else
        {
            tensor DR =
            (
                RfI[lowerNeiFaces[cellI]].T()
              & RfI[upperNeiFaces[cellI]]
            );

            vector DRotAngle = pseudoVector(DR);

            RI[cellI] =
            (
                RfI[lowerNeiFaces[cellI]]
              & rotationMatrix(0.5*DRotAngle)
            );
        }
    }
}

/*
    NOTE: currently unused - kept for reference.

    Face rotation angle gradient: the relative rotation vector between the
    owner and neighbour cells (or the boundary value) multiplied by the
    face deltaCoeffs. No caller has been found in the history of beamFoam.
*/
tmp<surfaceVectorField> rotationAngleDerivative
(
    const volTensorField& rotMat
)
{
    tmp<surfaceVectorField> tRotAngDeriv
    (
        new surfaceVectorField
        (
            IOobject
            (
                "rotAngleDerivative",
                rotMat.time().timeName(),
                rotMat.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotMat.mesh(),
            dimensionedVector("zero", dimless/dimLength, vector::zero)
        )
    );

    vectorField& snGradRotAngleI = tRotAngDeriv.ref().primitiveFieldRef();

    const tensorField& RI = rotMat.internalField();

    const fvMesh& mesh = rotMat.mesh();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs();

    forAll(snGradRotAngleI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];
    }

    forAll(tRotAngDeriv().boundaryField(), patchI)
    {
        const tensorField& RI = rotMat.boundaryField()[patchI];

        vectorField& snGradRotAngleI =
            tRotAngDeriv.ref().boundaryFieldRef()[patchI];

        const labelList& own =
            mesh.boundary()[patchI].faceCells();
        const scalarField& deltaCoeffsI =
            mesh.deltaCoeffs().boundaryField()[patchI];

        forAll(snGradRotAngleI, faceI)
        {
            tensor DR = (rotMat[own[faceI]].T() & RI[faceI]);

            vector DRotAngle = pseudoVector(DR);

            snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];
        }
    }

    return tRotAngDeriv;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
