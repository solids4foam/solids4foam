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

#include "spinTensor.H"
#include "emptyFvPatchFields.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<surfaceTensorField> spinTensor(const surfaceVectorField& axialVector)
{
    tmp<surfaceTensorField> tresult
    (
        new surfaceTensorField
        (
            IOobject
            (
                "spin("+axialVector.name()+")",
                axialVector.time().timeName(),
                axialVector.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            axialVector.mesh(),
            dimensionedTensor("0", axialVector.dimensions(), tensor::zero)
        )
    );

    surfaceTensorField& result = tresult.ref();
    tensorField& resultI = result.primitiveFieldRef();

    const vectorField& axialVectorI = axialVector.internalField();

    forAll(resultI, faceI)
    {
        resultI[faceI] = spinTensor(axialVectorI[faceI]);
    }

    forAll(result.boundaryField(), patchI)
    {
        if
        (
            axialVector.boundaryField()[patchI].type()
            != emptyFvPatchField<vector>::typeName
        )
        {
            tensorField& pResult = result.boundaryFieldRef()[patchI];
            const vectorField& pAxialVector =
                axialVector.boundaryField()[patchI];

            forAll(pResult, faceI)
            {
                pResult[faceI] = spinTensor(pAxialVector[faceI]);
            }
        }
    }

    return tresult;
}

tmp<volTensorField> spinTensor(const volVectorField& axialVector)
{
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "spin("+axialVector.name()+")",
                axialVector.time().timeName(),
                axialVector.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            axialVector.mesh(),
            dimensionedTensor("0", axialVector.dimensions(), tensor::zero)
        )
    );

    volTensorField& result = tresult.ref();
    tensorField& resultI = result.primitiveFieldRef();

    const vectorField& axialVectorI = axialVector.internalField();

    forAll(resultI, cellI)
    {
        resultI[cellI] = spinTensor(axialVectorI[cellI]);
    }

    forAll(result.boundaryField(), patchI)
    {
        if
        (
            axialVector.boundaryField()[patchI].type()
            != emptyFvPatchField<vector>::typeName
        )
        {
            tensorField& pResult = result.boundaryFieldRef()[patchI];
            const vectorField& pAxialVector =
                axialVector.boundaryField()[patchI];

            forAll(pResult, faceI)
            {
                pResult[faceI] = spinTensor(pAxialVector[faceI]);
            }
        }
    }

    return tresult;
}

tmp<tensorField> spinTensor(const vectorField& axialVector)
{
    tmp<tensorField> tresult
    (
        new tensorField(axialVector.size(), tensor::zero)
    );

    tensorField& result = tresult.ref();

    forAll(result, cellI)
    {
        result[cellI] = spinTensor(axialVector[cellI]);
    }

    return tresult;
}

tensor spinTensor(const vector& v)
{
    tensor result = tensor::zero;

    result.xy() = -v.z();
    result.xz() = v.y();

    result.yx() = v.z();
    result.yz() = -v.x();

    result.zx() = -v.y();
    result.zy() = v.x();

    return result;
}


tmp<surfaceVectorField> axialVector(const surfaceTensorField& spinTensor)
{
    tmp<surfaceVectorField> tresult
    (
        new surfaceVectorField
        (
            IOobject
            (
                "axial"+spinTensor.name()+")",
                spinTensor.time().timeName(),
                spinTensor.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            spinTensor.mesh(),
            dimensionedVector("0", spinTensor.dimensions(), vector::zero)
        )
    );

    surfaceVectorField& result = tresult.ref();
    vectorField& resultI = result.primitiveFieldRef();

    const tensorField& spinTensorI = spinTensor.internalField();

    forAll(resultI, faceI)
    {
        resultI[faceI] = axialVector(spinTensorI[faceI]);
    }

    forAll(result.boundaryField(), patchI)
    {
        if
        (
            spinTensor.boundaryField()[patchI].type()
            != emptyFvPatchField<tensor>::typeName
        )
        {
            vectorField& pResult = result.boundaryFieldRef()[patchI];
            const tensorField& pSpinTensor =
                spinTensor.boundaryField()[patchI];

            forAll(pResult, faceI)
            {
                pResult[faceI] = axialVector(pSpinTensor[faceI]);
            }
        }
    }

    return tresult;
}

tmp<volVectorField> axialVector(const volTensorField& spinTensor)
{
    tmp<volVectorField> tresult
    (
        new volVectorField
        (
            IOobject
            (
                "axial("+spinTensor.name()+")",
                spinTensor.time().timeName(),
                spinTensor.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            spinTensor.mesh(),
            dimensionedVector("0", spinTensor.dimensions(), vector::zero)
        )
    );

    volVectorField& result = tresult.ref();
    vectorField& resultI = result.primitiveFieldRef();

    const tensorField& spinTensorI = spinTensor.internalField();

    forAll(resultI, cellI)
    {
        resultI[cellI] = axialVector(spinTensorI[cellI]);
    }

    forAll(result.boundaryField(), patchI)
    {
        if
        (
            spinTensor.boundaryField()[patchI].type()
            != emptyFvPatchField<tensor>::typeName
        )
        {
            vectorField& pResult = result.boundaryFieldRef()[patchI];
            const tensorField& pSpinTensor =
                spinTensor.boundaryField()[patchI];

            forAll(pResult, faceI)
            {
                pResult[faceI] = axialVector(pSpinTensor[faceI]);
            }
        }
    }

    return tresult;
}

tmp<vectorField> axialVector(const tensorField& spinTensor)
{
    tmp<vectorField> tresult
    (
        new vectorField(spinTensor.size(), vector::zero)
    );

    vectorField& result = tresult.ref();

    forAll(result, cellI)
    {
        result[cellI] = axialVector(spinTensor[cellI]);
    }

    return tresult;
}

vector axialVector(const tensor& T)
{
    vector result = vector::zero;

    result.x() = T.zy();
    result.y() = T.xz();
    result.z() = T.yx();

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
