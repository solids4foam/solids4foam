/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Calculate the exponetial of a volSymmTensorField.

    To calculate the exp of a tensor, we must rotate the tensor to principal
    components and calculate the exp of the principals components, then
    rotate these principal exp components back to get the exp tensor.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "expVolFields.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


tmp<volSymmTensorField> exp(const volSymmTensorField& vf)
{
    tmp<volSymmTensorField> tresult
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "exp("+vf.name()+")",
                    vf.time().timeName(),
                    vf.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                vf.mesh(),
                dimensionedSymmTensor("zero", vf.dimensions(), symmTensor::zero)
            )
        );

#ifdef OPENFOAMESIORFOUNDATION
    volSymmTensorField& result = tresult.ref();
#else
    volSymmTensorField& result = tresult();
#endif

    // Calculate eigen values and eigen vectors
    // The OpenFOAM eigenValues/eigenVectors sometimes gives wrong results, when
    // eigenValues are repeated or zero, so I will use my own implementation.
    // The efficiency of the implementation may need to be revisited, however,
    // it is fine for creation of post processing fields e.g calculate true
    // strain etc.

    // Eigen value field
    // We will store the eigen values in a vector instead of a diagTensor
    // because the tranform function is not definite for diagTensors on a wedge
    // boundary
    volVectorField eigenVal
    (
        IOobject
        (
            "eigenVal(" + vf.name() + ")",
            vf.time().timeName(),
            vf.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vf.mesh(),
        dimensionedVector("zero", vf.dimensions(), vector::zero)
    );

    // Eigen vectors will be store in the rows i.e. the first eigen vector
    // is (eigenVec.xx() eigenVec.xy() eigenVec.xz())
    volTensorField eigenVec
    (
        IOobject
        (
            "eigenVec(" + vf.name() + ")",
            vf.time().timeName(),
            vf.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vf.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    // Calculate eigen values and eigen vectors of vf
    eig3Field(vf, eigenVec, eigenVal);

    // Now we will calculate the exp of the eigenValues and then rotate the
    // tensor back to the physcial configuration

    const vectorField& eigenValI = eigenVal.internalField();
    const tensorField& eigenVecI = eigenVec.internalField();
    symmTensor expEigenVal = symmTensor::zero;

#ifdef OPENFOAMESIORFOUNDATION
    symmTensorField& resultI = result.primitiveFieldRef();
#else
    symmTensorField& resultI = result.internalField();
#endif

    forAll(eigenValI, cellI)
    {
        // Calculate exp
        expEigenVal[symmTensor::XX] =
            Foam::exp(eigenValI[cellI][vector::X]);
        expEigenVal[symmTensor::YY] =
            Foam::exp(eigenValI[cellI][vector::Y]);
        expEigenVal[symmTensor::ZZ] =
            Foam::exp(eigenValI[cellI][vector::Z]);

        // Rotate back
        resultI[cellI] = transform(eigenVecI[cellI].T(), expEigenVal);
    }

    forAll(eigenVal.boundaryField(), patchI)
    {
        if
        (
            vf.boundaryField()[patchI].type()
         != emptyFvPatchField<symmTensor>::typeName
        )
        {
            const vectorField& eigenValB = eigenVal.boundaryField()[patchI];
            const tensorField& eigenVecB = eigenVec.boundaryField()[patchI];
#ifdef OPENFOAMESIORFOUNDATION
            symmTensorField& resultB = result.boundaryFieldRef()[patchI];
#else
            symmTensorField& resultB = result.boundaryField()[patchI];
#endif

            forAll(eigenValB, faceI)
            {
                // Calculate exp
                expEigenVal[symmTensor::XX] =
                    Foam::exp(eigenValB[faceI][vector::X]);
                expEigenVal[symmTensor::YY] =
                    Foam::exp(eigenValB[faceI][vector::Y]);
                expEigenVal[symmTensor::ZZ] =
                    Foam::exp(eigenValB[faceI][vector::Z]);

                // Rotate back
                resultB[faceI] = transform(eigenVecB[faceI].T(), expEigenVal);
            }
        }
    }

    return tresult;
}


} // End namespace Foam

// ************************************************************************* //
