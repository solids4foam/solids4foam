/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Calculate log of a volSymmTensorField
    To calculate log of a tensor, we must rotate the tensor to principal
    components and calculate the log of the principals components, then
    rotate this principal log components back to get the log tensor

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "logVolFields.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


tmp<volSymmTensorField> log(const volSymmTensorField& vf)
{
    tmp<volSymmTensorField> tresult
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "log("+vf.name()+")",
                    vf.time().timeName(),
                    vf.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                    ),
                vf
                )
            );

    volSymmTensorField& result = tresult();

    // Calculate eigen values and eigen vectors
    // The OpenFOAM eigenValues/eigenVectors sometimes give wrong results when
    // eigenValues are repeated or zero, so I will use my own implementation.
    // The efficiency of the implementation may need to be revisited, however,
    // it is fine for creation of post processing fields e.g calculate true
    // strain

    volDiagTensorField eigenVal
        (
            IOobject
            (
                "eigenVal("+vf.name()+")",
                vf.time().timeName(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            dimensionedDiagTensor("zero", vf.dimensions(), diagTensor::zero)
        );

    volTensorField eigenVec
        (
            IOobject
            (
                "eigenVec("+vf.name()+")",
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

    // Now we will calculate the log og the eigenValues and then rotate the
    // tensor back to the physcial configuration

    const diagTensorField& eigenValI = eigenVal.internalField();
    const tensorField& eigenVecI = eigenVec.internalField();
    symmTensor logEigenVal = symmTensor::zero;
    symmTensorField& resultI = result.internalField();

    forAll(eigenValI, cellI)
    {
        // We cannot have log of zeros
        if
        (
            eigenValI[cellI][diagTensor::XX] < SMALL
            || eigenValI[cellI][diagTensor::YY] < SMALL
            || eigenValI[cellI][diagTensor::ZZ] < SMALL
        )
        {
            FatalError
                << "log of zero is not allowed"
                << abort(FatalError);
        }

        // Calculate log
        logEigenVal[symmTensor::XX] =
            Foam::log(eigenValI[cellI][diagTensor::XX]);
        logEigenVal[symmTensor::YY] =
            Foam::log(eigenValI[cellI][diagTensor::YY]);
        logEigenVal[symmTensor::ZZ] =
            Foam::log(eigenValI[cellI][diagTensor::ZZ]);

        // Rotate back
        resultI[cellI] = transform(eigenVecI[cellI].T(), logEigenVal);
    }

    forAll(eigenVal.boundaryField(), patchI)
    {
        if
        (
            !vf.boundaryField()[patchI].coupled()
            && vf.boundaryField()[patchI].type()
                != emptyFvPatchField<symmTensor>::typeName
        )
        {
            const diagTensorField& eigenValB = eigenVal.boundaryField()[patchI];
            const tensorField& eigenVecB = eigenVec.boundaryField()[patchI];
            symmTensorField& resultB = result.boundaryField()[patchI];

            forAll(eigenValB, faceI)
            {
                // We cannot have log of zeros
                if
                (
                    eigenValB[faceI][diagTensor::XX] < SMALL
                    || eigenValB[faceI][diagTensor::YY] < SMALL
                    || eigenValB[faceI][diagTensor::ZZ] < SMALL
                )
                {
                    FatalError
                        << "log of zero is not allowed"
                        << abort(FatalError);
                }

                // Calculate log
                logEigenVal[symmTensor::XX] =
                    Foam::log(eigenValB[faceI][diagTensor::XX]);
                logEigenVal[symmTensor::YY] =
                    Foam::log(eigenValB[faceI][diagTensor::YY]);
                logEigenVal[symmTensor::ZZ] =
                    Foam::log(eigenValB[faceI][diagTensor::ZZ]);

                // Rotate back
                resultB[faceI] = transform(eigenVecB[faceI].T(), logEigenVal);
            }
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


} // End namespace Foam


// ************************************************************************* //
