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


Class
    doubleDotProduct

Description
    Helper function to mimic the double dot product between a 4th order tensor
    and a 2nd order tensor, as natively implemented in foam-extend.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef doubleDotProduct_H
#define doubleDotProduct_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef OPENFOAMESIORFOUNDATION

namespace Foam
{

    //- Calculate the stress = C && epsilon, where C is stored as a tensor
    void doubleDotProduct(symmTensor& s, const tensor& C, const symmTensor& e)
    {
        s[symmTensor::XX] =
            C[0]*e[symmTensor::XX]
          + C[1]*e[symmTensor::YY]
          + C[2]*e[symmTensor::ZZ];

        s[symmTensor::YY] =
            C[1]*e[symmTensor::XX]
          + C[3]*e[symmTensor::YY]
          + C[4]*e[symmTensor::ZZ];

        s[symmTensor::ZZ] =
            C[2]*e[symmTensor::XX]
          + C[4]*e[symmTensor::YY]
          + C[5]*e[symmTensor::ZZ];

        s[symmTensor::XY] = C[6]*e[symmTensor::XY];
        s[symmTensor::YZ] = C[7]*e[symmTensor::YZ];
        s[symmTensor::XZ] = C[8]*e[symmTensor::XZ];
    }


    //- Calculate the stress = C && epsilon, where C is stored as a tensor
    void doubleDotProduct
    (
        volSymmTensorField& sigma,
        const volTensorField& C,
        const volSymmTensorField& epsilon
    )
    {
        // Take references
        symmTensorField& sigmaI = sigma.primitiveFieldRef();
        const tensorField& CI = C.primitiveField();
        const symmTensorField& epsilonI = epsilon.primitiveField();

        // Calculate internal field
        forAll(sigmaI, cellI)
        {
            doubleDotProduct(sigmaI[cellI], CI[cellI], epsilonI[cellI]);
        }

        // Calculate boundary field
        forAll(sigma.boundaryField(), patchI)
        {
            // Take references
            symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
            const tensorField& CP = C.boundaryField()[patchI];
            const symmTensorField& epsilonP = epsilon.boundaryField()[patchI];

            // Calculate internal field
            forAll(sigmaP, faceI)
            {
                doubleDotProduct(sigmaP[faceI], CP[faceI], epsilonP[faceI]);
            }
        }
    }


    //- Calculate the stress = C && epsilon, where C is stored as a tensor
    void doubleDotProduct
    (
        surfaceSymmTensorField& sigma,
        const surfaceTensorField& C,
        const surfaceSymmTensorField& epsilon
    )
    {
        // Take references
        symmTensorField& sigmaI = sigma.primitiveFieldRef();
        const tensorField& CI = C.primitiveField();
        const symmTensorField& epsilonI = epsilon.primitiveField();

        // Calculate internal field
        forAll(sigmaI, faceI)
        {
            doubleDotProduct(sigmaI[faceI], CI[faceI], epsilonI[faceI]);
        }

        // Calculate boundary field
        forAll(sigma.boundaryField(), patchI)
        {
            // Take references
            symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
            const tensorField& CP = C.boundaryField()[patchI];
            const symmTensorField& epsilonP = epsilon.boundaryField()[patchI];

            // Calculate internal field
            forAll(sigmaP, faceI)
            {
                doubleDotProduct(sigmaP[faceI], CP[faceI], epsilonP[faceI]);
            }
        }
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // ifdef OPENFOAMESIORFOUNDATION

#endif

// ************************************************************************* //
