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

#include "eig3Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void eig3Field
(
    const volTensorField& A, volTensorField& V, volVectorField& d
)
{
    const tensorField& AI = A.internalField();
#ifdef OPENFOAMESIORFOUNDATION
    tensorField& VI = V.primitiveFieldRef();
    vectorField& dI = d.primitiveFieldRef();
#else
    tensorField& VI = V.internalField();
    vectorField& dI = d.internalField();
#endif

    forAll(AI, cellI)
    {
        eig3().eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if (A.boundaryField()[patchI].type() != "empty")
        {
            const tensorField& AB = A.boundaryField()[patchI];
#ifdef OPENFOAMESIORFOUNDATION
            tensorField& VB = V.boundaryFieldRef()[patchI];
            vectorField& dB = d.boundaryFieldRef()[patchI];
#else
            tensorField& VB = V.boundaryField()[patchI];
            vectorField& dB = d.boundaryField()[patchI];
#endif

            forAll(AB, faceI)
            {
                eig3().eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }
}


void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volVectorField& d
)
{
    const symmTensorField& AI = A.internalField();
#ifdef OPENFOAMESIORFOUNDATION
    tensorField& VI = V.primitiveFieldRef();
    vectorField& dI = d.primitiveFieldRef();
#else
    tensorField& VI = V.internalField();
    vectorField& dI = d.internalField();
#endif

    forAll(AI, cellI)
    {
        eig3().eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if (A.boundaryField()[patchI].type() != "empty")
        {
            const symmTensorField& AB = A.boundaryField()[patchI];
#ifdef OPENFOAMESIORFOUNDATION
            tensorField& VB = V.boundaryFieldRef()[patchI];
            vectorField& dB = d.boundaryFieldRef()[patchI];
#else
            tensorField& VB = V.boundaryField()[patchI];
            vectorField& dB = d.boundaryField()[patchI];
#endif

            forAll(AB, faceI)
            {
                eig3().eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }
}


#ifndef OPENFOAMESIORFOUNDATION
void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volDiagTensorField& d
)
{
    const symmTensorField& AI = A.internalField();
    tensorField& VI = V.internalField();
    diagTensorField& dI = d.internalField();

    forAll(AI, cellI)
    {
        eig3().eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if (A.boundaryField()[patchI].type() != "empty")
        {
            const symmTensorField& AB = A.boundaryField()[patchI];
            tensorField& VB = V.boundaryField()[patchI];
            diagTensorField& dB = d.boundaryField()[patchI];

            forAll(AB, faceI)
            {
                eig3().eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }
}
#endif


void eig3Field
(
    const surfaceSymmTensorField& A,
    surfaceTensorField& V,
    surfaceVectorField& d
)
{
    const symmTensorField& AI = A.internalField();
#ifdef OPENFOAMESIORFOUNDATION
    tensorField& VI = V.primitiveFieldRef();
    vectorField& dI = d.primitiveFieldRef();
#else
    tensorField& VI = V.internalField();
    vectorField& dI = d.internalField();
#endif

    forAll(AI, cellI)
    {
        eig3().eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if (A.boundaryField()[patchI].type() != "empty")
        {
            const symmTensorField& AB = A.boundaryField()[patchI];
#ifdef OPENFOAMESIORFOUNDATION
            tensorField& VB = V.boundaryFieldRef()[patchI];
            vectorField& dB = d.boundaryFieldRef()[patchI];
#else
            tensorField& VB = V.boundaryField()[patchI];
            vectorField& dB = d.boundaryField()[patchI];
#endif

            forAll(AB, faceI)
            {
                eig3().eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
