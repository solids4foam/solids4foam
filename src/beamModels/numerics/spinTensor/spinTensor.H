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

Description
    This code is used for:
    1. Construction of a 3 x 3 skew-symmetric matrix, i.e., A.T() = - A,
       from a 3 x 1 vector (member function 'spinTensor').
    2. Construction of a 3 x 1 vector (also called as axial vector) from
       a skew-symmetric matrix (member function 'axialVector')

Author
    Zeljko Tukovic, FSB Zagreb,
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef spinTensor_H
#define spinTensor_H

#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<surfaceTensorField> spinTensor(const surfaceVectorField& axialVector);
tmp<volTensorField> spinTensor(const volVectorField& axialVector);
tmp<tensorField> spinTensor(const vectorField& axialVector);

tensor spinTensor(const vector& axialVector);

tmp<surfaceVectorField> axialVector(const surfaceTensorField& spinTensor);
tmp<volVectorField> axialVector(const volTensorField& spinTensor);
tmp<vectorField> axialVector(const tensorField& spintTensor);

vector axialVector(const tensor& spinTensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
