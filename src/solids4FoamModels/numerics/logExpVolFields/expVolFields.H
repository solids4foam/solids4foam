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

Description
    Calculate the exponetial of a volSymmTensorField.

    To calculate the exp of a tensor, we must rotate the tensor to principal
    components and calculate the exp of the principals components, then
    rotate these principal exp components back to get the exp tensor.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef expVolSymmTensorField_H
#define expVolSymmTensorField_H

#include "eig3Field.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Exponetial of a symmTensor field
    tmp<volSymmTensorField> exp(const volSymmTensorField& vf);

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
