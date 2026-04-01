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

#include "symmTensor3rdOrder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const symmTensor3rdOrder::typeName = "symmTensor3rdOrder";

template<>
const char* symmTensor3rdOrder::componentNames[] =
{
    "xxx",
    "xxy",
    "xxz",
    "xyy",
    "xyz",
    "xzz",
    "yyy",
    "yyz",
    "yzz",
    "zzz"
};

template<>
const symmTensor3rdOrder symmTensor3rdOrder::zero
(
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0
);

template<>
const symmTensor3rdOrder symmTensor3rdOrder::one
(
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
);

template<>
const symmTensor3rdOrder symmTensor3rdOrder::max
(
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT,
    VGREAT, VGREAT, VGREAT, VGREAT, VGREAT
);

template<>
const symmTensor3rdOrder symmTensor3rdOrder::min
(
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT,
    -VGREAT, -VGREAT, -VGREAT, -VGREAT, -VGREAT
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
