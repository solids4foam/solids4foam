/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "labelSymmTensor4thOrder.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const labelSymmTensor4thOrder::typeName = "labelSymmTensor4thOrder";

template<>
const char* labelSymmTensor4thOrder::componentNames[] =
{
    "xxxx", "xxyy", "xxzz",
            "yyyy", "yyzz",
                    "zzzz",
                           "xyxy",
                                  "yzyz",
                                         "zxzx"
};

template<>
const labelSymmTensor4thOrder labelSymmTensor4thOrder::zero
(
    0, 0, 0,
       0, 0,
          0,
            0,
              0,
                0
);

template<>
const labelSymmTensor4thOrder labelSymmTensor4thOrder::one
(
    1, 1, 1,
       1, 1,
          1,
            1,
              1,
                1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
