/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

\*---------------------------------------------------------------------------*/

#include "symmTensor4thOrder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const symmTensor4thOrder::typeName = "symmTensor4thOrder";

template<>
const char* symmTensor4thOrder::componentNames[] =
{
    "xxxx", "xxyy", "xxzz",
            "yyyy", "yyzz",
                    "zzzz",
                           "xyxy",
                                  "yzyz",
                                         "zxzx"
};

template<>
const symmTensor4thOrder symmTensor4thOrder::zero
(
    0, 0, 0,
       0, 0,
          0,
            0,
              0,
                0
);

template<>
const symmTensor4thOrder symmTensor4thOrder::one
(
    1, 1, 1,
       1, 1,
          1,
            1,
              1,
                1
);

template<>
const symmTensor4thOrder symmTensor4thOrder::max
(
    VGREAT, VGREAT, VGREAT,
            VGREAT, VGREAT,
                    VGREAT,
                           VGREAT,
                                  VGREAT,
                                         VGREAT
);

template<>
const symmTensor4thOrder symmTensor4thOrder::min
(
    -VGREAT, -VGREAT, -VGREAT,
             -VGREAT, -VGREAT,
                      -VGREAT,
                              -VGREAT,
                                      -VGREAT,
                                              -VGREAT
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
