
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

Class
    WendlandC0Function

Description
    Wendland C0 radial basis function.

Author
    David Blom, TU Delft. All rights reserved.
    NOTE: This file is distributed under the GNU General Public License (GPL).
    The phrase "All rights reserved" above is retained from the original
    source but does not limit the permissions granted by the GPL.

\*---------------------------------------------------------------------------*/

#ifndef WendlandC0Function_H
#define WendlandC0Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC0Function : public RBFFunctionInterface
    {
        public:
            explicit WendlandC0Function( scalar radius );

            virtual ~WendlandC0Function();

            virtual scalar evaluate( scalar value );

            scalar radius;
    };
}

#endif
