
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
    RBFFunctionInterface

Description
    Abstract interface for radial basis functions.

Author
    David Blom, TU Delft. All rights reserved.
    NOTE: This file is distributed under the GNU General Public License (GPL).
    The phrase "All rights reserved" above is retained from the original
    source but does not limit the permissions granted by the GPL.

\*---------------------------------------------------------------------------*/

#ifndef RBFFunctionInterface_H
#define RBFFunctionInterface_H

#include <math.h>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include "fvCFD.H"

namespace rbf
{
    class RBFFunctionInterface
    {
        public:
            virtual ~RBFFunctionInterface(){}

            virtual scalar evaluate( scalar value ) = 0;
    };
}

#endif
