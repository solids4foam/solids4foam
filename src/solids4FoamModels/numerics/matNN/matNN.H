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
    Foam::matNN

  Description
      Lightweight fixed-size matrix types (row-major storage) for material
      tangent operators or similar data in Voigt notation.

      - Provides (i,j) element accessors.
      - Provides m() and n() methods for compatibility with existing code.
      - Includes helpers to zero, fill with a constant, or set to identity.
      - Supports initialisation from an initializer_list.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef matNN_H
#define matNN_H

#include "mat39.H"
#include "mat66.H"

#endif

// ************************************************************************* //
