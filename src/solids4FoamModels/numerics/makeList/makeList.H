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
    Foam::makeList

Description
      Provides a helper function Foam::makeList<T>() for constructing
      Foam::List<T> objects from an initializer_list.

      This utility enables code such as

          labelList ids = Foam::makeList<label>({0,1,2});

      which is supported in modern OpenFOAM via public initializer_list
      constructors, but not in foam-extend. Using makeList() centralises
      the difference and provides a single, portable syntax across all
      OpenFOAM variants.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef makeList_H
#define makeList_H

#include "List.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class T>
    inline List<T> makeList(std::initializer_list<T> il)
    {
        List<T> L(il.size());
        label i = 0;

        for (const T& v : il)
        {
            L[i++] = v;
        }

        return L; // NRVO/move on modern forks; fine for small lists regardless
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
