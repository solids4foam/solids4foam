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
    Foam::ListMaxOp

Description
    Op function class to perform operations on a list.

Author
    Philip Cardiff, UCD.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef ListMaxOp_H
#define ListMaxOp_H

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class Type>
    struct ListMaxOp
    {
        public:

        // Op function
        List<Type> operator()
        (
            const List<Type>& x, const List<Type>& y
        ) const
        {
            if (x.size() != y.size())
            {
                FatalErrorIn
                (
                    "List<Type> ListMaxOp()"
                    "(const List<Type>& x, const List<Type>& y) const"
                )   << "The list sizes do not match! "
                    << "x: " << x.size() << ", y: " << y.size()
                    << abort(FatalError);
            }

            List<Type> result(x.size());

            forAll(result, i)
            {
                result[i] = max(x[i], y[i]);
            }

            return result;
        }
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
