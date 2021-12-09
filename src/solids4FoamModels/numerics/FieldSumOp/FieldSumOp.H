/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Class
    Foam::FieldSumOp

Description
    Op function class to perform operations on a list.

    This is required for clang, where reduce operations on lists are ambiguous.

Author
    Philip Cardiff, UCD.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef FieldSumOp_H
#define FieldSumOp_H

#include "Field.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class Type>
    class FieldSumOp
    {
        public:

        // Op function
        tmp< Field<Type> > operator()
        (
            const Field<Type>& x, const Field<Type>& y
        ) const
        {
            tmp<Field<Type> > tresult(new Field<Type>(x));
#ifdef OPENFOAMESIORFOUNDATION
            Field<Type>& result = tresult.ref();
#else
            Field<Type>& result = tresult();
#endif

            forAll(result, i)
            {
                result[i] += y[i];
            }

            return tresult;
        }
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
