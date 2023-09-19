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
    BlockCholeskyPrecon

Description
    Template specialisation for tensor block Cholesky preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    tensorBlockCholeskyPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef tensorMyBlockCholeskyPrecon_H
#define tensorMyBlockCholeskyPrecon_H

#include "myBlockCholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::myBlockCholeskyPrecon<tensor>::calcPreconDiag();

template<>
void Foam::myBlockCholeskyPrecon<tensor>::precondition
(
    tensorField& x,
    const tensorField& b
) const;


template<>
void Foam::myBlockCholeskyPrecon<tensor>::preconditionT
(
    tensorField& xT,
    const tensorField& bT
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
