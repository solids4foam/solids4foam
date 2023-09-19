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
    MyBlockCholeskyPrecon

Description
    Template specialisation for scalar block Cholesky preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    scalarMyBlockCholeskyPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef scalarMyBlockCholeskyPrecon_H
#define scalarMyBlockCholeskyPrecon_H

#include "myBlockCholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::myBlockCholeskyPrecon<scalar>::calcPreconDiag();

template<>
void Foam::myBlockCholeskyPrecon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const;


template<>
void Foam::myBlockCholeskyPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
