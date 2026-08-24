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

Typedef
    Foam::symmTensor3rdOrder

Description
    SymmTensor3rdOrder of scalars.

SourceFiles
    symmTensor3rdOrder.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensor3rdOrder_H
#define symmTensor3rdOrder_H

#include "SymmTensor3rdOrderTemplate.H"
#include "vector.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef SymmTensor3rdOrder<scalar> symmTensor3rdOrder;

#ifdef OPENFOAM_COM

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<SymmTensor3rdOrder<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<SymmTensor3rdOrder<Cmpt>> :
    is_contiguous_label<Cmpt>
{};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<SymmTensor3rdOrder<Cmpt>> :
    is_contiguous_scalar<Cmpt>
{};

#else

//- Data associated with symmTensor3rdOrder type are contiguous
template<>
inline bool contiguous<symmTensor3rdOrder>() {return true;}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
