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
    BlockLduMatrix

Description
    Template specialisation for vector block matrix

Author
    Philip Cardiff

SourceFiles
    vectorBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef vectorBlockLduMatrix_H
#define vectorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// template<>
// void BlockLduMatrix<vector>::sumDiag();

// template<>
// void BlockLduMatrix<vector>::negSumDiag();

// template<>
// void BlockLduMatrix<vector>::check() const;

// template<>
// void BlockLduMatrix<vector>::relax
// (
//     const vectorField& x,
//     vectorField& b,
//     const scalar alpha
// );

// template<>
// void BlockLduMatrix<vector>::operator*=(const scalarField& sf);

// template<>
// void BlockLduMatrix<vector>::AmulCore
// (
//     vectorField& mul,
//     const vectorField& x
// ) const;

// template<>
// void BlockLduMatrix<vector>::TmulCore
// (
//     vectorField& mul,
//     const vectorField& x
// ) const;

// template<>
// void BlockLduMatrix<vector>::segregateB
// (
//     vectorField& mul,
//     const vectorField& x
// ) const;


template<>
void BlockLduMatrix<vector>::Amul
(
    Field<vector>& Ax,
    const Field<vector>& x
) const;


// template<>
// void BlockLduMatrix<vector>::Tmul
// (
//     Field<vector>& Ax,
//     const Field<vector>& x
// ) const;


// template<>
// tmp<vectorField> BlockLduMatrix<vector>::H(const vectorField& x) const;

// template<>
// tmp<vectorField> BlockLduMatrix<vector>::faceH(const vectorField& x) const;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
