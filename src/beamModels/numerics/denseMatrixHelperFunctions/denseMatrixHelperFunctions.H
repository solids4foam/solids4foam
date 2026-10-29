/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Helper functions for the dense matrix class

Author
    Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#ifndef denseMatrixHelperFunctions_H
#define denseMatrixHelperFunctions_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    tmp<scalarField> magSqr(const Field<scalarRectangularMatrix>& mats);

// tmp<scalarRectangularMatrix> cmptDivide(const Field<scalarRectangularMatrix>& mat1 , const Field<scalarRectangularMatrix>& mat2 )
// {
//     tmp<scalarRectangularMatrix> outputMat(new scalarRectangularMatrix(mat1.size(), 0.0));
//     Field<scalarRectangularMatrix> result = outputMat.ref();

//     forAll(mat1, i)
//     {
//         result(i,0) = mat1(i,0)/mat2(i,0);
//     }

//     return outputMat;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
