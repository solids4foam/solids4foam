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

Function
    multiplyCoeff

Description
    Helper function used when assembling finite volume solid mechanics block
    coupled stiffness matrices.

SourceFile
    multiplyCoeff.C

Author
    Federico Mazzanti, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef multiplyCoeff_H
#define multiplyCoeff_H

#include "fieldTypes.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "scalarMatrices.H"
#else
    #include "scalarSquareMatrix.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // This function multiplies a vector by a 4th order tensor by a third
    // order tensor as encountered in the assembling of block coefficients
    // in a finite volume solid mechanics procedure. The result is placed
    // in the coeff variable.
    // The 4th order material tangent matrix is represented by a 6x6 matrix.
    // The multiplication evaluates the following in index notation:
    //     coeff_ij = Sf_m C_mikl g_k delta_lj
    // where delta is kronecker delta
    void multiplyCoeff
    (
        tensor& coeff,
        const vector& Sf,
        const scalarSquareMatrix& C,
        const vector& g
    );

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
