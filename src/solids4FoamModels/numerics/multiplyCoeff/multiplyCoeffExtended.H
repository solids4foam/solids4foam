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
    multiplyCoeffExtended

Description
    Helper function used when assembling finite volume solid mechanics block
    coupled stiffness matrices. Extended to include sensitivity term for
    large deformations.

SourceFile
    multiplyCoeffExtended.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Federico Mazzanti, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef multiplyCoeffExtended_H
#define multiplyCoeffExtended_H

#include "fieldTypes.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // The multiplication evaluates the following in index notation:
    //     coeff_ij = (Sf_m C_mikl + G_mkl sigma_mi)*g_k delta_lj
    // where delta is kronecker delta
    void multiplyCoeffExtended
    (
        tensor& coeff,
        const vector& Sf,
        const RectangularMatrix<scalar>& C,
        const RectangularMatrix<scalar>& G,
        const symmTensor& sigma,
        const vector& g
    );

    // The multiplication evaluates the following in index notation:
    //     coeff_i = Sf_m P_mi where P = d(s - pI)/dp
    void multiplyCoeffPressure
	(
		vector& coeff,
		const vector& Sf,
		const symmTensor& P
	);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
