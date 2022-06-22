/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Function
    multiplyCoeff

Description
    Helper function used when assembling finite volume solid mechanics block
    coupled stiffness matrices.

SourceFile
    multiplyCoeff.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef multiplyCoeff_H
#define multiplyCoeff_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // This function multiplies a vector by a 4th order tensor by a third
    // order tensor as encountered in the assembling of block coefficients
    // in a finite volume solid mechanics procedure. The result is placed
    // in the coeff variable.
    // The multiplication evaluates the following in index notation:
    //     coeff_ij = Sf_m C_mikl g_k delta_lj
    // where delta is kronecker delta
    void multiplyCoeff
    (
        tensor& coeff,
        const vector& Sf,
        const symmTensor4thOrder& C,
        const vector& g
    );

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
