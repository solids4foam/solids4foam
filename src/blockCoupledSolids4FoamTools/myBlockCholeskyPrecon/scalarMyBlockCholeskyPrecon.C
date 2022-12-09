/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#ifndef scalarMyBlockCholeskyPrecon_H
#define scalarMyBlockCholeskyPrecon_H

#include "myBlockCholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::myBlockCholeskyPrecon<scalar>::calcPreconDiag()
{
    // Precondition the diagonal
    if (matrix_.symmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();

        forAll (upper, coeffI)
        {
            preconDiag_[upperAddr[coeffI]] -=
                sqr(upper[coeffI])/preconDiag_[lowerAddr[coeffI]];
        }
    }
    else if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        forAll (upper, coeffI)
        {
            preconDiag_[upperAddr[coeffI]] -=
                upper[coeffI]*lower[coeffI]/preconDiag_[lowerAddr[coeffI]];
        }
    }

    // Invert the diagonal for future use
    forAll (preconDiag_, i)
    {
        preconDiag_[i] = 1.0/preconDiag_[i];
    }
}


template<>
void Foam::myBlockCholeskyPrecon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const
{
    forAll(x, i)
    {
        x[i] = b[i]*preconDiag_[i];
    }

    if (matrix_.symmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();

        forAll (upper, coeffI)
        {
            x[upperAddr[coeffI]] -=
                preconDiag_[upperAddr[coeffI]]*
                upper[coeffI]*x[lowerAddr[coeffI]];
        }

        forAllReverse (upper, coeffI)
        {
            x[lowerAddr[coeffI]] -=
                preconDiag_[lowerAddr[coeffI]]*
                upper[coeffI]*x[upperAddr[coeffI]];
        }
    }
    else if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
        const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        label losortCoeff;

        forAll (lower, coeffI)
        {
            losortCoeff = losortAddr[coeffI];

            x[upperAddr[losortCoeff]] -=
                preconDiag_[upperAddr[losortCoeff]]*
                lower[losortCoeff]*x[lowerAddr[losortCoeff]];
        }

        forAllReverse (upper, coeffI)
        {
            x[lowerAddr[coeffI]] -=
                preconDiag_[lowerAddr[coeffI]]*
                upper[coeffI]*x[upperAddr[coeffI]];
        }
    }
}


template<>
void Foam::myBlockCholeskyPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const
{
    if (matrix_.symmetric())
    {
        precondition(xT, bT);
    }

    forAll(xT, i)
    {
        xT[i] = bT[i]*preconDiag_[i];
    }

    if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
        const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        label losortCoeff;

        forAll (lower, coeffI)
        {
            xT[upperAddr[coeffI]] -=
                preconDiag_[upperAddr[coeffI]]*
                upper[coeffI]*xT[lowerAddr[coeffI]];
        }

        forAllReverse (upper, coeffI)
        {
            losortCoeff = losortAddr[coeffI];

            xT[lowerAddr[losortCoeff]] -=
                preconDiag_[lowerAddr[losortCoeff]]*
                lower[losortCoeff]*xT[upperAddr[losortCoeff]];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
