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

#include "error.H"
#include "myBlockCholeskyPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::myBlockCholeskyPrecon<Type>::calcPreconDiag()
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
    }

    // Invert the diagonal
    preconDiag_ = inv(preconDiag_);
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::diagMultiply
(
    Field<DiagType>& dDiag,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                upper[coeffI],
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::diagMultiplyCoeffT
(
    Field<DiagType>& dDiag,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                upper[coeffI].T(),        // Upper coefficient transposed
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::diagMultiply
(
    Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                lower[coeffI],
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::ILUmultiply
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    forAll (upper, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI], x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::ILUmultiplyCoeffT
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    forAll (upper, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI].T(), x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::ILUmultiply
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

    register label losortCoeff;

    forAll (lower, coeffI)
    {
        losortCoeff = losortAddr[coeffI];

        x[upperAddr[losortCoeff]] -=
            mult
            (
                dDiag[upperAddr[losortCoeff]],
                mult(lower[losortCoeff], x[lowerAddr[losortCoeff]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::myBlockCholeskyPrecon<Type>::ILUmultiplyTranspose
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

    register label losortCoeff;

    //HJ Not sure if the coefficient itself needs to be transposed.
    // HJ, 30/Oct/2007
    forAll (lower, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI], x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        losortCoeff = losortAddr[coeffI];

        x[lowerAddr[losortCoeff]] -=
            mult
            (
                dDiag[lowerAddr[losortCoeff]],
                mult(lower[losortCoeff], x[upperAddr[losortCoeff]])
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::myBlockCholeskyPrecon<Type>::myBlockCholeskyPrecon
(
    const BlockLduMatrix<Type>& matrix
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag())
{
    calcPreconDiag();
}


template<class Type>
Foam::myBlockCholeskyPrecon<Type>::myBlockCholeskyPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag())
{
    calcPreconDiag();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::myBlockCholeskyPrecon<Type>::~myBlockCholeskyPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::myBlockCholeskyPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
    }
}


template<class Type>
void Foam::myBlockCholeskyPrecon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        precondition(xT, bT);
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
    }
}


// ************************************************************************* //
