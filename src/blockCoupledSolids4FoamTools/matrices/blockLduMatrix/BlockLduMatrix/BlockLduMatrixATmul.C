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

Description
    Vector-matrix multiplication operations for a block matrix

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"
#include "vectorBlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockLduMatrix<Type>::Amul
(
    TypeField& Ax,
    const TypeField& x
) const
{
    Ax = pTraits<Type>::zero;

    // Initialise the update of coupled interfaces
    initInterfaces(coupleUpper_, Ax, x);

    AmulCore(Ax, x);

    // Update coupled interfaces
    updateInterfaces(coupleUpper_, Ax, x);
}


template<class Type>
void Foam::BlockLduMatrix<Type>::AmulCore
(
    TypeField& Ax,
    const TypeField& x
) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const unallocLabelList& u = lduAddr().upperAddr();
    const unallocLabelList& l = lduAddr().lowerAddr();

    const TypeCoeffField& Diag = this->diag();
    const TypeCoeffField& Upper = this->upper();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Diagonal multiplication, no indirection
    multiply(Ax, Diag, x);

    // Lower multiplication

    if (symmetric())
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Ax[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeUpper = Upper.asLinear();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Ax[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeUpper = Upper.asSquare();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                // Use transpose upper coefficient
                Ax[u[coeffI]] +=
                    mult(activeUpper[coeffI].T(), x[l[coeffI]]);
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& Lower = this->lower();

        if (Lower.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeLower = Lower.asScalar();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Ax[u[coeffI]] += mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
        else if (Lower.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeLower = Lower.asLinear();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Ax[u[coeffI]] += mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
        else if (Lower.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeLower = Lower.asSquare();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Ax[u[coeffI]] += mult(activeLower[coeffI], x[l[coeffI]]);
            }
        }
    }


    // Upper multiplication

    if (Upper.activeType() == blockCoeffBase::SCALAR)
    {
        const scalarTypeField& activeUpper = Upper.asScalar();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            Ax[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
        }
    }
    else if (Upper.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& activeUpper = Upper.asLinear();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            Ax[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
        }
    }
    else if (Upper.activeType() == blockCoeffBase::SQUARE)
    {
        const squareTypeField& activeUpper = Upper.asSquare();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            Ax[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
        }
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::Tmul
(
    TypeField& Ax,
    const TypeField& x
) const
{
    Ax = pTraits<Type>::zero;

    // Initialise the update of coupled interfaces
    initInterfaces(coupleLower_, Ax, x);

    TmulCore(Ax, x);

    // Update coupled interfaces
    updateInterfaces(coupleLower_, Ax, x);
}


template<class Type>
void Foam::BlockLduMatrix<Type>::TmulCore
(
    TypeField& Tx,
    const TypeField& x
) const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const unallocLabelList& u = lduAddr().upperAddr();
    const unallocLabelList& l = lduAddr().lowerAddr();

    const TypeCoeffField& Diag = this->diag();
    const TypeCoeffField& Upper = this->upper();

    // Diagonal multiplication, no indirection
    multiply(Tx, Diag, x);

    // Upper multiplication

    if (Upper.activeType() == blockCoeffBase::SCALAR)
    {
        const scalarTypeField& activeUpper = Upper.asScalar();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            Tx[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
        }
    }
    else if (Upper.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& activeUpper = Upper.asLinear();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            Tx[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
        }
    }
    else if (Upper.activeType() == blockCoeffBase::SQUARE)
    {
        const squareTypeField& activeUpper = Upper.asSquare();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            Tx[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
        }
    }

    // Lower multiplication

    if (symmetric())
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Tx[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeUpper = Upper.asLinear();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Tx[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeUpper = Upper.asSquare();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                // Use transpose upper coefficient
                Tx[l[coeffI]] +=
                    mult(activeUpper[coeffI].T(), x[u[coeffI]]);
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& Lower = this->lower();

        if (Lower.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeLower = Lower.asScalar();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Tx[l[coeffI]] += mult(activeLower[coeffI], x[u[coeffI]]);
            }
        }
        else if (Lower.activeType() == blockCoeffBase::LINEAR)
        {
            const linearTypeField& activeLower = Lower.asLinear();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Tx[l[coeffI]] += mult(activeLower[coeffI], x[u[coeffI]]);
            }
        }
        else if (Lower.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeLower = Lower.asSquare();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                Tx[l[coeffI]] += mult(activeLower[coeffI], x[u[coeffI]]);
            }
        }
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::segregateB
(
    TypeField& sMul,
    const TypeField& x
) const
{
    typedef typename TypeCoeffField::linearType linearType;
    typedef typename TypeCoeffField::squareType squareType;

    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const unallocLabelList& u = lduAddr().upperAddr();
    const unallocLabelList& l = lduAddr().lowerAddr();

    // Diagonal multiplication
    if (thereIsDiag())
    {
        if (diag().activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeDiag = this->diag().asSquare();
            linearTypeField lf(activeDiag.size());
            squareTypeField sf(activeDiag.size());

            // Expand and contract
            contractLinear(lf, activeDiag);
            expandLinear(sf, lf);

            sMul -= (activeDiag - sf) & x;
        }
    }

    // Lower multiplication

    if (thereIsLower())
    {
        if (lower().activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeLower = this->lower().asSquare();

            // Auxiliary variables used in expand/contract
            linearType lt;
            squareType st;

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                contractLinear(lt, activeLower[coeffI]);
                expandLinear(st, lt);

                sMul[u[coeffI]] -= (activeLower[coeffI] - st) & x[l[coeffI]];
            }
        }
    }

    // Upper multiplication

    if (thereIsUpper())
    {
        if (upper().activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& activeUpper = this->upper().asSquare();

            // Auxiliary variables used in expand/contract
            linearType lt;
            squareType st;

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                contractLinear(lt, activeUpper[coeffI]);
                expandLinear(st, lt);

                sMul[l[coeffI]] -= (activeUpper[coeffI] - st) & x[u[coeffI]];
            }

            // If the matrix is symmetric, the lower triangular product
            // is also needed
            if (symmetric())
            {
                for (register label coeffI = 0; coeffI < u.size(); coeffI++)
                {
                    // Use transpose upper coefficient
                    contractLinear(lt, activeUpper[coeffI]);
                    expandLinear(st, lt);
                    sMul[u[coeffI]] -=
                        (activeUpper[coeffI].T() - st) & x[l[coeffI]];
                }
            }
        }
    }
}


// ************************************************************************* //
