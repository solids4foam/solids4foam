/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    myBlockCholeskyPrecon

Description
    Incomplete Cholesky preconditioning with no fill-in.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    myBlockCholeskyPrecon.C
    myBlockCholeskyPreconDecoupled.C

\*---------------------------------------------------------------------------*/

#ifndef myBlockCholeskyPrecon_H
#define myBlockCholeskyPrecon_H

#include "BlockLduPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class myBlockCholeskyPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class myBlockCholeskyPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Data

        //- Preconditioned diagonal
        mutable CoeffField<Type> preconDiag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        myBlockCholeskyPrecon(const myBlockCholeskyPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const myBlockCholeskyPrecon&);

        //- Precondition the diagonal
        void calcPreconDiag();

        // Diagonal multiplication, symmetric matrix
        template<class DiagType, class ULType>
        void diagMultiply
        (
            Field<DiagType>& dDiag,
            const Field<ULType>& upper
        );

        //- Diagonal multiplication with transpose upper square coeff
        //  for the symmetric matrix
        template<class DiagType, class ULType>
        void diagMultiplyCoeffT
        (
            Field<DiagType>& dDiag,
            const Field<ULType>& upper
        );

        //- Diagonal multiplication, asymmetric matrix
        template<class DiagType, class ULType>
        void diagMultiply
        (
            Field<DiagType>& dDiag,
            const Field<ULType>& lower,
            const Field<ULType>& upper
        );

        //- ILU multiplication, symmetric matrix
        template<class DiagType, class ULType>
        void ILUmultiply
        (
            Field<Type>& x,
            const Field<DiagType>& dDiag,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;

        //- ILU multiplication, with transpose upper square coeff
        //  for a symmetric matrix
        template<class DiagType, class ULType>
        void ILUmultiplyCoeffT
        (
            Field<Type>& x,
            const Field<DiagType>& dDiag,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;

        //- ILU multiplication, asymmetric matrix
        template<class DiagType, class ULType>
        void ILUmultiply
        (
            Field<Type>& x,
            const Field<DiagType>& dDiag,
            const Field<ULType>& lower,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;

        //- ILU multiplication transposed asymmetric matrix
        template<class DiagType, class ULType>
        void ILUmultiplyTranspose
        (
            Field<Type>& x,
            const Field<DiagType>& dDiag,
            const Field<ULType>& lower,
            const Field<ULType>& upper,
            const Field<Type>& b
        ) const;


        // Decoupled operations, used in template specialisation

            //- Precondition the diagonal, decoupled version
            void calcDecoupledPreconDiag();

            //- Execute preconditioning, decoupled version
            void decoupledPrecondition
            (
                Field<Type>& x,
                const Field<Type>& b
            ) const;

            //- Execute preconditioning with matrix transpose,
            //  decoupled version
            void decoupledPreconditionT
            (
                Field<Type>& xT,
                const Field<Type>& bT
            ) const;


public:

    //- Runtime type information
    TypeName("myCholesky");


    // Constructors

        //- Construct from matrix for smoother use
        myBlockCholeskyPrecon
        (
            const BlockLduMatrix<Type>& matrix
        );

        //- Construct from components
        myBlockCholeskyPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~myBlockCholeskyPrecon();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "myBlockCholeskyPrecon.C"
#   include "myBlockCholeskyPreconDecoupled.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
