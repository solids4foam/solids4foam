/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    BlockEigenILUTPrecon

Description
    Incomplete LUT preconditioning as implemented in the Eigen library.
    This class is just a wrapper for the Eigen library function.

    The foam matrix is copied into the Eigen vectors.

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    BlockEigenILUTPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenILUTPrecon_H
#define BlockEigenILUTPrecon_H

#include "BlockLduPrecon.H"
#include <Eigen/Sparse>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockEigenILUTPrecon Declaration
\*---------------------------------------------------------------------------*/

class BlockEigenILUTPrecon
:
    public BlockLduPrecon<vector>
{
    // Private Data

        //- Preconditioned diagonal
        //mutable CoeffField<vector> preconDiag_;

        //- Spare matrix stored in the Eigen format
        mutable autoPtr< Eigen::SparseMatrix<scalar> > APtr_;

        //- Eigne incomplete LU preconditoner with dial thresholding
        mutable autoPtr< Eigen::IncompleteLUT<scalar> > preconPtr_;

        //- Is the matrix for a 2-D problem
        bool twoD_;

        //- Fill-in factor (default is 10)
        const label fillFactor_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenILUTPrecon(const BlockEigenILUTPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenILUTPrecon&);

        //- Calculate the preconditioner
        void calcPrecon(const BlockLduMatrix<vector>& matrix);

        //- Return the matrix
        const Eigen::SparseMatrix<scalar>& A() const;

        //- Return the preconditioner
        const Eigen::IncompleteLUT<scalar>& precon() const;


public:

    //- Runtime type information
    TypeName("EigenILUT");


    // Constructors

        //- Construct from matrix for smoother use
        BlockEigenILUTPrecon
        (
            const BlockLduMatrix<vector>& matrix
        );

        //- Construct from components
        BlockEigenILUTPrecon
        (
            const BlockLduMatrix<vector>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockEigenILUTPrecon();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<vector>& x,
            const Field<vector>& b
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            Field<vector>& xT,
            const Field<vector>& bT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "BlockEigenILUTPrecon.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
