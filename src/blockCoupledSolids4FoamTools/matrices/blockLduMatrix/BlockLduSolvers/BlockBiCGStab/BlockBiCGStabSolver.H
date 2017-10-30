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
    BlockBiCGStabSolver

Description
    Preconditioned Bi-Conjugate Gradient stabilised solver.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockBiCGStabSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockBiCGStabSolver_H
#define BlockBiCGStabSolver_H

#include "blockLduSolvers.H"
#include "BlockIterativeSolver.H"
#include "blockLduPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BlockBiCGStabSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockBiCGStabSolver
:
    public BlockIterativeSolver<Type>
{
    // Private data

        //- Preconditioner
        autoPtr<BlockLduPrecon<Type> > preconPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockBiCGStabSolver(const BlockBiCGStabSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockBiCGStabSolver<Type>&);


public:

    //- Runtime type information
    TypeName("BiCGStab");


    // Constructors

        //- Construct from matrix components and solver data stream
        BlockBiCGStabSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


        // Destructor

            virtual ~BlockBiCGStabSolver()
            {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<Type> solve
        (
            Field<Type>& x,
            const Field<Type>& b
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockBiCGStabSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
