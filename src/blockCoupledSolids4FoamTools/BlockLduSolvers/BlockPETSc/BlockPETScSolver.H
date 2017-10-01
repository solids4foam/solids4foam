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
    BlockPETScSolver

Description
    A wrapper class to call linear solvers from the PETSc library.

    WIP

    Notes
    - The solver does not run in parallel
    - The solver is currently only implemented for vectors

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    BlockPETScSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockPETScSolver_H
#define BlockPETScSolver_H

#include "Field.H"
#include "typeInfo.H"

#include "runTimeSelectionTables.H"

#include "Switch.H"

#include "BlockLduMatrix.H"
#include "BlockLduSolver.H"
#include "blockLduSolvers.H"
#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// template<class Type>
// class BlockLduMatrix;

// template<class Type>
// class BlockLduSolver;

// template<class vector>
// class BlockSolverPerformance;

/*---------------------------------------------------------------------------*\
                      Class BlockPETScSolver Declaration
\*---------------------------------------------------------------------------*/

class BlockPETScSolver
:
    public blockVectorSolver
{
        //- PETSc command line arguments for choosing the solver and settings
        const wordList PETScArgs_;

        //- Tolerance for iterative solvers
        const scalar tol_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockPETScSolver(const BlockPETScSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockPETScSolver&);


public:

    //- Runtime type information
    TypeName("PETSc");


    // Constructors

        //- Construct from matrix
        BlockPETScSolver
        (
            const word& fieldName,
            const BlockLduMatrix<vector>& matrix,
            const dictionary& dict
        );


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<vector> solve
        (
            Field<vector>& x,
            const Field<vector>& b
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "BlockPETScSolver.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
