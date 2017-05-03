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
    BlockMUMPS

Description
    A wrapper class to call the multi-frontal sparse direct solver from the
    MUMPS library; more details can be found at http://mumps.enseeiht.fr

    Notes
    - The solver does not run in parallel
    - The solver is currently only implemented for vectors

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    BlockMUMPS.C

\*---------------------------------------------------------------------------*/

#ifndef BlockMUMPS_H
#define BlockMUMPS_H

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
                      Class BlockMUMPS Declaration
\*---------------------------------------------------------------------------*/

class BlockMUMPS
:
    public blockVectorSolver //BlockLduSolver<vector>
{
        //- Write matrix and source to file for solution in Matlab
        // System can be read and solved like this:
        //     % Matlab '.m' script
        //     load matlabSparseMatrix.txt;
        //     A = spconvert(matrix);
        //     B = dlmread('matlabSource.txt', ' ');
        //     x = A\B;
        const Switch exportSystemMatlab_;

        // If the system should be scaled before solving to equilibrate rows
        // and column norms in matrices
        const Switch scale_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockMUMPS(const BlockMUMPS&);

        //- Disallow default bitwise assignment
        void operator=(const BlockMUMPS&);


public:

    //- Runtime type information
    TypeName("MUMPS");


    // Constructors

        //- Construct from matrix
        BlockMUMPS
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
//#   include "BlockMUMPS.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
