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
    BlockTrilinosSolver

Description
    A wrapper to call linear solvers from the Trilinos package.

    The foam matix is first converted to the cmpressed-row-format and then
    passed to the Trilinos AztecOO linear solver, currently set to GMRes.

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    BlockTrilinosSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockTrilinosSolver_H
#define BlockTrilinosSolver_H

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

/*---------------------------------------------------------------------------*\
                      Class BlockTrilinosSolver Declaration
\*---------------------------------------------------------------------------*/

class BlockTrilinosSolver
:
    public blockVectorSolver //BlockLduSolver<vector>
{
        // Tolerance for iterative solution
        const scalar tol_;

        // Number of directions for GMRes
        const label nDir_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockTrilinosSolver(const BlockTrilinosSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockTrilinosSolver&);


public:

    //- Runtime type information
    TypeName("Trilinos");


    // Constructors

        //- Construct from matrix
        BlockTrilinosSolver
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
//#   include "BlockTrilinosSolver.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
