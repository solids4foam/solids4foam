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
    BlockEigenGMResSolver

Description
    The GMRes iterative sparse linear equation solver from package Eigen is
    called to solve the system.

    Notes
    - The solver does not run in parallel
    - The solver is currently only implemented for vectors
    - I would recommend removing '-Wold-style-cast' from
      $WM_PROJECT_DIR/wmake/rules/$WM_ARCH/c++ as the terminal will be polluted
      with 'old-stlye-casts' warnings from the compilation of the Eigen header
      files.

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    BlockEigenGMResSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenGMResSolver_H
#define BlockEigenGMResSolver_H

#include "BlockEigenSolver.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockEigenGMResSolver Declaration
\*---------------------------------------------------------------------------*/

class BlockEigenGMResSolver
:
    public BlockEigenSolver
{
        // Tolerance for iterative solution
        const scalar tol_;

        // Restart frequency
        const int restart_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenGMResSolver(const BlockEigenGMResSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenGMResSolver&);


public:

    //- Runtime type information
    TypeName("EigenGMRes");


    // Constructors

        //- Construct from matrix
        BlockEigenGMResSolver
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
  
        //- Re-initialise solver after matrix coefficient update
        virtual void initMatrix()
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
