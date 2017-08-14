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
    BlockEigenProblemSolver

Description
    The SelfAdjointEigenSolver eigen problem solver from the package Eigen is
    called to calculate the eigenvalues and eigenvectors for the system.

    This solver is used to solve, for example, linear buckling problems, free
    vibration problems, etc.

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
    BlockEigenProblemSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenProblemSolver_H
#define BlockEigenProblemSolver_H

#include "BlockEigenSolver.H"

#include "Field.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"
#include "BlockLduMatrix.H"
#include "BlockLduSolver.H"
#include "blockLduSolvers.H"
#include "BlockSolverPerformance.H"
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class BlockEigenProblemSolver Declaration
\*---------------------------------------------------------------------------*/

class BlockEigenProblemSolver
:
    public BlockEigenSolver
{
        // Pointer to the eigenvalues
        autoPtr
        <
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf>::RealVectorType
        > eigenvaluesPtr_;

        // Pointer to the eigenvectors
        autoPtr
        <
            Eigen::MatrixXf
        > eigenvectorsPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenProblemSolver(const BlockEigenProblemSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenProblemSolver&);


public:

    //- Runtime type information
    TypeName("EigenProblem");


    // Constructors

        //- Construct from matrix
        BlockEigenProblemSolver
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

    // Access Functions

        //- Retrieve a mode
        void mode(const int modeI, Field<vector>& x) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
