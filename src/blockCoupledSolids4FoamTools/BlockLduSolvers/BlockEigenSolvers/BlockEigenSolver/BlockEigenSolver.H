/*---------------------------------------------------------------------------*\
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

Class
    BlockEigenSolver

Description
    An abstract base class wrapper for linear solvers from the Eigen
    (http://eigen.tuxfamily.org) library.

    Notes
    - The solver is currently only implemented for vectors
    - I would recommend removing '-Wold-style-cast' from
      $WM_PROJECT_DIR/wmake/rules/$WM_ARCH/c++ as the terminal will be polluted
      with 'old-stlye-casts' warnings from the compilation of the Eigen header
      files.

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    BlockEigenSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenSolver_H
#define BlockEigenSolver_H

#include "typeInfo.H"
#include "Switch.H"
#include <Eigen/Sparse>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class BlockEigenSolver Declaration
\*---------------------------------------------------------------------------*/

class BlockEigenSolver
:
    public blockVectorSolver
{
        //- Write matrix and source to file for solution in Matlab
        const Switch writeMatlabFiles_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenSolver(const BlockEigenSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenSolver&);


protected:

        //- Check if the matrix is for a 2-D case
        bool checkTwoD() const;

        //- Calculate the number of degrees of freedom
        int calcDegreesOfFreedom
        (
            const BlockLduMatrix<vector>& matrix,
            const bool twoD
        ) const;

        //- Convert foam matrix to the Eigen matrix format
        void convertFoamMatrixToEigenMatrix
        (
            const BlockLduMatrix<vector>& matrix,
            Eigen::SparseMatrix<scalar>& A
        ) const;

        //- Switch to enable exporting of the linear system to the Matlab format
        const Switch writeMatlabFiles() const
        {
            return writeMatlabFiles_;
        }

        //- Export of the linear system to the Matlab format
        void writeLinearSystemToMatlabFiles
        (
            const Eigen::SparseMatrix<scalar>& A,
            const Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b
        ) const;

public:

    //- Runtime type information
    TypeName("EigenSolver");


    // Constructors

        //- Construct from matrix
        BlockEigenSolver
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
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
