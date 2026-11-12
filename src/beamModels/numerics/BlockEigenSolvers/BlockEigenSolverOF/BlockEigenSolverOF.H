/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    BlockEigenSolverOF

Description
    Wrapper class to solve a 6x6 block coupled linear system using the Eigen
    header library.

    This class compiles with OpenFOAM-v2012.

Author
    Philip Cardiff, UCD
    Amirhossein Taran, UCD
    Based on class from Zeljko Tukovic, FMENA, Zagreb

SourceFiles
    BlockEigenSolverOF.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenSolverOF_H
#define BlockEigenSolverOF_H

// #include <functional>
// #include <numeric>
#include "typeInfo.H"
#include "Switch.H"
#include <Eigen/Sparse>
#include "scalarMatrices.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockEigenSolverOF Declaration
\*---------------------------------------------------------------------------*/

class BlockEigenSolverOF
{
    // Private Data

        //- Const-reference to the sparse matrix diagonal
        const Field<scalarSquareMatrix>& d_;

        //- Const-reference to the sparse matrix upper
        const Field<scalarSquareMatrix>& l_;

        //- Const-reference to the sparse matrix lower
        const Field<scalarSquareMatrix>& u_;

        //- Const-reference to owner list
        const labelList& own_;

        //- Const-reference to neighbour list
        const labelList& nei_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenSolverOF(const BlockEigenSolverOF&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenSolverOF&);

protected:

        //- Convert foam matrix to the Eigen matrix format
        void convertFoamMatrixToEigenMatrix
        (
            const Field<scalarSquareMatrix>& d,
            const Field<scalarSquareMatrix>& l,
            const Field<scalarSquareMatrix>& u,
            const labelList& own,
            const labelList& nei,
            Eigen::SparseMatrix<scalar>& A
        );

public:

    //- Runtime type information
    TypeName("Eigen");

    // Constructors

        //- Construct from matrix
        BlockEigenSolverOF
        (
            const Field<scalarSquareMatrix>& d,
            const Field<scalarSquareMatrix>& l,
            const Field<scalarSquareMatrix>& u,
            const labelList& own,
            const labelList& nei
            // const dictionary& dict
        );

    // Destructor

        virtual ~BlockEigenSolverOF()
        {}

    // Member Functions

        //- Solve the matrix with this source and solution vector and return the
        //  initial residual magnitude
        scalar solve
        (
            Field<scalarRectangularMatrix>& x,
            const Field<scalarRectangularMatrix>& b
        );

        // //- Re-initialise solver after matrix coefficient update
        // virtual void initMatrix()
        // {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
