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
    blockFvPatchVectorField

Description
    Abstract base class for boundary conditions for the solid block coupled
    finite volume solver.
    Defines virtual insertBlockCoeffs function to insert boundary condition
    equation into the block coupled linear system.

Author
   Philip Cardiff UCD

SourceFiles
    blockFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef blockFvPatchVectorField_H
#define blockFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidPolyMesh.H"
#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class blockFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class blockFvPatchVectorField
{

public:


    // Constructors

    //- Runtime type information
    TypeName("blockFvPatchVectorField");

        //- Null Constructor
        blockFvPatchVectorField();

    // Destructors

        virtual ~blockFvPatchVectorField()
        {}


    // Member functions

        //- Insert coefficients into the block matrix
        virtual void insertBlockCoeffs
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& muf,
            const surfaceScalarField& lambdaf,
            const GeometricField<vector, fvPatchField, volMesh>& U,
            Field<vector>& blockB,
            BlockLduMatrix<vector>& blockM
        ) const = 0;

        //- Set snGrad on the patch after solution
        //virtual void setSnGrad() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
