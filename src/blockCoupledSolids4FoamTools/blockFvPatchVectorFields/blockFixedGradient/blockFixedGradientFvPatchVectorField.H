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
    Foam::blockFixedGradientFvPatchVectorField

Description
    Fixed gradient boundary condition for coupled solid solver.

    Implements
        gb = n & gradU

    where gb is the user specified normal gradient.


SourceFiles
    blockFixedGradientFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockFixedGradientFvPatchVectorField_H
#define blockFixedGradientFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "blockFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blockFixedGradientFvPatch Declaration
\*---------------------------------------------------------------------------*/

class blockFixedGradientFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public blockFvPatchVectorField
{
    // Private Data

        //- Gradient
        vectorField gradient_;

public:

    //- Runtime type information
    TypeName("blockFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        blockFixedGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blockFixedGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blockFixedGradientFvPatchVectorField
        //  onto a new patch
        blockFixedGradientFvPatchVectorField
        (
            const blockFixedGradientFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blockFixedGradientFvPatchVectorField
        (
            const blockFixedGradientFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blockFixedGradientFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blockFixedGradientFvPatchVectorField
        (
            const blockFixedGradientFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new blockFixedGradientFvPatchVectorField(*this, iF)
            );
        }



    // Member functions


        // Access

            virtual const vectorField& gradient() const
            {
                return gradient_;
            }

            virtual vectorField& gradient()
            {
                return gradient_;
            }

        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const
            {
                return gradient_;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        // Block functions

            //- Insert coefficients into the block matrix
            virtual void insertBlockCoeffs
            (
                const solidPolyMesh& solidMesh,
                const surfaceScalarField& muf,
                const surfaceScalarField& lambdaf,
                const GeometricField<vector, fvPatchField, volMesh>& U,
                Field<vector>& blockB,
                BlockLduMatrix<vector>& blockM
            ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
