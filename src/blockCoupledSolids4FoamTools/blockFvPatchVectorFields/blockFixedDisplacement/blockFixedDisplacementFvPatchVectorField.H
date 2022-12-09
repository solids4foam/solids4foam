/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::blockFixedDisplacementFvPatchVectorField

Description
    BlockFixedDisplacement boundary condition for coupled solid solver.
    This is equivalent to a fixedValue boundary condition.

    Implements
        U = I & Ub

    where Ub is the user specified blockFixedDisplacement

SourceFiles
    blockFixedDisplacementFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockFixedDisplacementFvPatchVectorField_H
#define blockFixedDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "blockFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blockFixedDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class blockFixedDisplacementFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public blockFvPatchVectorField
{
    // Private Data

        //- Constant total displacement
        vector totalDisp_;

        //- Displacement time series
        interpolationTable<vector> dispSeries_;

public:

    //- Runtime type information
    TypeName("blockFixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        blockFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blockFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blockFixedDisplacementFvPatchVectorField
        //  onto a new patch
        blockFixedDisplacementFvPatchVectorField
        (
            const blockFixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blockFixedDisplacementFvPatchVectorField
        (
            const blockFixedDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blockFixedDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blockFixedDisplacementFvPatchVectorField
        (
            const blockFixedDisplacementFvPatchVectorField&,
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
                new blockFixedDisplacementFvPatchVectorField(*this, iF)
            );
        }



    // Member functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

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
