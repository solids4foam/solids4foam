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
    Foam::fixedValueCorrectedFvPatchVectorField

Description
    Fixed value boundary condition with non-orthogonal corrections.

Author
    Philip Cardiff, UCD.

SourceFiles
    fixedValueCorrectedFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueCorrectedFvPatchVectorField_H
#define fixedValueCorrectedFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedValueCorrectedFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedValueCorrectedFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

        //- Enable non-orthogonal corrections
        const Switch nonOrthogonalCorrections_;

public:

    //- Runtime type information
    TypeName("fixedValueCorrected");


    // Constructors

        //- Construct from patch and internal field
        fixedValueCorrectedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueCorrectedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueCorrectedFvPatchVectorField
        //  onto a new patch
        fixedValueCorrectedFvPatchVectorField
        (
            const fixedValueCorrectedFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueCorrectedFvPatchVectorField
        (
            const fixedValueCorrectedFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedValueCorrectedFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        fixedValueCorrectedFvPatchVectorField
        (
            const fixedValueCorrectedFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedValueCorrectedFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
