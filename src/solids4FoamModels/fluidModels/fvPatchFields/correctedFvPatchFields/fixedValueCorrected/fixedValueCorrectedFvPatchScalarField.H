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
    Foam::fixedValueCorrectedFvPatchScalarField

Description
    Fixed value boundary condition with non-orthogonal corrections.

Author
    Philip Cardiff, UCD.

SourceFiles
    fixedValueCorrectedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueCorrectedFvPatchScalarField_H
#define fixedValueCorrectedFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedValueCorrectedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedValueCorrectedFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

        //- Enable non-orthogonal corrections
        const Switch nonOrthogonalCorrections_;

public:

    //- Runtime type information
    TypeName("fixedValueCorrected");


    // Constructors

        //- Construct from patch and internal field
        fixedValueCorrectedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueCorrectedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueCorrectedFvPatchScalarField
        //  onto a new patch
        fixedValueCorrectedFvPatchScalarField
        (
            const fixedValueCorrectedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedValueCorrectedFvPatchScalarField
        (
            const fixedValueCorrectedFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new fixedValueCorrectedFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        fixedValueCorrectedFvPatchScalarField
        (
            const fixedValueCorrectedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new fixedValueCorrectedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<scalar> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
