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
    fixedGradientCorrectedFvPatchScalarField

Description
    Fixed gradient boundary condition with non-orthogonal corrections.

Author
    Philip Cardiff, UCD.

SourceFiles
    fixedGradientCorrectedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedGradientCorrectedFvPatchScalarField_H
#define fixedGradientCorrectedFvPatchScalarField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedGradientCorrectedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedGradientCorrectedFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

    // Private Data

        //- Enable non-orthogonal corrections
        const Switch nonOrthogonalCorrections_;

        //- Second order accuracy correction
        const Switch secondOrder_;

        //- Extrapolate values to the patch using the patch internal gradient
        //  field. Alternatively, if false, use the gradient at the patch
        Switch extrapolateValue_;


public:

    //- Runtime type information
    TypeName("fixedGradientCorrected");


    // Constructors

        //- Construct from patch and internal field
        fixedGradientCorrectedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedGradientCorrectedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedGradientCorrectedFvPatchScalarField onto a new patch
        fixedGradientCorrectedFvPatchScalarField
        (
            const fixedGradientCorrectedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedGradientCorrectedFvPatchScalarField
        (
            const fixedGradientCorrectedFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedGradientCorrectedFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        fixedGradientCorrectedFvPatchScalarField
        (
            const fixedGradientCorrectedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedGradientCorrectedFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        // Access

            // Non-const access to the extrapolateValue flag
            Switch& extrapolateValue()
            {
                return extrapolateValue_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
#ifdef OPENFOAM_NOT_EXTEND
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
