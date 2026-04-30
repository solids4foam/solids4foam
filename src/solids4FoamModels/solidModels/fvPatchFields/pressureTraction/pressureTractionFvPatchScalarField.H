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
    pressureTractionFvPatchScalarField

Description
    Pressure-increment boundary condition paired with
    tractionPressureDisplacement in coupledPressureDisplacementSolid cases.

SourceFiles
    pressureTractionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureTractionFvPatchScalarField_H
#define pressureTractionFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pressureTractionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class pressureTractionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

public:

    //- Runtime type information
    TypeName("pressureTraction");

    // Constructors

        //- Construct from patch and internal field
        pressureTractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureTractionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureTractionFvPatchScalarField onto a new patch
        pressureTractionFvPatchScalarField
        (
            const pressureTractionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureTractionFvPatchScalarField
        (
            const pressureTractionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureTractionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureTractionFvPatchScalarField
        (
            const pressureTractionFvPatchScalarField&,
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
                new pressureTractionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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


            //- Return patch-normal gradient
            virtual tmp<Field<scalar> > snGrad() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

        //- Normal traction
        tmp<scalarField> normalTraction() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
