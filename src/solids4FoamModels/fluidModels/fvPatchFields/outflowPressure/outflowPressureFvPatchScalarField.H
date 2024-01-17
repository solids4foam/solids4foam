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
    Foam::outflowPressureFvPatchScalarField

Description
    Foam::outflowPressureFvPatchScalarField

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    outflowPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef outflowPressureFvPatchScalarField_H
#define outflowPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class outflowPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class outflowPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

//     //To keep track when to copy patch face centres
//     label myTimeIndex_;

//     //Old face centres of patch. Needed for (second order) time integration
//     scalarField Fc_;
//     scalarField oldFc_;
//     scalarField oldoldFc_;


public:

    //- Runtime type information
    TypeName("outflowPressure");


    // Constructors

        //- Construct from patch and internal field
        outflowPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outflowPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outflowPressureFvPatchScalarField
        //  onto a new patch
        outflowPressureFvPatchScalarField
        (
            const outflowPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outflowPressureFvPatchScalarField
        (
            const outflowPressureFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new outflowPressureFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        outflowPressureFvPatchScalarField
        (
            const outflowPressureFvPatchScalarField&,
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
                new outflowPressureFvPatchScalarField(*this, iF)
            );
        }



    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

    // Evaluation functions

//         //- Return patch-normal gradient
//         virtual tmp<Field<scalar> > snGrad() const;

//         //- Return the matrix source coefficients corresponding to the
//         //  evaluation of the gradient of this patchField
//         virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
