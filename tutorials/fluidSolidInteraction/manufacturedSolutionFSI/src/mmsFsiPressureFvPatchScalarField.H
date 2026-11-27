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
    Foam::mmsFsiPressureFvPatchScalarField

Description
    Fixed value kinematic fluid pressure of the manufactured FSI solution,
    evaluated at the current face centres and time.

Usage
    \verbatim
    outlet
    {
        type        mmsFsiPressure;
        value       uniform 0;
    }
    \endverbatim

SourceFiles
    mmsFsiPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mmsFsiPressureFvPatchScalarField_H
#define mmsFsiPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class mmsFsiPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mmsFsiPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("mmsFsiPressure");


    // Constructors

        //- Construct from patch and internal field
        mmsFsiPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mmsFsiPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mmsFsiPressureFvPatchScalarField onto
        //  a new patch
        mmsFsiPressureFvPatchScalarField
        (
            const mmsFsiPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mmsFsiPressureFvPatchScalarField
        (
            const mmsFsiPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mmsFsiPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mmsFsiPressureFvPatchScalarField
        (
            const mmsFsiPressureFvPatchScalarField&,
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
                new mmsFsiPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
