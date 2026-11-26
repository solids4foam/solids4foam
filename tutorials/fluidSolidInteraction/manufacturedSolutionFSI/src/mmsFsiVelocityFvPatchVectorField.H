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
    Foam::mmsFsiVelocityFvPatchVectorField

Description
    Fixed value fluid velocity of the manufactured FSI solution, evaluated at
    the current face centres and time.

Usage
    \verbatim
    inlet
    {
        type        mmsFsiVelocity;
        value       uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    mmsFsiVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef mmsFsiVelocityFvPatchVectorField_H
#define mmsFsiVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class mmsFsiVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class mmsFsiVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("mmsFsiVelocity");


    // Constructors

        //- Construct from patch and internal field
        mmsFsiVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mmsFsiVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mmsFsiVelocityFvPatchVectorField onto
        //  a new patch
        mmsFsiVelocityFvPatchVectorField
        (
            const mmsFsiVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mmsFsiVelocityFvPatchVectorField
        (
            const mmsFsiVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new mmsFsiVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mmsFsiVelocityFvPatchVectorField
        (
            const mmsFsiVelocityFvPatchVectorField&,
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
                new mmsFsiVelocityFvPatchVectorField(*this, iF)
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
