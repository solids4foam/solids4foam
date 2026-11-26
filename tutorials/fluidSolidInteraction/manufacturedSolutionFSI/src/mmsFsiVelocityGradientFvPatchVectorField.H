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
    Foam::mmsFsiVelocityGradientFvPatchVectorField

Description
    Fixed gradient fluid velocity of the manufactured FSI solution: the
    normal gradient is n & grad(U) of the exact velocity, evaluated at the
    current face centres and time.

Usage
    \verbatim
    outlet
    {
        type        mmsFsiVelocityGradient;
        value       uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    mmsFsiVelocityGradientFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef mmsFsiVelocityGradientFvPatchVectorField_H
#define mmsFsiVelocityGradientFvPatchVectorField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class mmsFsiVelocityGradientFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class mmsFsiVelocityGradientFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("mmsFsiVelocityGradient");


    // Constructors

        //- Construct from patch and internal field
        mmsFsiVelocityGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mmsFsiVelocityGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mmsFsiVelocityGradientFvPatchVectorField onto a new patch
        mmsFsiVelocityGradientFvPatchVectorField
        (
            const mmsFsiVelocityGradientFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mmsFsiVelocityGradientFvPatchVectorField
        (
            const mmsFsiVelocityGradientFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new mmsFsiVelocityGradientFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mmsFsiVelocityGradientFvPatchVectorField
        (
            const mmsFsiVelocityGradientFvPatchVectorField&,
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
                new mmsFsiVelocityGradientFvPatchVectorField(*this, iF)
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
