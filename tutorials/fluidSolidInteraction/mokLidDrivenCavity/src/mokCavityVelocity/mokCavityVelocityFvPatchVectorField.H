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
    mokCavityVelocityFvPatchVectorField

Description
    Velocity of the lid and the inlet of the Mok lid-driven cavity: a velocity
    in the x direction with the periodic magnitude

        vBar(t) = 1 - cos(2 pi t/period),

    optionally scaled by a profile that rises linearly in y from zero at
    yBottom to one at yTop:

        u = vBar(t) (y - yBottom)/(yTop - yBottom).

    It replaces codedFixedValue conditions, which foam-extend does not
    provide.

Usage
    \verbatim
    lid
    {
        type        mokCavityVelocity;
        period      5;
        value       uniform (0 0 0);
    }
    inlet
    {
        type        mokCavityVelocity;
        period      5;
        yBottom     0.875;
        yTop        1;
        value       uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    mokCavityVelocityFvPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mokCavityVelocityFvPatchVectorField_H
#define mokCavityVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mokCavityVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class mokCavityVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Period of the velocity magnitude
        scalar period_;

        //- Whether the velocity is scaled by the linear profile in y
        bool linearProfile_;

        //- Height at which the linear profile is zero
        scalar yBottom_;

        //- Height at which the linear profile is one
        scalar yTop_;


public:

    //- Runtime type information
    TypeName("mokCavityVelocity");


    // Constructors

        //- Construct from patch and internal field
        mokCavityVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mokCavityVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mokCavityVelocityFvPatchVectorField
        //  onto a new patch
        mokCavityVelocityFvPatchVectorField
        (
            const mokCavityVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mokCavityVelocityFvPatchVectorField
        (
            const mokCavityVelocityFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new mokCavityVelocityFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        mokCavityVelocityFvPatchVectorField
        (
            const mokCavityVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new mokCavityVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
