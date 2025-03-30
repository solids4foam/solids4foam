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
    Foam::symmetricLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a symmetric boundary patch
    according to the Riemann solver. External boundary traction can
    be specified together with a linear ramp function.

Usage
    \table
        Property     | Description              | Required    | Default value
        traction     | boundary traction        | no          | vector(0,0,0)
        rampEndTime  | ramp time for traction   | no          | VSMALL
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            symmetricLinearMomentum;
        traction        (10 10 10);
        rampEndTime     1.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    symmetricLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetricLinearMomentumFvPatchVectorField_H
#define symmetricLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class symmetricLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class symmetricLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

private:

    // Boundary traction
    vector t_P_;

    // Traction ramp end time
    scalar tRamp_;


public:

    // Runtime type information
    TypeName("symmetricLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        symmetricLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetricLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  symmetricLinearMomentumFvPatchVectorField onto a new patch
        symmetricLinearMomentumFvPatchVectorField
        (
            const symmetricLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        symmetricLinearMomentumFvPatchVectorField
        (
            const symmetricLinearMomentumFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new symmetricLinearMomentumFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        symmetricLinearMomentumFvPatchVectorField
        (
            const symmetricLinearMomentumFvPatchVectorField&,
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
                new symmetricLinearMomentumFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        // Access

        // Mapping functions

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //