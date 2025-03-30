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

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::movingTractionFvPatchVectorField

Description
    Specifies contact traction for a traction boundary patch according to
    the Riemann solver. External boundary traction can also be specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        linearMomentum  | boundary momentum    | no          | vector(0,0,0)
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        TypeName            movingTraction;
        linearMomentum      (10,0 0);
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    movingTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingTractionFvPatchVectorField_H
#define movingTractionFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fvsPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class movingTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingTractionFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:

    // Boundary linear momentum
    vector lm_P_;

public:

    // Runtime type information
    TypeName("movingTraction");


    // Constructors

        //- Construct from patch and internal field
        movingTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingTractionFvPatchVectorField
        //  onto a new patch
        movingTractionFvPatchVectorField
        (
            const movingTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingTractionFvPatchVectorField
        (
            const movingTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingTractionFvPatchVectorField
        (
            const movingTractionFvPatchVectorField&,
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
                new movingTractionFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        // Access

            virtual vector& linearMomentum()
            {
                return lm_P_;
            }

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