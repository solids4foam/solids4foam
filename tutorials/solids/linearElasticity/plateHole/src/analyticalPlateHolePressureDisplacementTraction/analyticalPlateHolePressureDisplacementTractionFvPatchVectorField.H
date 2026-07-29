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
    analyticalPlateHolePressureDisplacementTractionFvPatchVectorField

Description
    Pressure-displacement traction boundary condition for the plate-hole test
    case to enforce analytical traction at the boundaries.

SourceFiles
    analyticalPlateHolePressureDisplacementTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalPlateHolePressureDisplacementTractionFvPatchVectorField_H
#define analyticalPlateHolePressureDisplacementTractionFvPatchVectorField_H

#include "tractionPressureDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class analyticalPlateHolePressureDisplacementTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
:
    public tractionPressureDisplacementFvPatchVectorField
{
    // Private Data

        //- Far field traction
        //  It is assumed that it acts in the positive X direction
        const scalar T_;

        //- Hole radius
        const scalar holeR_;


public:

    //- Runtime type information
    TypeName("analyticalPlateHolePressureDisplacementTraction");


    // Constructors

        //- Construct from patch and internal field
        analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
        //  onto a new patch
        analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
        (
            const analyticalPlateHolePressureDisplacementTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
        (
            const analyticalPlateHolePressureDisplacementTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
        (
            const analyticalPlateHolePressureDisplacementTractionFvPatchVectorField&,
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
                new analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
