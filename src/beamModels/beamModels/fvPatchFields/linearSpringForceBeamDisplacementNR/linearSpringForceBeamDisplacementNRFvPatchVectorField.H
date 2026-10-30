/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    linearSpringForceBeamDisplacementNRFvPatchVectorField

Description
    A linear spring boundary condition for the beam solver which takes the
    spring stiffness and unit vector direction of spring to the
    boundary patch
    -- Derived from forceBeamDisplacementNR boundary condition

SourceFiles
    linearSpringForceBeamDisplacementNRFvPatchVectorField.C

Authors
    Seevani Bali, UCD.
\*---------------------------------------------------------------------------*/

#ifndef linearSpringForceBeamDisplacementNRFvPatchVectorField_H
#define linearSpringForceBeamDisplacementNRFvPatchVectorField_H

#include "fvPatchFields.H"
#include "forceBeamDisplacementNRFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class linearSpringForceBeamDisplacementNRFvPatch Declaration
\*---------------------------------------------------------------------------*/

class linearSpringForceBeamDisplacementNRFvPatchVectorField
:
    public forceBeamDisplacementNRFvPatchVectorField
{
    // Private Data

    // Spring Stiffness
    scalar springStiffness_;

    // Spring Direction
    vector springDirection_;

    // Spring Stiffness tensor
    tensorField Ks_;

    //- Offset force at spring BC end of beam
    vectorField offsetForce_;

    //- Offset Force time series if present
    interpolationTable<vector> offsetForceSeries_;


public:

    //- Runtime type information
    TypeName("linearSpringForceBeamDisplacementNR");


    // Constructors

        //- Construct from patch and internal field
        linearSpringForceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        linearSpringForceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  linearSpringForceBeamDisplacementNRFvPatchVectorField onto a new patch
        linearSpringForceBeamDisplacementNRFvPatchVectorField
        (
            const linearSpringForceBeamDisplacementNRFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        linearSpringForceBeamDisplacementNRFvPatchVectorField
        (
            const linearSpringForceBeamDisplacementNRFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new linearSpringForceBeamDisplacementNRFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        linearSpringForceBeamDisplacementNRFvPatchVectorField
        (
            const linearSpringForceBeamDisplacementNRFvPatchVectorField&,
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
                new linearSpringForceBeamDisplacementNRFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets Updated to false
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
