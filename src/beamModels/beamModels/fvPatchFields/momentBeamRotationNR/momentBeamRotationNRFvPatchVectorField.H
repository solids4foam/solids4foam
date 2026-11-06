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
    momentBeamRotationNRFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient rotation equation (stressedFoam).

SourceFiles
    momentBeamRotationNRFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef momentBeamRotationNRFvPatchVectorField_H
#define momentBeamRotationNRFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class momentBeamRotationNRFvPatch Declaration
\*---------------------------------------------------------------------------*/

class momentBeamRotationNRFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Specified moment
        vectorField moment_;

        //- Traction time series
        interpolationTable<vector> momentSeries_;

public:

    //- Runtime type information
    TypeName("momentBeamRotationNR");


    // Constructors

        //- Construct from patch and internal field
        momentBeamRotationNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        momentBeamRotationNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  momentBeamRotationNRFvPatchVectorField onto a new patch
        momentBeamRotationNRFvPatchVectorField
        (
            const momentBeamRotationNRFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        momentBeamRotationNRFvPatchVectorField
        (
            const momentBeamRotationNRFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new momentBeamRotationNRFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        momentBeamRotationNRFvPatchVectorField
        (
            const momentBeamRotationNRFvPatchVectorField&,
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
                new momentBeamRotationNRFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return current moment vector
            const vectorField& moment() const
            {
                return moment_;
            }

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

            //const Pstream::commsTypes commsType = Pstream::blocking
            //const Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
