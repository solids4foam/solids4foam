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
    myFixedRotationFvPatchVectorField

Description
    Specifying a beam rotation value at the patch

SourceFiles
    fixedBeamRotationFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedBeamRotationFvPatchVectorField_H
#define fixedBeamRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fixedBeamRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedBeamRotationFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Constant total rotation field
        vectorField totalTheta_;

        //- Theta time series
        interpolationTable<vector> thetaSeries_;

        //- Current time index
        label timeIndex_;

public:

    //- Runtime type information
    TypeName("fixedBeamRotation");


    // Constructors

        //- Construct from patch and internal field
        fixedBeamRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedBeamRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedBeamRotationFvPatchVectorField onto a new patch
        fixedBeamRotationFvPatchVectorField
        (
            const fixedBeamRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedBeamRotationFvPatchVectorField
        (
            const fixedBeamRotationFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedBeamRotationFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedBeamRotationFvPatchVectorField
        (
            const fixedBeamRotationFvPatchVectorField&,
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
                new fixedBeamRotationFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        //- Update the coefficients associated with the patch field
            virtual vectorField& totalTheta()
            {
                return totalTheta_;
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

        // //- Evaluate the patch field
        // virtual void evaluate
        // (
        //     const Pstream::commsTypes commsType = Pstream::blocking
        // );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
