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
    forceBeamDisplacementNRFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient rotation equation (stressedFoam).

SourceFiles
    forceBeamDisplacementNRFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef forceBeamDisplacementNRFvPatchVectorField_H
#define forceBeamDisplacementNRFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class forceBeamDisplacementNRFvPatch Declaration
\*---------------------------------------------------------------------------*/

class forceBeamDisplacementNRFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Specified force
        vectorField force_;
        vectorField curForce_;

        //- Traction time series
        interpolationTable<vector> forceSeries_;

        //- Transition function time series
        interpolationTable<scalar> transitionSeries_;

public:

    //- Runtime type information
    TypeName("forceBeamDisplacementNR");


    // Constructors

        //- Construct from patch and internal field
        forceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        forceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  forceBeamDisplacementNRFvPatchVectorField onto a new patch
        forceBeamDisplacementNRFvPatchVectorField
        (
            const forceBeamDisplacementNRFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        forceBeamDisplacementNRFvPatchVectorField
        (
            const forceBeamDisplacementNRFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new forceBeamDisplacementNRFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        forceBeamDisplacementNRFvPatchVectorField
        (
            const forceBeamDisplacementNRFvPatchVectorField&,
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
                new forceBeamDisplacementNRFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

           // Return current force vector
           const vectorField& force() const
           {
               return curForce_;
           }

           // Return current force vector
           vectorField& force()
           {
               return curForce_;
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
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
