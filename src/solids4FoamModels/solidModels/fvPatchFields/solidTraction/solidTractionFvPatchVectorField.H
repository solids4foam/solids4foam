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
    solidTractionFvPatchVectorField

Description
    Fixed traction boundary condition for solid solvers.

Author
    Philip Cardiff, UCD.
    Zeljko Tukovic, FSB Zagreb.

SourceFiles
    solidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionFvPatchVectorField_H
#define solidTractionFvPatchVectorField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        //- Enable non-orthogonal corrections
        const Switch nonOrthogonalCorrections_;

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Traction time series
        interpolationTable<vector> tractionSeries_;

        //- Pressure time series
        interpolationTable<scalar> pressureSeries_;

        //- Traction field given as a vol field
        autoPtr<volVectorField> tractionFieldPtr_;

        //- Pressure field given as a vol field
        autoPtr<volScalarField> pressureFieldPtr_;

        //- Second order accuracy correction
        const Switch secondOrder_;

        //- Switch so that the effective traction is enforced rather than the
        //  total traction. This only applies to porous flow based solids
        const Switch setEffectiveTraction_;

        //- Optional: under-relaxation factor
        const scalar relaxFac_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("solidTraction");


    // Constructors

        //- Construct from patch and internal field
        solidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidTractionFvPatchVectorField onto a new patch
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&
#ifdef OPENFOAMFOUNDATION
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidTractionFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidTractionFvPatchVectorField
        (
            const solidTractionFvPatchVectorField&,
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
                new solidTractionFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual scalarField& pressure()
            {
                return pressure_;
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

        //- Evaluate the patch field
        virtual void evaluate
        (
#ifdef OPENFOAMESIORFOUNDATION
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
