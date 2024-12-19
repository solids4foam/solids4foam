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
    explicitSolidTractionLinearMomentumFvPatchVectorField

Description
    Fixed traction boundary condition for explicitGodunovCC solid solver
    for linear momentum field.

Original Authors
    Philip Cardiff, UCD.
    Zeljko Tukovic, FSB Zagreb.

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

SourceFiles
    explicitSolidTractionLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef explicitSolidTractionLinearMomentumFvPatchVectorField_H
#define explicitSolidTractionLinearMomentumFvPatchVectorField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class explicitSolidTractionLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class explicitSolidTractionLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
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
    TypeName("explicitSolidTractionLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        explicitSolidTractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        explicitSolidTractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  explicitSolidTractionLinearMomentumFvPatchVectorField onto a new patch
        explicitSolidTractionLinearMomentumFvPatchVectorField
        (
            const explicitSolidTractionLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        explicitSolidTractionLinearMomentumFvPatchVectorField
        (
            const explicitSolidTractionLinearMomentumFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new explicitSolidTractionLinearMomentumFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        explicitSolidTractionLinearMomentumFvPatchVectorField
        (
            const explicitSolidTractionLinearMomentumFvPatchVectorField&,
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
                new explicitSolidTractionLinearMomentumFvPatchVectorField(*this, iF)
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
#ifdef OPENFOAM_NOT_EXTEND
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
