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
#include "CompactListList.H"
#include "Switch.H"
#include "interpolationTable.H"
#include "fieldWithHistory.H"

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

        //- Use undeformed areas flag
        //  Note that this flag is simply stored here and it is the
        //  responsibility of the solid model to lookup and use this
        const Switch useUndeformedArea_;

        //- Traction field and, optionally, old-time traction fields
        mutable fieldWithHistory traction_;

        //- Optional: traction-rate field and old-time traction rate fields
        mutable autoPtr<fieldWithHistory> tractionRatePtr_;

        //- Pressure
        mutable scalarField pressure_;

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

        //- Extrapolate values to the patch using the patch internal gradient
        //  field. Alternatively, if false, use the gradient at the patch
        Switch extrapolateValue_;

        //- Switch so that the effective traction is enforced rather than the
        //  total traction. This only applies to porous flow based solids
        const Switch setEffectiveTraction_;

        //- Optional: under-relaxation factor
        const scalar relaxFac_;

        //- Current time index
        mutable label curTimeIndex_;


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
#ifdef OPENFOAM_ORG
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

            virtual const fieldWithHistory& tractionHistory() const
            {
                return traction_;
            }

            virtual fieldWithHistory& tractionHistory()
            {
                return traction_;
            }

            virtual const fieldWithHistory& tractionRateHistory() const
            {
                if (!tractionRatePtr_.valid())
                {
                    tractionRatePtr_.reset
                    (
                        new fieldWithHistory(traction_.size())
                    );
                }

                return tractionRatePtr_();
            }

            virtual fieldWithHistory& tractionRateHistory()
            {
                if (!tractionRatePtr_.valid())
                {
                    tractionRatePtr_.reset
                    (
                        new fieldWithHistory(traction_.size())
                    );
                }

                return tractionRatePtr_();
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual scalarField& pressure()
            {
                return pressure_;
            }

            // Non-const access to the extrapolateValue flag
            Switch& extrapolateValue()
            {
                return extrapolateValue_;
            }

            //- Const access to the useUndeformedArea_ flag
            const Switch& useUndeformedArea() const
            {
                return useUndeformedArea_;
            }

#ifndef FOAMEXTEND
            //- Return face quadrature values
            //  Used for high order solvers, argument is position of points
            virtual autoPtr<CompactListList<vector>>
            evaluateQuadrature() const;
#endif

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
