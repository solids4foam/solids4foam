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
    tractionPressureDisplacementFvPatchVectorField

Description
    Traction boundary condition for the displacement-increment field in
    coupledPressureDisplacementSolid cases.

    The matching scalar pressure boundary condition is pressureTraction.

SourceFiles
    tractionPressureDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionPressureDisplacementFvPatchVectorField_H
#define tractionPressureDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "interpolationTable.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionPressureDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tractionPressureDisplacementFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private Data

        vectorField traction_;
        scalarField pressure_;

        //- Traction time series
        interpolationTable<vector> tractionSeries_;

        //- Pressure time series
        interpolationTable<scalar> pressureSeries_;

        //- Secondary gradient
        vectorField secondaryGradient_;

        //- Rhie-Chow displacement correction
        scalarField rhieChowDisplCorr_;

        //- Block coupled solution procedure
        Switch blockCoupled_;

        //- Interpolation points indices
        PtrList<DynamicList<label> > iPointsIndices_;

        //- Interpolation points
        PtrList<DynamicList<vector> > iPoints_;

        //- Interpolation matrices
        PtrList<scalarRectangularMatrix> iMatrices_;

        //- Far neighbour cell
        Field<label> PP_;

        //- Frozen traction
        bool frozenTraction_;

    // Private Functions

        //- Create interpolation points
        void makeInterpolationPoints();

        //- Patch extrapolated derivative
        tmp<vectorField> extrapolatedDerivative
        (
            const volTensorField& gradU,
            const vectorField& dir
        );

        //- Patch extrapolated derivative
        tmp<tensorField> extrapolatedGradient
        (
            const volTensorField& gradU
        );

        //- Patch extrapolated derivative
        tmp<vectorField> extrapolatedDisplacement();

        //- Deformation gradient corresponding to mean configuration
        tmp<tensorField> FM(const bool refConfig) const;
        tmp<tensorField> intFM(const bool refConfig) const;

        //- Displacement correction
        tmp<vectorField> dispCorr(const bool refConfig) const;

        //- Displacement correction
        tmp<vectorField> patchInternalSnGrad(const bool refConfig) const;

        //- Patch gradient
        tmp<vectorField> patchGradient() const;

        //- Patch gradient
        tmp<vectorField> patchGradient2() const;

        //- Patch value
        tmp<vectorField> patchValue() const;

        //- Patch value
        tmp<vectorField> valueBndCoeffs() const;

        //- Transformed normal derivative
        tmp<vectorField> transformedGradient() const;

public:

    //- Runtime type information
    TypeName("tractionPressureDisplacement");


    // Constructors

        //- Construct from patch and internal field
        tractionPressureDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionPressureDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  tractionPressureDisplacementFvPatchVectorField onto a new patch
        tractionPressureDisplacementFvPatchVectorField
        (
            const tractionPressureDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionPressureDisplacementFvPatchVectorField
        (
            const tractionPressureDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionPressureDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionPressureDisplacementFvPatchVectorField
        (
            const tractionPressureDisplacementFvPatchVectorField&,
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
                new tractionPressureDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            bool frozenTraction()
            {
                return frozenTraction_;
            }

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


            const scalarField& rhieChowDisplCorr() const
            {
                return rhieChowDisplCorr_;
            }

            scalarField& rhieChowDisplCorr()
            {
                return rhieChowDisplCorr_;
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
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector> > valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector> > valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
