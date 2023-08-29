/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::flowRateOutletPressureFvPatchScalarField

Description
    Foam::flowRateOutletPressureFvPatchScalarField

SourceFiles
    flowRateOutletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateOutletPressureFvPatchScalarFields_H
#define flowRateOutletPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "Switch.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class flowRateOutletPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class flowRateOutletPressureFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to
        //  normalise the mass flux if neccessary
        word rhoName_;

        //- Flow rate fraction of the total inlet flow rate
        scalar flowRateFraction_;

        //- Flow rate fraction is not specified and
        //  will be calculated based on patch area
        bool calculatedFlowRateFraction_;
    
        //- Flow rate fractions is checked?
        bool checkedFlowRateFractions_;

        //- Inlet flow rate patch index
        labelList inletPatchIndices_;

        //- Current phi correction
        scalarField phiCorr_;

public:

    //- Runtime type information
    TypeName("flowRateOutletPressure");


    // Constructors

        //- Construct from patch and internal field
        flowRateOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateOutletPressureFvPatchScalarField onto a new patch
        flowRateOutletPressureFvPatchScalarField
        (
            const flowRateOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateOutletPressureFvPatchScalarField
        (
            const flowRateOutletPressureFvPatchScalarField&
#ifdef OPENFOAMFOUNDATION
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new flowRateOutletPressureFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        flowRateOutletPressureFvPatchScalarField
        (
            const flowRateOutletPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new flowRateOutletPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return constant reference to flow rate fraction
        const scalar& flowRateFraction() const
        {
            return flowRateFraction_;
        }
    
        //- Return constant reference to flow rate fraction
        bool calculatedFlowRateFraction() const
        {
            return calculatedFlowRateFraction_;
        }
    
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
    
        //- Manipulate matrix
        virtual void manipulateMatrix(fvMatrix<scalar>& matrix);
    
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
