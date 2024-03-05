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

\*---------------------------------------------------------------------------*/

#include "flowRateOutletPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

#include "flowRateInletVelocityFvPatchVectorField.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateOutletPressureFvPatchScalarField::
flowRateOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    flowRateFraction_(1.0),
    calculatedFlowRateFraction_(false),
    checkedFlowRateFractions_(false),
    inletPatchIndices_(),
    phiCorr_(p.size(), 0)
{}


Foam::flowRateOutletPressureFvPatchScalarField::
flowRateOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    flowRateFraction_(dict.lookupOrDefault<scalar>("flowRateFraction", -1)),
    calculatedFlowRateFraction_(false),
    checkedFlowRateFractions_(false),
    inletPatchIndices_(),
    phiCorr_(p.size(), 0)
{
    fvPatchField<scalar>::operator=(patchInternalField());

    if (flowRateFraction_ < 0)
    {
        calculatedFlowRateFraction_ = true;
    }
    else
    {
        checkedFlowRateFractions_ = true;
    }
}


Foam::flowRateOutletPressureFvPatchScalarField::
flowRateOutletPressureFvPatchScalarField
(
    const flowRateOutletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    flowRateFraction_(ptf.flowRateFraction_),
    calculatedFlowRateFraction_(ptf.calculatedFlowRateFraction_),
    checkedFlowRateFractions_(false),
    inletPatchIndices_(ptf.inletPatchIndices_),
    phiCorr_(ptf.phiCorr_)
{}


#ifndef OPENFOAM_ORG
Foam::flowRateOutletPressureFvPatchScalarField::
flowRateOutletPressureFvPatchScalarField
(
    const flowRateOutletPressureFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    flowRateFraction_(ptf.flowRateFraction_),
    calculatedFlowRateFraction_(ptf.calculatedFlowRateFraction_),
    checkedFlowRateFractions_(ptf.checkedFlowRateFractions_),
    inletPatchIndices_(ptf.inletPatchIndices_),
    phiCorr_(ptf.phiCorr_)
{}
#endif


Foam::flowRateOutletPressureFvPatchScalarField::
flowRateOutletPressureFvPatchScalarField
(
    const flowRateOutletPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    flowRateFraction_(ptf.flowRateFraction_),
    calculatedFlowRateFraction_(ptf.calculatedFlowRateFraction_),
    checkedFlowRateFractions_(ptf.checkedFlowRateFractions_),
    inletPatchIndices_(ptf.inletPatchIndices_),
    phiCorr_(ptf.phiCorr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowRateOutletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!checkedFlowRateFractions_)
    {
        // Calculate flow-rate fractions if not specified

        const volScalarField& p =
            db().lookupObject<volScalarField>
            (
#ifdef OPENFOAM_NOT_EXTEND
                this->internalField().name()
#else
                this->dimensionedInternalField().name()
#endif
            );

        scalar totOutletArea = 0;
        scalar specFlowRateFraction = 0;

        forAll(p.boundaryField(), patchI)
        {
            if
            (
                isA<flowRateOutletPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
            )
            {
                const flowRateOutletPressureFvPatchScalarField& frop =
                    refCast<const flowRateOutletPressureFvPatchScalarField>
                    (
                        p.boundaryField()[patchI]
                    );

                if (frop.calculatedFlowRateFraction())
                {
                    totOutletArea +=
                        gSum(p.mesh().boundary()[patchI].magSf());
                }
                else
                {
                    specFlowRateFraction += frop.flowRateFraction();
                }
            }
        }

        if (totOutletArea < SMALL)
        {
            FatalErrorIn
            (
                "flowRateOutletPressureFvPatchScalarField::updateCoeffs()"
            ) << "Zero total outlet area.\n" << exit(FatalError);
        }

        if (flowRateFraction_ < 0)
        {
            scalar specArea =
                totOutletArea*specFlowRateFraction/
                (1.0-specFlowRateFraction);

            flowRateFraction_ =
                gSum(this->patch().magSf())/
                (totOutletArea + specArea);
        }

        checkedFlowRateFractions_ = true;
    }

    const volVectorField& U =
        db().lookupObject<volVectorField>(UName_);

    // Find inlet patch indices
    if (inletPatchIndices_.size() == 0)
    {
        DynamicList<label> ips;

        forAll(U.boundaryField(), patchI)
        {
            if
            (
                isA<flowRateInletVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
            {
                ips.append(patchI);
                // inletPatchIndex_ = patchI;
                // break;
            }
        }

        inletPatchIndices_ = ips;

        if (inletPatchIndices_.size() == 0)
        {
            FatalErrorIn
            (
                "flowRateOutletPressureFvPatchScalarField"
                "::updateCoeffs()"
            )   << "Can not find flowRateInletVelocityFvPatchVectorField"
                << exit(FatalError);
        }
    }

    // Scale outlet flow rates
    surfaceScalarField& phi =
        const_cast<surfaceScalarField&>
        (
            db().lookupObject<surfaceScalarField>(phiName_)
        );

    scalar inletFlowRate = 0;
    forAll(inletPatchIndices_, pI)
    {
        inletFlowRate +=
            gSum(phi.boundaryField()[inletPatchIndices_[pI]]);
    }

    scalar reqOutletFlowRate = -inletFlowRate*flowRateFraction_;

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& phip =
        phi.boundaryFieldRef()[this->patch().index()];
#else
    scalarField& phip =
        phi.boundaryField()[this->patch().index()];
#endif

    label size = this->patch().size();
    reduce(size, sumOp<label>());

    scalar adjustableFlowRateOut = 0;
    scalar adjustableFlowRateIn = 0;

    forAll(phip, faceI)
    {
        if (phip[faceI] > 0)
        {
            adjustableFlowRateOut += phip[faceI];
        }
        else
        {
            adjustableFlowRateIn -= phip[faceI];
        }
    }

    reduce(adjustableFlowRateOut, sumOp<scalar>());
    reduce(adjustableFlowRateIn, sumOp<scalar>());

    scalar netFlowRate =
        adjustableFlowRateOut
      - adjustableFlowRateIn
      - reqOutletFlowRate;

    scalar totAdjustableFlowRate =
        adjustableFlowRateOut + adjustableFlowRateIn;

    phiCorr_ = 0;

    if (totAdjustableFlowRate < SMALL)
    {
        if (mag(netFlowRate) > SMALL)
        {
            // Uniform flow-rate
            phiCorr_ = netFlowRate/size;
        }
    }
    else
    {
        forAll(phiCorr_, faceI)
        {
            phiCorr_[faceI] =
               -netFlowRate*mag(phip[faceI])/
                totAdjustableFlowRate;
        }
    }

    // Correct phi
    phip += phiCorr_;

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        FatalErrorIn
        (
            "flowRateOutletPressureFvPatchScalarField"
            "::updateCoeffs()"
        )   << this->type() << " not implemented for compressible flow"
            << exit(FatalError);
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


void Foam::flowRateOutletPressureFvPatchScalarField::
manipulateMatrix(fvMatrix<scalar>& matrix)
{
    const labelList& faceCells = this->patch().faceCells();

    scalarField& source = matrix.source();

    forAll(faceCells, faceI)
    {
        source[faceCells[faceI]] += phiCorr_[faceI];
    }

    fvPatchScalarField::manipulateMatrix(matrix);
}


void Foam::flowRateOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

#ifdef OPENFOAM_COM
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
#else
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
#endif

    os.writeKeyword("flowRateFraction")
        << flowRateFraction_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        flowRateOutletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
