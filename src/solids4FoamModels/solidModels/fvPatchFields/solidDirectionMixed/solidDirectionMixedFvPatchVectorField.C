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

#include "solidDirectionMixedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "patchCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    nonOrthogonalCorrections_(true),
    secondOrder_(false),
    limitCoeff_(1.0)
{
    // solidModel may not exist here, so we cannot look it up
    // Lookup the solidModel object
    // const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // if (solMod.solidModelDict().found("snGradLimitCoeff"))
    // {
    //     limitCoeff_ =
    //         readScalar
    //         (
    //             solMod.solidModelDict().lookup("snGradLimitCoeff")
    //         );

    //     Info<< "snGradLimitCoeff: " << limitCoeff_ << endl;
    // }
}


solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const solidDirectionMixedFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    limitCoeff_(pvf.limitCoeff_)
{}


solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF, dict),
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    ),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0))
{
    Info<< "Creating " << type() << " boundary condition" << endl;
    directionMixedFvPatchVectorField::evaluate();

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    if (solMod.solidModelDict().found("snGradLimitCoeff"))
    {
        limitCoeff_ =
            readScalar
            (
                solMod.solidModelDict().lookup("snGradLimitCoeff")
            );

        Info<< "snGradLimitCoeff: " << limitCoeff_ << endl;
    }

    Info<< "    Limiter coefficient: " << limitCoeff_ << nl
        << "    Second order correction: " << secondOrder_ << endl;
}

solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const solidDirectionMixedFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    limitCoeff_(pvf.limitCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidDirectionMixedFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidDirectionMixedFvPatchVectorField::rmap
(
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(pvf, addr);
}


void solidDirectionMixedFvPatchVectorField::updateCoeffs()
{
    directionMixedFvPatchVectorField::updateCoeffs();
}

void solidDirectionMixedFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (nonOrthogonalCorrections_)
    {
        // Lookup the gradient field
        const fvPatchField<tensor>& gradD =
            patch().lookupPatchField<volTensorField, tensor>
            (
            #ifdef OPENFOAM_NOT_EXTEND
                "grad(" + internalField().name() + ")"
            #else
                "grad(" + dimensionedInternalField().name() + ")"
            #endif
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        // Calc limited snGrad correction

        vectorField snGradCorrection
        (
          - (k & gradD.patchInternalField())*patch().deltaCoeffs()
        );

        if (limitCoeff_ < (1.0 - SMALL))
        {
            const vectorField uncorrectedSnGrad
            (
                (
                   *this - patchInternalField()
                )*patch().deltaCoeffs()
            );

            const scalarField limiter
            (
                min
                (
                    limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
                   /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
                    1.0
                )
            );

            snGradCorrection *= limiter;
        }

        const Field<vector> normalValue(transform(valueFraction(), refValue()));

        Field<vector> gradValue(patch().size());

        if (secondOrder_)
        {
            // Unit normals
            const vectorField n(patch().nf());

            // Normal gradient at internal cell
            const vectorField nGradDP(n & gradD.patchInternalField());

            gradValue =
                patchInternalField()
              + (k & gradD.patchInternalField())
              + 0.5*(nGradDP + refGrad())/patch().deltaCoeffs();
        }
        else
        {
            gradValue =
                patchInternalField()
              - snGradCorrection/patch().deltaCoeffs()
              + refGrad()/patch().deltaCoeffs();
        }

        const Field<vector> transformGradValue
        (
            transform(I - valueFraction(), gradValue)
        );

        Field<vector>::operator=(normalValue + transformGradValue);
    }
    else
    {
        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        const Field<vector> normalValue(transform(valueFraction(), refValue()));

        const Field<vector> gradValue
        (
            patchInternalField() + refGrad()/patch().deltaCoeffs()
        );

        const Field<vector> transformGradValue
        (
            transform(I - valueFraction(), gradValue)
        );

        Field<vector>::operator=(normalValue + transformGradValue);
    }

    fvPatchField<vector>::evaluate();
}


Foam::tmp<Foam::Field<vector> >
solidDirectionMixedFvPatchVectorField::snGrad() const
{
    if (nonOrthogonalCorrections_)
    {
        // Lookup the gradient field
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
            #ifdef OPENFOAM_NOT_EXTEND
                "grad(" + internalField().name() + ")"
            #else
                "grad(" + dimensionedInternalField().name() + ")"
            #endif
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        // Patch internal field
        const Field<vector> pif(this->patchInternalField());

        // Patch normal value
        const Field<vector> normalValue
        (
            transform(this->valueFraction(), this->refValue())
        );

        vectorField snGradCorrection
        (
          - (k & gradField.patchInternalField())
           *patch().deltaCoeffs()
        );

        if (limitCoeff_ < (1.0 - SMALL))
        {
            const vectorField uncorrectedSnGrad
            (
                (
                   *this
                  - patchInternalField()
                )*patch().deltaCoeffs()
            );

            const scalarField limiter
            (
                (
                    min
                    (
                        limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
                       /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
                        1.0
                    )
                )
            );

            snGradCorrection *= limiter;
        }

        Field<vector> gradValue(patch().size(), vector::zero);

        if (secondOrder_)
        {
            // Unit normals
            const vectorField n(patch().nf());

            // Normal gradient at internal cell
            const vectorField nGradDP(n & gradField.patchInternalField());

            gradValue =
                patchInternalField()
              + (k & gradField.patchInternalField())
              + 0.5*(nGradDP + refGrad())/patch().deltaCoeffs();
        }
        else
        {
            gradValue =
                pif
              - snGradCorrection/patch().deltaCoeffs()
              + refGrad()/patch().deltaCoeffs();
        }

        const Field<vector> transformGradValue
        (
            transform(I - this->valueFraction(), gradValue)
        );

        if (secondOrder_)
        {
            // Unit normals
            const vectorField n(patch().nf());

            // Normal gradient at internal cell
            const vectorField nGradDP(n & gradField.patchInternalField());

            return
                2.0
               *(
                   normalValue + transformGradValue
                 - (pif + (k & gradField.patchInternalField()))
               )*patch().deltaCoeffs()
             - nGradDP;
        }
        else
        {
            return
            (
                normalValue + transformGradValue
              - (
                    pif - snGradCorrection/this->patch().deltaCoeffs()
                )
            )*patch().deltaCoeffs();
        }
    }

    // else, no non-orthogonal correction

    // Patch internal field
    const Field<vector> pif(this->patchInternalField());

    // Patch normal value
    const Field<vector> normalValue
    (
        transform(this->valueFraction(), this->refValue())
    );

    const Field<vector> gradValue
    (
        pif + refGrad()/patch().deltaCoeffs()
    );

    const Field<vector> transformGradValue
    (
        transform(I - this->valueFraction(), gradValue)
    );

    return
    (
        normalValue + transformGradValue - pif
    )*patch().deltaCoeffs();
}

// Write
void solidDirectionMixedFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);

    os.writeKeyword("nonOrthogonalCorrections")
        << nonOrthogonalCorrections_ << token::END_STATEMENT << nl;
    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
    os.writeKeyword("limitCoeff")
        << limitCoeff_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidDirectionMixedFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
