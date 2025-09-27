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

#include "fixedGradientCorrectedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "patchCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedGradientCorrectedFvPatchScalarField::
fixedGradientCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    nonOrthogonalCorrections_(true),
    secondOrder_(false),
    extrapolateValue_(false)
{
    fvPatchScalarField::operator=(patchInternalField());
    gradient() = 0.0;
}


fixedGradientCorrectedFvPatchScalarField::
fixedGradientCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    ),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    extrapolateValue_(dict.lookupOrDefault<Switch>("extrapolateValue", false))
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        gradient() = 0.0;
    }

    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }

    if (secondOrder_)
    {
        Info<< "    second order correction" << endl;
    }

    if (extrapolateValue_)
    {
        Info<< "    extrapolating patch values" << endl;
    }
}


fixedGradientCorrectedFvPatchScalarField::
fixedGradientCorrectedFvPatchScalarField
(
    const fixedGradientCorrectedFvPatchScalarField& pvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    extrapolateValue_(pvf.extrapolateValue_)
{}

#ifndef OPENFOAM_ORG
fixedGradientCorrectedFvPatchScalarField::
fixedGradientCorrectedFvPatchScalarField
(
    const fixedGradientCorrectedFvPatchScalarField& pvf
)
:
    fixedGradientFvPatchScalarField(pvf),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    extrapolateValue_(pvf.extrapolateValue_)
{}
#endif

fixedGradientCorrectedFvPatchScalarField::
fixedGradientCorrectedFvPatchScalarField
(
    const fixedGradientCorrectedFvPatchScalarField& pvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    extrapolateValue_(pvf.extrapolateValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedGradientCorrectedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedGradientCorrectedFvPatchScalarField::rmap
(
    const fvPatchScalarField& pvf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(pvf, addr);

    // const fixedGradientCorrectedFvPatchScalarField& rpvf =
    //     refCast<const fixedGradientCorrectedFvPatchScalarField>(pvf);
}


// Update the coefficients associated with the patch field
void fixedGradientCorrectedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fixedGradientCorrectedFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if
    (
        nonOrthogonalCorrections_
     && db().foundObject<volVectorField>("grad(" + internalField().name() + ")")
    )
    {
        // Lookup the gradient field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>
            (
            #ifdef OPENFOAM_NOT_EXTEND
                "grad(" + internalField().name() + ")"
            #else
                "grad(" + dimensionedInternalField().name() + ")"
            #endif
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        if (secondOrder_)
        {
            // Face unit normals
            const vectorField n(patch().nf());

            // Correction to internal field cells
            const scalarField dUP(k & gradField.patchInternalField());

            // Normal gradient at internal field cells
            const scalarField nGradUP(n & gradField.patchInternalField());

            Field<scalar>::operator=
            (
                patchInternalField()
              + dUP
              + 0.5*(gradient() + nGradUP)/patch().deltaCoeffs()
            );
        }
        else if (extrapolateValue_)
        {
            // Face unit normals
            const vectorField n(patch().nf());

            Field<scalar>::operator=
            (
                patchInternalField()
              + (k & gradField.patchInternalField())
              + (n & gradField.patchInternalField())/patch().deltaCoeffs()
            );
        }
        else
        {
            Field<scalar>::operator=
            (
                patchInternalField()
              + (k & gradField.patchInternalField())
              + gradient()/patch().deltaCoeffs()
            );
        }
    }
    else
    {
        WarningInFunction
            << "Non-orthogonal corrections not applied on " << patch().name()
            << endl;

        Field<scalar>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<scalar>::evaluate();
}


void fixedGradientCorrectedFvPatchScalarField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //fixedGradientFvPatchScalarField::write(os);
    fvPatchScalarField::write(os);

    os.writeKeyword("nonOrthogonalCorrections")
        << nonOrthogonalCorrections_ << token::END_STATEMENT << nl;
    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
    os.writeKeyword("extrapolateValue")
        << extrapolateValue_ << token::END_STATEMENT << nl;

#ifdef OPENFOAM_ORG
    writeEntry(os, "value", *this);
    writeEntry(os, "gradient", gradient());
#else
    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedGradientCorrectedFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
