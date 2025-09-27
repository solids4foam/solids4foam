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

#include "fixedGradientCorrectedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "patchCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedGradientCorrectedFvPatchVectorField::
fixedGradientCorrectedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    nonOrthogonalCorrections_(true),
    secondOrder_(false),
    extrapolateValue_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


fixedGradientCorrectedFvPatchVectorField::
fixedGradientCorrectedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
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
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
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


fixedGradientCorrectedFvPatchVectorField::
fixedGradientCorrectedFvPatchVectorField
(
    const fixedGradientCorrectedFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    extrapolateValue_(pvf.extrapolateValue_)
{}

#ifndef OPENFOAM_ORG
fixedGradientCorrectedFvPatchVectorField::
fixedGradientCorrectedFvPatchVectorField
(
    const fixedGradientCorrectedFvPatchVectorField& pvf
)
:
    fixedGradientFvPatchVectorField(pvf),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    extrapolateValue_(pvf.extrapolateValue_)
{}
#endif

fixedGradientCorrectedFvPatchVectorField::
fixedGradientCorrectedFvPatchVectorField
(
    const fixedGradientCorrectedFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    secondOrder_(pvf.secondOrder_),
    extrapolateValue_(pvf.extrapolateValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedGradientCorrectedFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedGradientCorrectedFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(pvf, addr);

    // const fixedGradientCorrectedFvPatchVectorField& rpvf =
    //     refCast<const fixedGradientCorrectedFvPatchVectorField>(pvf);
}


// Update the coefficients associated with the patch field
void fixedGradientCorrectedFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void fixedGradientCorrectedFvPatchVectorField::evaluate
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
     && db().foundObject<volTensorField>("grad(" + internalField().name() + ")")
    )
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

        if (secondOrder_)
        {
            // Face unit normals
            const vectorField n(patch().nf());

            // Correction to internal field cells
            const vectorField dUP(k & gradField.patchInternalField());

            // Normal gradient at internal field cells
            const vectorField nGradUP(n & gradField.patchInternalField());

            Field<vector>::operator=
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

            Field<vector>::operator=
            (
                patchInternalField()
              + (k & gradField.patchInternalField())
              + (n & gradField.patchInternalField())/patch().deltaCoeffs()
            );
        }
        else
        {
            Field<vector>::operator=
            (
                patchInternalField()
              + (k & gradField.patchInternalField())
              + gradient()/patch().deltaCoeffs()
            );
        }
    }
    else
    {
        // No non-orthogonal correction
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


void fixedGradientCorrectedFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //fixedGradientFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

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

makePatchTypeField(fvPatchVectorField, fixedGradientCorrectedFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
