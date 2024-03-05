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

#include "fixedTemperatureGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedTemperatureGradientFvPatchScalarField::
fixedTemperatureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    fieldName_("undefined")
{}


fixedTemperatureGradientFvPatchScalarField::
fixedTemperatureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    fieldName_(dimensionedInternalField().name())
{}


fixedTemperatureGradientFvPatchScalarField::
fixedTemperatureGradientFvPatchScalarField
(
    const fixedTemperatureGradientFvPatchScalarField& stpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(stpvf, p, iF, mapper),
    fieldName_(stpvf.fieldName_)
{}

#ifndef OPENFOAM_ORG
fixedTemperatureGradientFvPatchScalarField::
fixedTemperatureGradientFvPatchScalarField
(
    const fixedTemperatureGradientFvPatchScalarField& stpvf
)
:
    fixedGradientFvPatchScalarField(stpvf),
    fieldName_(stpvf.fieldName_)
{}
#endif

fixedTemperatureGradientFvPatchScalarField::
fixedTemperatureGradientFvPatchScalarField
(
    const fixedTemperatureGradientFvPatchScalarField& stpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(stpvf, iF),
    fieldName_(stpvf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedTemperatureGradientFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedTemperatureGradientFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
// void fixedTemperatureGradientFvPatchScalarField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }

//     fixedGradientFvPatchScalarField::updateCoeffs();
// }


void fixedTemperatureGradientFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup grad from solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + fieldName_ + ")"
        );

    vectorField n(patch().nf());
    vectorField delta(patch().delta());

    // Correction vectors
    vectorField k = delta - n*(n&delta);

    Field<scalar>::operator=
    (
        this->patchInternalField()
      + (k&gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<scalar>::evaluate();
}

// Write
void fixedTemperatureGradientFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedTemperatureGradientFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
