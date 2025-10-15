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

#include "movingTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
	lm_P_(vector::zero)
{}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
	lm_P_(dict.lookupOrDefault<vector>
    (
        "linearMomentum",
        vector::zero)
    )
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
	updateCoeffs();
}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const movingTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
	lm_P_(ptf.lm_P_)
{}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const movingTractionFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    lm_P_(rifvpvf.lm_P_)
{}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const movingTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    lm_P_(rifvpvf.lm_P_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void movingTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    const fvsPatchField<vector>& t_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("t_M");

    const fvsPatchField<tensor>& S_lm_ =
        patch().lookupPatchField<surfaceTensorField, tensor>("S_lm");

    fvsPatchField<vector> t_C(lm_M_);
    t_C = t_M_ + (S_lm_ & (lm_P_ - lm_M_));

    this->operator==(t_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void movingTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("lm_P_") << lm_P_ << token::END_STATEMENT << nl;
    // writeEntry(os, "value", *this);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
