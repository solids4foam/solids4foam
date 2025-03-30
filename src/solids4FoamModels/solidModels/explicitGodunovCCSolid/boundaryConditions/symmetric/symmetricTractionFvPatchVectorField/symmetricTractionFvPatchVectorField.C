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

#include "symmetricTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    t_P_(vector::zero),
    tRamp_(VSMALL)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    t_P_(vector::zero),
    tRamp_(VSMALL)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    t_P_ = dict.lookupOrDefault<vector>("traction", vector::zero);
    tRamp_ = dict.lookupOrDefault("rampEndTime", VSMALL);

    updateCoeffs();
}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    t_P_(ptf.t_P_),
    tRamp_(ptf.tRamp_)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    t_P_(rifvpvf.t_P_),
    tRamp_(rifvpvf.tRamp_)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    t_P_(rifvpvf.t_P_),
    tRamp_(rifvpvf.tRamp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void symmetricTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void symmetricTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    const fvsPatchField<vector>& t_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("t_M");

    const fvsPatchField<vector>& n_ =
        patch().lookupPatchField<surfaceVectorField, vector>("n");

    const fvPatchField<scalar>& Up_ =
        patch().lookupPatchField<volScalarField, scalar>("Up");

    fvsPatchField<vector> t_C(lm_M_);

    scalar ramp = this->db().time().value()/tRamp_;
    if (this->db().time().value() >= tRamp_)
    {
        ramp = 1.0;
    }

    t_C = ((n_*n_) & (t_M_ - Up_*lm_M_)) + ((tensor::I - n_*n_) & (ramp*t_P_));

    this->operator==(t_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void symmetricTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("traction") << t_P_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
    //writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    symmetricTractionFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
