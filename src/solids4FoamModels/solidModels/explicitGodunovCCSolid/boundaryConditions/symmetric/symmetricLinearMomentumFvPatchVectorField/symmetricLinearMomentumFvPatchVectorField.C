/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "symmetricLinearMomentumFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

symmetricLinearMomentumFvPatchVectorField::
symmetricLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    t_P_(vector::zero),
    tRamp_(VSMALL)
{}


symmetricLinearMomentumFvPatchVectorField::
symmetricLinearMomentumFvPatchVectorField
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


symmetricLinearMomentumFvPatchVectorField::
symmetricLinearMomentumFvPatchVectorField
(
    const symmetricLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    t_P_(ptf.t_P_),
    tRamp_(ptf.tRamp_)
{}


symmetricLinearMomentumFvPatchVectorField::
symmetricLinearMomentumFvPatchVectorField
(
    const symmetricLinearMomentumFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    t_P_(rifvpvf.t_P_),
    tRamp_(rifvpvf.tRamp_)
{}


symmetricLinearMomentumFvPatchVectorField::
symmetricLinearMomentumFvPatchVectorField
(
    const symmetricLinearMomentumFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    t_P_(rifvpvf.t_P_),
    tRamp_(rifvpvf.tRamp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void symmetricLinearMomentumFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void symmetricLinearMomentumFvPatchVectorField::updateCoeffs()
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

    const fvPatchField<scalar>& Us_ =
        patch().lookupPatchField<volScalarField, scalar>("Us");

    fvsPatchField<vector> lm_C(lm_M_);

    scalar ramp = this->db().time().value()/tRamp_;
    if (this->db().time().value() >= tRamp_)
    {
        ramp = 1.0;
    }

    lm_C = (tensor::I - n_*n_) & (lm_M_ + (((ramp*t_P_) - t_M_)/Us_));

    this->operator==(lm_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void symmetricLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("traction") << t_P_ << token::END_STATEMENT << nl;
    // writeEntry(os, "value", *this);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    symmetricLinearMomentumFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
