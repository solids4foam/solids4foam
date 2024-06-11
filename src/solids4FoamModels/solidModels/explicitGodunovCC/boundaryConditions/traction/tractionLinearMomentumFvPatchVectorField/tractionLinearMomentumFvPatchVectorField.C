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

#include "tractionLinearMomentumFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    loadingType_("none"),
    t_P_(vector::zero),
    p_P_(0.0),
    pressureSeries_()
{}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    loadingType_(dict.lookupOrDefault<word>("loadingType", "none")),
    t_P_(dict.lookupOrDefault("traction", vector::zero)),
    p_P_(dict.lookupOrDefault("pressure", 0.0)),
    pressureSeries_(dict.subDict("pressureSeries"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (loadingType_ == "traction" && !dict.found("traction") == 1)
    {
        FatalErrorIn("tractionLinearMomentumFvPatchVectorField")
            << "Keyword 'traction' is undefined in dictionary '"
            << iF.name() << ".boundaryField." << p.name()
            << "' for patch type '" << this->type() << "'." << nl
            << exit(FatalError);
    }

    if (loadingType_ == "pressure" && !dict.found("pressure") == 1)
    {
        FatalErrorIn ("tractionLinearMomentumFvPatchVectorField")
            << "Keyword 'pressure' is undefined in dictionary '"
            << iF.name() << ".boundaryField." << p.name()
            << "' for patch type '" << this->type() << "'." << nl
            << exit(FatalError);
    }

    updateCoeffs();
}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const tractionLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    loadingType_(ptf.loadingType_),
    t_P_(ptf.t_P_),
    p_P_(ptf.p_P_),
    pressureSeries_(ptf.pressureSeries_)
{}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const tractionLinearMomentumFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    loadingType_(rifvpvf.loadingType_),
    t_P_(rifvpvf.t_P_),
    p_P_(rifvpvf.p_P_),
    pressureSeries_(rifvpvf.pressureSeries_)
{}


tractionLinearMomentumFvPatchVectorField::
tractionLinearMomentumFvPatchVectorField
(
    const tractionLinearMomentumFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    loadingType_(rifvpvf.loadingType_),
    t_P_(rifvpvf.t_P_),
    p_P_(rifvpvf.p_P_),
    pressureSeries_(rifvpvf.pressureSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionLinearMomentumFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void tractionLinearMomentumFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    const fvsPatchField<vector>& t_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("t_M");

    const fvsPatchField<tensor>& S_t_ =
        patch().lookupPatchField<surfaceTensorField, tensor>("S_t");

    fvsPatchField<vector> lm_C(lm_M_);

    if (loadingType_ == "none")
    {
        lm_C = lm_M_ + (S_t_ & (-t_M_));
    }

    else if (loadingType_ == "traction")
    {
        lm_C = lm_M_ + (S_t_ & ((t_P_) - t_M_));
    }

    else if (loadingType_ == "pressure")
    {
        p_P_ = pressureSeries_(db().time().timeOutputValue());

        const fvsPatchField<vector>& n_ =
            patch().lookupPatchField<surfaceVectorField, vector>("n");

        lm_C = lm_M_ + (S_t_ & ((-p_P_*n_) - t_M_));
    }

    else
    {
        FatalErrorIn
        (
            "tractionLinearMomentumFvPatchVectorField::updateCoeffs()"
        )   << "Unknown traction type '" << loadingType_ << "'" << nl
            << "Valid types are: 'none', traction' and 'pressure'" << nl
            << exit(FatalError);
    }

    this->operator==(lm_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void tractionLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("loadingType") << loadingType_ << token::END_STATEMENT
        << nl;
    os.writeKeyword("traction") << t_P_ << token::END_STATEMENT << nl;
    os.writeKeyword("pressure") << p_P_ << token::END_STATEMENT << nl;
    // writeEntry(os, "value", *this);
    writeEntry("value", os);

    os.writeKeyword("pressureSeries") << nl;
    os << token::BEGIN_BLOCK << nl;
    pressureSeries_.write(os);
    os << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionLinearMomentumFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
