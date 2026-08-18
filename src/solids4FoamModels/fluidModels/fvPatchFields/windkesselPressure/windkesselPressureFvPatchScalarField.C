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

#include "windkesselPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windkesselPressureFvPatchScalarField::
windkesselPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    po_(0),
    poo_(0),
    Qo_(0),
    Qoo_(0),
    R_(1),
    Rch_(1),
    C_(1),
    rho_(1),
    timeIndex_(-1)
{}


Foam::windkesselPressureFvPatchScalarField::
windkesselPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    po_(dict.lookupOrDefault<scalar>("po", 0)),
    poo_(dict.lookupOrDefault<scalar>("poo", 0)),
    Qo_(dict.lookupOrDefault<scalar>("Qo", 0)),
    Qoo_(dict.lookupOrDefault<scalar>("Qoo", 0)),
    R_(dict.lookupOrDefault<scalar>("R", 1)),
    Rch_(dict.lookupOrDefault<scalar>("Rch", 1)),
    C_(dict.lookupOrDefault<scalar>("C", 1)),
    rho_(dict.lookupOrDefault<scalar>("rho", 1)),
    timeIndex_(-1)
{
    Info<< "Properties of Windkessel model, R = " << R_
        << ", Rch = " << Rch_ << ", C = " << C_
        << ", rho = " << rho_ << endl;
}

Foam::windkesselPressureFvPatchScalarField::
windkesselPressureFvPatchScalarField
(
    const windkesselPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    po_(ptf.po_),
    poo_(ptf.poo_),
    Qo_(ptf.Qo_),
    Qoo_(ptf.Qoo_),
    R_(ptf.R_),
    Rch_(ptf.Rch_),
    C_(ptf.C_),
    rho_(ptf.rho_),
    timeIndex_(ptf.timeIndex_)
{}


Foam::windkesselPressureFvPatchScalarField::
windkesselPressureFvPatchScalarField
(
    const windkesselPressureFvPatchScalarField& wkpsf
)
:
    fixedValueFvPatchScalarField(wkpsf),
    po_(wkpsf.po_),
    poo_(wkpsf.poo_),
    Qo_(wkpsf.Qo_),
    Qoo_(wkpsf.Qoo_),
    R_(wkpsf.R_),
    Rch_(wkpsf.Rch_),
    C_(wkpsf.C_),
    rho_(wkpsf.rho_),
    timeIndex_(wkpsf.timeIndex_)
{}


Foam::windkesselPressureFvPatchScalarField::
windkesselPressureFvPatchScalarField
(
    const windkesselPressureFvPatchScalarField& wkpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wkpsf, iF),
    po_(wkpsf.po_),
    poo_(wkpsf.poo_),
    Qo_(wkpsf.Qo_),
    Qoo_(wkpsf.Qoo_),
    R_(wkpsf.R_),
    Rch_(wkpsf.Rch_),
    C_(wkpsf.C_),
    rho_(wkpsf.rho_),
    timeIndex_(wkpsf.timeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::windkesselPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar dt = db().time().deltaTValue();

    // Volumetric flow rate through the patch
    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const scalar Qn = gSum(phip);
    const label timeIndex = db().time().timeIndex();

    if (timeIndex != timeIndex_)
    {
        timeIndex_ = timeIndex;

        poo_ = po_;
        po_ = gAverage(*this)*rho_;

        Qoo_ = Qo_;
        Qo_ = Qn;
    }

    scalar pn = po_;

    if (timeIndex > 1)
    {
        pn = (3*Qn + Qoo_ - 4*Qo_)*C_*R_*Rch_/(2*dt)
          + (Rch_ + R_)*Qn
          - (poo_ - 4*po_)*C_*R_/(2*dt);

        pn /= (1.0 + 3*C_*R_/(2*dt)) + SMALL;
    }
    else
    {
        pn = (Qn - Qo_)*C_*R_*Rch_/dt
          + (Rch_ + R_)*Qn
          + po_*C_*R_/dt;

        pn /= (1.0 + C_*R_/dt) + SMALL;
    }

    operator==(pn/(rho_ + SMALL));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::windkesselPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);

#ifdef OPENFOAM_COM
    os.writeEntry("Qo", Qo_);
    os.writeEntry("Qoo", Qoo_);
    os.writeEntry("po", po_);
    os.writeEntry("poo", poo_);

    os.writeEntryIfDifferent<scalar>("R", 1, R_);
    os.writeEntryIfDifferent<scalar>("Rch", 1, Rch_);
    os.writeEntryIfDifferent<scalar>("C", 1, C_);
    os.writeEntryIfDifferent<scalar>("rho", 1, rho_);
#else
    writeEntry(os, "Qo", Qo_);
    writeEntry(os, "Qoo", Qoo_);
    writeEntry(os, "po", po_);
    writeEntry(os, "poo", poo_);

    writeEntryIfDifferent<scalar>(os, "R",    1, R_);
    writeEntryIfDifferent<scalar>(os, "Rch",  1, Rch_);
    writeEntryIfDifferent<scalar>(os, "C",    1, C_);
    writeEntryIfDifferent<scalar>(os, "rho",  1, rho_);

    writeEntry(os, "value", *this);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        windkesselPressureFvPatchScalarField
    );
}

// ************************************************************************* //
