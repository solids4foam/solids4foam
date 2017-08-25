/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_()
{
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

    // Check if traction is time-varying
    if (dict.found("tractionSeries"))
    {
        Info<< "traction is time-varying" << endl;
        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        Info<< "pressure is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    traction_(stpvf.traction_, mapper),
    pressure_(stpvf.pressure_, mapper),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const solidTractionFvPatchVectorField& dmptf =
        refCast<const solidTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    // Lookup the solidModel object
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const solidModel* solModPtr = NULL;
    if (mesh.foundObject<solidModel>("solidProperties"))
    {
        solModPtr = &mesh.lookupObject<solidModel>("solidProperties");
    }
    else
    {
        solModPtr = &mesh.parent().lookupObject<solidModel>("solidProperties");
    }

    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    gradient() =
        solModPtr->tractionBoundarySnGrad(traction_, pressure_, patch());

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    vectorField n = patch().nf();
    vectorField delta = patch().delta();

    //- non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k&gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<vector>::evaluate();
}

// Write
void solidTractionFvPatchVectorField::write(Ostream& os) const
{
    //fvPatchVectorField::write(os);
    fixedGradientFvPatchVectorField::write(os);

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        traction_.writeEntry("traction", os);
    }

    if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        pressure_.writeEntry("pressure", os);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
