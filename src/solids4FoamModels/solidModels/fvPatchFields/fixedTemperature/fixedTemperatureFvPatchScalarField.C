/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "fixedTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperatureSeries_()
{}


fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fixedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    temperatureSeries_(ptf.temperatureSeries_)
{}


fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    temperatureSeries_()
{
    // Check if temperature is time-varying
    if (dict.found("temperatureSeries"))
    {
        Info<< "temperature is time-varying" << endl;
        temperatureSeries_ =
            interpolationTable<scalar>(dict.subDict("temperatureSeries"));

        fvPatchField<scalar>::operator==
        (
            temperatureSeries_(this->db().time().timeOutputValue())
        );
    }
}

#ifndef OPENFOAM_ORG
fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fixedTemperatureFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    temperatureSeries_(pivpvf.temperatureSeries_)
{}
#endif

fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fixedTemperatureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    temperatureSeries_(pivpvf.temperatureSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (temperatureSeries_.size())
    {
        fvPatchField<scalar>::operator==
        (
            temperatureSeries_(db().time().timeOutputValue())
        );
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<Foam::Field<scalar> > fixedTemperatureFvPatchScalarField::snGrad() const
{
    // Lookup gradient of temperature from the solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
#ifdef OPENFOAM_NOT_EXTEND
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vectors
    const vectorField k((I - sqr(n)) & delta);

    // Correction sngrad
    return
    (
        *this - (patchInternalField() + (k & gradField.patchInternalField()))
    )*patch().deltaCoeffs();
}

tmp<Field<scalar> > fixedTemperatureFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    // Lookup gradoent of temperature from the solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
#ifdef OPENFOAM_NOT_EXTEND
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vectors
    const vectorField k((I - sqr(n)) & delta);

    return patch().deltaCoeffs()*(*this - (k & gradField.patchInternalField()));
}

void fixedTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
