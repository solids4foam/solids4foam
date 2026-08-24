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

#include "pressureTractionFvPatchScalarField.H"
#include "tractionPressureDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "lookupSolidModel.H"
#include "coupledPressureDisplacementSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureTractionFvPatchScalarField::
pressureTractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


pressureTractionFvPatchScalarField::
pressureTractionFvPatchScalarField
(
    const pressureTractionFvPatchScalarField& tdpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(tdpvf, p, iF, mapper)
{}


pressureTractionFvPatchScalarField::
pressureTractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator==(patchInternalField());
}


pressureTractionFvPatchScalarField::
pressureTractionFvPatchScalarField
(
    const pressureTractionFvPatchScalarField& tdpvf
)
:
    fixedValueFvPatchScalarField(tdpvf)
{}


pressureTractionFvPatchScalarField::
pressureTractionFvPatchScalarField
(
    const pressureTractionFvPatchScalarField& tdpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tdpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureTractionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void pressureTractionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    // const pressureTractionFvPatchScalarField& dmptf =
    //     refCast<const pressureTractionFvPatchScalarField>(ptf);
}


// Return gradient at boundary
Foam::tmp<Foam::Field<scalar>>
pressureTractionFvPatchScalarField::snGrad() const
{
    // // Face unit normals
    // const vectorField n = this->patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Uncorrected patch internal field
    scalarField DpP = patchInternalField(); // ZT: for uncorrected snGrad

    // This is related to satisfaction of continuity equation
    return (*this - DpP)/mag(delta);
}


// Update the coefficients associated with the patch field
void pressureTractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the solidModel object
    const solidModel& solMod =
        lookupSolidModel(patch().boundaryMesh().mesh());

    scalarField patchPressure(this->size(), 0);

    const solidModels::coupledPressureDisplacementSolid& mixedSolMod =
        refCast<const solidModels::coupledPressureDisplacementSolid>
        (
            solMod
        );

    const fvPatchField<vector>& pDD =
        patch().lookupPatchField<volVectorField, vector>("DD");

    const tractionPressureDisplacementFvPatchVectorField& tedDD =
        refCast
        <
            const tractionPressureDisplacementFvPatchVectorField
        >(pDD);

    patchPressure =
        mixedSolMod.tractionBoundaryHydPressure
        (
            tedDD.traction(),
            tedDD.pressure(),
            patch()
        );

    this->operator==(patchPressure);

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<Field<scalar> >
pressureTractionFvPatchScalarField::gradientBoundaryCoeffs() const
{
    // Delta vectors
    vectorField delta = patch().delta();

    return
    (
        // ZT: without correction since uncorrected snGrad in
        // Laplacian is use for pressure equation
        (*this)/mag(delta)
    );
}


tmp<Field<scalar> >
pressureTractionFvPatchScalarField::
gradientInternalCoeffs() const
{
    // Delta vectors
    vectorField delta = patch().delta();

    return -(1.0/mag(delta));
}


tmp<scalarField> pressureTractionFvPatchScalarField::normalTraction() const
{
    tmp<scalarField> tNormalTraction
    (
        new scalarField(this->size(), 0)
    );

    vectorField n = patch().nf();

    const fvPatchField<vector>& pDD =
        patch().lookupPatchField<volVectorField, vector>("DD");

    const tractionPressureDisplacementFvPatchVectorField& tpdDD =
        refCast
        <
            const tractionPressureDisplacementFvPatchVectorField
        >(pDD);


    tNormalTraction() = ((tpdDD.traction() - tpdDD.pressure()*n) & n);

    return tNormalTraction;
}


void pressureTractionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, pressureTractionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
