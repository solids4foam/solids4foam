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

#include "fixedValueCorrectedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "patchCorrectionVectors.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedValueCorrectedFvPatchScalarField::
fixedValueCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    nonOrthogonalCorrections_(true)
{}


Foam::fixedValueCorrectedFvPatchScalarField::
fixedValueCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    )
{}


Foam::fixedValueCorrectedFvPatchScalarField::
fixedValueCorrectedFvPatchScalarField
(
    const fixedValueCorrectedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    nonOrthogonalCorrections_(true)
{}


#ifndef OPENFOAM_ORG
Foam::fixedValueCorrectedFvPatchScalarField::
fixedValueCorrectedFvPatchScalarField
(
    const fixedValueCorrectedFvPatchScalarField& mwvpvf
)
:
    fixedValueFvPatchScalarField(mwvpvf)
{}
#endif

Foam::fixedValueCorrectedFvPatchScalarField::
fixedValueCorrectedFvPatchScalarField
(
    const fixedValueCorrectedFvPatchScalarField& mwvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar> >
Foam::fixedValueCorrectedFvPatchScalarField::snGrad() const
{
    const word gradName(internalFieldName(*this));

    if (nonOrthogonalCorrections_ && db().foundObject<volVectorField>(gradName))
    {
        // Lookup the gradient field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>(gradName);

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            *this
          - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        // fixedValue snGrad with no correction
        return (*this - patchInternalField())*patch().deltaCoeffs();
    }
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::fixedValueCorrectedFvPatchScalarField::gradientBoundaryCoeffs() const
{
    const word gradName(internalFieldName(*this));

    if (nonOrthogonalCorrections_ && db().foundObject<volVectorField>(gradName))
    {
        // Lookup the gradient field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>(gradName);

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            patch().deltaCoeffs()
           *(*this - (k & gradField.patchInternalField()))
        );
    }
    else
    {
        return (patch().deltaCoeffs()*(*this));
    }
}


void Foam::fixedValueCorrectedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("nonOrthogonalCorrections")
        << nonOrthogonalCorrections_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedValueCorrectedFvPatchScalarField
    );
}

// ************************************************************************* //
