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

#include "mmsFsiVelocityGradientFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mmsFsiFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mmsFsiVelocityGradientFvPatchVectorField::
mmsFsiVelocityGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


Foam::mmsFsiVelocityGradientFvPatchVectorField::
mmsFsiVelocityGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
    // The gradient is set by updateCoeffs: start from zero gradient
    gradient() = Zero;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


Foam::mmsFsiVelocityGradientFvPatchVectorField::
mmsFsiVelocityGradientFvPatchVectorField
(
    const mmsFsiVelocityGradientFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::mmsFsiVelocityGradientFvPatchVectorField::
mmsFsiVelocityGradientFvPatchVectorField
(
    const mmsFsiVelocityGradientFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf)
{}


Foam::mmsFsiVelocityGradientFvPatchVectorField::
mmsFsiVelocityGradientFvPatchVectorField
(
    const mmsFsiVelocityGradientFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mmsFsiVelocityGradientFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();

    const vectorField n(patch().nf());

    forAll(Cf, faceI)
    {
        gradient()[faceI] =
            n[faceI] & mmsFsiFields::fluid(Cf[faceI], t).gradU;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        mmsFsiVelocityGradientFvPatchVectorField
    );
}


// ************************************************************************* //
