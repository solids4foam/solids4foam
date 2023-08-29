/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "extrapolatedPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedPressureFvPatchScalarField::
extrapolatedPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF)
{}


extrapolatedPressureFvPatchScalarField::
extrapolatedPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(patchInternalField());
}


extrapolatedPressureFvPatchScalarField::
extrapolatedPressureFvPatchScalarField
(
    const extrapolatedPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


#ifndef OPENFOAMFOUNDATION
extrapolatedPressureFvPatchScalarField::
extrapolatedPressureFvPatchScalarField
(
    const extrapolatedPressureFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf)
{}
#endif


extrapolatedPressureFvPatchScalarField::
extrapolatedPressureFvPatchScalarField
(
    const extrapolatedPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedPressureFvPatchScalarField::evaluate
(
#ifdef OPENFOAMESIORFOUNDATION
    const Pstream::commsTypes commsType
#else
    const Pstream::commsTypes commsType
#endif
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<vector>& gradP =
        patch().lookupPatchField<volVectorField, vector>("grad(p)");

    const vectorField delta(this->patch().delta());

    Field<scalar>::operator=
    (
        this->patchInternalField() + (delta & gradP.patchInternalField())
    );

    fvPatchField<scalar>::evaluate();
}


// void extrapolatedPressureFvPatchScalarField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }

//     const uniformDimensionedVectorField& g =
//         db().lookupObject<uniformDimensionedVectorField>("g");

//     const fvPatchField<scalar>& rho =
//         patch().lookupPatchField<volScalarField, scalar>(rhoName_);

//     // If the variable name is "p_rgh" or "pd" assume it is p - rho*g.h
//     // and set the gradient appropriately.
//     // Otherwise assume the variable is the static pressure.
//     if
//     (
// #ifdef OPENFOAMESIORFOUNDATION
//         internalField().name() == "p_rgh"
//      || internalField().name() == "pd"
// #else
//         dimensionedInternalField().name() == "p_rgh"
//      || dimensionedInternalField().name() == "pd"
// #endif
//     )
//     {
//         gradient() = -rho.snGrad()*(g.value() & patch().Cf());
//     }
//     else
//     {
//         gradient() = rho*(g.value() & patch().nf());
//     }

//     zeroGradientFvPatchScalarField::updateCoeffs();
// }


void extrapolatedPressureFvPatchScalarField::write(Ostream& os) const
{
    zeroGradientFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    extrapolatedPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
