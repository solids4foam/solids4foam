/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extrapolatedVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF)
{}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF)
{
    fvPatchField<vector>::operator=(patchInternalField());
}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const extrapolatedVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const extrapolatedVelocityFvPatchVectorField& ptf
)
:
    zeroGradientFvPatchVectorField(ptf)
{}


extrapolatedVelocityFvPatchVectorField::
extrapolatedVelocityFvPatchVectorField
(
    const extrapolatedVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedVelocityFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("grad(U)");

    vectorField delta = this->patch().delta();

    Field<vector>::operator=
    (
        this->patchInternalField() + (delta & gradU.patchInternalField())
    );

    fvPatchField<vector>::evaluate();
}


// void extrapolatedVelocityFvPatchVectorField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }

//     const uniformDimensionedVectorField& g =
//         db().lookupObject<uniformDimensionedVectorField>("g");

//     const fvPatchField<vector>& rho =
//         patch().lookupPatchField<volVectorField, vector>(rhoName_);

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

//     zeroGradientFvPatchVectorField::updateCoeffs();
// }


void extrapolatedVelocityFvPatchVectorField::write(Ostream& os) const
{
    zeroGradientFvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    extrapolatedVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
