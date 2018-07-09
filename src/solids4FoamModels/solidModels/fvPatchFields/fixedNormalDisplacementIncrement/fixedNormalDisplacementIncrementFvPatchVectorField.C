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

Class
    fixedNormalDisplacementIncrementFvPatchVectorField

Description

\*---------------------------------------------------------------------------*/

#include "fixedNormalDisplacementIncrementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedNormalDisplacementIncrementFvPatchVectorField
::fixedNormalDisplacementIncrementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF)
{}


fixedNormalDisplacementIncrementFvPatchVectorField
::fixedNormalDisplacementIncrementFvPatchVectorField
(
    const fixedNormalDisplacementIncrementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, p, iF, mapper)
{}


fixedNormalDisplacementIncrementFvPatchVectorField
::fixedNormalDisplacementIncrementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
//     directionMixedDisplacementFvPatchVectorField(p, iF, dict)
    directionMixedDisplacementFvPatchVectorField(p, iF)
{
    this->refValue() = vectorField("refValue", dict, p.size());

    if (dict.found("refGradient"))
    {
        this->refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = vector::zero;
    }

    // Set normal direction to fixed value
    this->valueFraction() = sqr(patch().nf());

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector> normalValue = transform(valueFraction(), refValue());

        Field<vector> gradValue =
            this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

        Field<vector> transformGradValue =
            transform(I - valueFraction(), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }
}


fixedNormalDisplacementIncrementFvPatchVectorField
::fixedNormalDisplacementIncrementFvPatchVectorField
(
    const fixedNormalDisplacementIncrementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedNormalDisplacementIncrementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedDisplacementFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedNormalDisplacementIncrementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedDisplacementFvPatchVectorField::rmap(ptf, addr);

    // const fixedNormalDisplacementIncrementFvPatchVectorField& dmptf =
    // refCast<const fixedNormalDisplacementIncrementFvPatchVectorField>(ptf);
}


void fixedNormalDisplacementIncrementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Set tangential grdaient to zzero
    refGrad() = vector::zero;

    directionMixedDisplacementFvPatchVectorField::updateCoeffs();
}


// Write
void fixedNormalDisplacementIncrementFvPatchVectorField
::write(Ostream& os) const
{
    directionMixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedNormalDisplacementIncrementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
