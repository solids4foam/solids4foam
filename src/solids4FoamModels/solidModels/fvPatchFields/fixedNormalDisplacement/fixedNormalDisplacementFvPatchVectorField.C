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
    fixedNormalDisplacementFvPatchVectorField

Description

\*---------------------------------------------------------------------------*/

#include "fixedNormalDisplacementFvPatchVectorField.H"
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

fixedNormalDisplacementFvPatchVectorField
::fixedNormalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF)
{}


fixedNormalDisplacementFvPatchVectorField
::fixedNormalDisplacementFvPatchVectorField
(
    const fixedNormalDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, p, iF, mapper)
{}


fixedNormalDisplacementFvPatchVectorField
::fixedNormalDisplacementFvPatchVectorField
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

//     if (dict.found("refValue"))
//     {
//         this->refValue() = vectorField("refValue", dict, p.size());
//     }
//     else
//     {
//         this->refValue() = vector::zero;
//     }

    if (dict.found("refGradient"))
    {
        this->refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = vector::zero;
    }

    // Patch normal
    vectorField n = patch().nf();

    this->valueFraction() = sqr(n);

//     if (dict.found("valueFraction"))
//     {
//         this->valueFraction() =
//             symmTensorField("valueFraction", dict, p.size());
//     }
//     else
//     {
//         this->valueFraction() = sqr(n);
//     }

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


fixedNormalDisplacementFvPatchVectorField
::fixedNormalDisplacementFvPatchVectorField
(
    const fixedNormalDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedNormalDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedDisplacementFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedNormalDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedDisplacementFvPatchVectorField::rmap(ptf, addr);

//     const fixedNormalDisplacementFvPatchVectorField& dmptf =
//         refCast<const fixedNormalDisplacementFvPatchVectorField>(ptf);
}


void fixedNormalDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Patch unit normals
    const vectorField n = patch().nf();

    // Set normal component to fixedValue
    valueFraction() = sqr(n);

    // Set tangential gradient to zero
    refGrad() = vector::zero;

    directionMixedDisplacementFvPatchVectorField::updateCoeffs();
}


// Write
void fixedNormalDisplacementFvPatchVectorField::write(Ostream& os) const
{
    directionMixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedNormalDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
