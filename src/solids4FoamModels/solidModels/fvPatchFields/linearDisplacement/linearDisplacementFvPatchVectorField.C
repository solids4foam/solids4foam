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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "linearDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearDisplacementFvPatchVectorField::linearDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    rate_(vector::zero),
    initialC_(p.size(), vector::zero)
{}


linearDisplacementFvPatchVectorField::linearDisplacementFvPatchVectorField
(
    const linearDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    rate_(ptf.rate_),
    initialC_(ptf.initialC_)
{}


linearDisplacementFvPatchVectorField::linearDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    rate_(dict.lookup("rate")),
    initialC_(p.size(), vector::zero)
{
    if (dict.found("initialC"))
    {
        initialC_ = vectorField("initialC", dict, p.size());
    }
    else
    {
        initialC_ = patch().Cf();
    }
}

#ifndef OPENFOAMFOUNDATION
linearDisplacementFvPatchVectorField::linearDisplacementFvPatchVectorField
(
    const linearDisplacementFvPatchVectorField& pivpvf
)
:
    fixedDisplacementFvPatchVectorField(pivpvf),
    rate_(pivpvf.rate_),
    initialC_(pivpvf.initialC_)
{}
#endif

linearDisplacementFvPatchVectorField::linearDisplacementFvPatchVectorField
(
    const linearDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(pivpvf, iF),
    rate_(pivpvf.rate_),
    initialC_(pivpvf.initialC_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void linearDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    totalDisp() = db().time().value()*cmptMultiply(rate_, initialC_);

    fixedDisplacementFvPatchVectorField::updateCoeffs();
}


void linearDisplacementFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("rate")
        << rate_ << token::END_STATEMENT << nl;
    initialC_.writeEntry("initialC", os);

    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    linearDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
