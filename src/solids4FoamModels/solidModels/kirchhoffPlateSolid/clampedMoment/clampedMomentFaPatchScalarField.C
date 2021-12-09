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

#include "clampedMomentFaPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "kirchhoffPlateSolid.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clampedMomentFaPatchScalarField::
clampedMomentFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedValueFaPatchField<scalar>(p, iF),
    relaxFac_(1.0)
{}


Foam::clampedMomentFaPatchScalarField::
clampedMomentFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<scalar>(p, iF),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 0.1))

{
   if (dict.found("value"))
   {
       faPatchField<scalar>::operator==(Field<scalar>("value", dict, p.size()));
   }
   else
   {
       faPatchField<scalar>::operator==(0.0);
   }
}


Foam::clampedMomentFaPatchScalarField::
clampedMomentFaPatchScalarField
(
    const clampedMomentFaPatchScalarField& ptf,
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<scalar>(ptf, p, iF, mapper),
    relaxFac_(ptf.relaxFac_)
{}


Foam::clampedMomentFaPatchScalarField::
clampedMomentFaPatchScalarField
(
    const clampedMomentFaPatchScalarField& ptf
)
:
    fixedValueFaPatchField<scalar>(ptf),
    relaxFac_(ptf.relaxFac_)
{}


Foam::clampedMomentFaPatchScalarField::
clampedMomentFaPatchScalarField
(
    const clampedMomentFaPatchScalarField& ptf,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedValueFaPatchField<scalar>(ptf, iF),
    relaxFac_(ptf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::clampedMomentFaPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Lookup angle of rotation field
    const faPatchField<vector>& theta =
        patch().lookupPatchField<areaVectorField, vector>("theta");

    // Lookup the gradient of rotation field
    const faPatchField<tensor>& gradTheta =
        patch().lookupPatchField<areaTensorField, tensor>("grad(theta)");

    // Calculate correction vectors
    const vectorField n = patch().edgeNormals();
    const vectorField delta(patch().delta());
    const vectorField k = (I - sqr(n)) & delta;

    // Calculate the patch internal field and correction for non-orthogonality
    const scalarField nDotThetaPif =
        n & (theta.patchInternalField() + (k & gradTheta.patchInternalField()));

    // Lookup fvMesh
    // Fix for FSI: this is only correct for
    const fvMesh* vmeshPtr = NULL;
    if (db().parent().foundObject<fvMesh>("solid"))
    {
        vmeshPtr = &db().parent().lookupObject<fvMesh>("solid");
    }
    else
    {
        vmeshPtr = &db().parent().lookupObject<fvMesh>("region0");
    }
    const fvMesh& vmesh = *vmeshPtr;

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(vmesh);

    // Cast the solid model to a Kirchhoff plate solid model
    const solidModels::kirchhoffPlateSolid& plateSolid =
        refCast<const solidModels::kirchhoffPlateSolid>(solMod);

    // Lookup flexural stiffness from the solver
    const scalar D = plateSolid.bendingStiffness().value();

    // Set the boundary moment sum to force a clamped edge
    const scalarField prevMomentSum = *this;
    faPatchField<scalar>::operator==
    (
        relaxFac_*(-D*nDotThetaPif*patch().deltaCoeffs())
      + (1.0 - relaxFac_)*prevMomentSum
    );

    fixedValueFaPatchField<scalar>::updateCoeffs();
}


void Foam::clampedMomentFaPatchScalarField::write
(
    Ostream& os
) const
{
    faPatchField<scalar>::write(os);

    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << endl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeFaPatchTypeField
    (
        faPatchScalarField,
        clampedMomentFaPatchScalarField
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
