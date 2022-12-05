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

#include "freeEdgeDisplacementFaPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "kirchhoffPlateSolid.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(1.0)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))

{
   if (dict.found("value"))
   {
       faPatchField<scalar>::operator==(Field<scalar>("value", dict, p.size()));
   }
   else
   {
       faPatchField<scalar>::operator==(0.0);
   }

   gradient() = 0.0;
}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf,
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedGradientFaPatchField<scalar>(ptf, p, iF, mapper),
    relaxFac_(ptf.relaxFac_)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf
)
:
    fixedGradientFaPatchField<scalar>(ptf),
    relaxFac_(ptf.relaxFac_)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(ptf, iF),
    relaxFac_(ptf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeEdgeDisplacementFaPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Info<< nl << "------------------------------------------" << nl << endl;

    // Lookup angle of rotation field
    // const faPatchField<vector>& theta =
    //     patch().lookupPatchField<areaVectorField, vector>("theta");

    // Lookup the gradient of rotation field
    // const faPatchField<tensor>& gradTheta =
    //     patch().lookupPatchField<areaTensorField, tensor>("grad(theta)");

    // Calculate correction vectors
    // const vectorField n = patch().edgeNormals();
    // const vectorField delta(patch().delta());
    // const vectorField k = (I - sqr(n)) & delta;

    // Calculate the patch internal field and correction for non-orthogonality
    // const scalarField nDotThetaPif =
    //  n & (theta.patchInternalField() + (k & gradTheta.patchInternalField()));

    // Info<< "theta: "
    //     << vectorField(theta) << nl << endl;
    // Info<< "theta pif: "
    //     << (theta.patchInternalField()) << nl << endl;

    // gradTheta patch internal field in normal direction
    // const vectorField nDotGradThetaPif = n & gradTheta.patchInternalField();

    // Lookup moment sum field
    const faPatchField<scalar>& M =
        patch().lookupPatchField<areaScalarField, scalar>("M");

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
    const scalar nu = plateSolid.nu().value();

    // Update the gradient
    gradient() +=
        relaxFac_*(1.0 + nu)*M/(D*(1.0 - pow(nu, 2))*patch().deltaCoeffs());

    // Info<< "prev gradient(): " << gradient() << nl << endl;

    // gradient() =
    //     relaxFac_
    //    *(
    //        nDotThetaPif
    //      + (delta & nDotGradThetaPif)
    //      + M/(D*(1.0 - pow(nu, 2))*patch().deltaCoeffs())
    //     )
    //   + (1.0 - relaxFac_)*gradient();

    // Info<< "new gradient(): " << gradient() << nl << endl;

    // Info<< "this (w): " << scalarField(*this) << nl << endl;

    // Info<< "nDotThetaPif: " << nDotThetaPif << endl;
    // Info<< "(delta & nDotGradThetaPif) "
    //     << (delta & nDotGradThetaPif) << endl;
    // Info<< "3: "
    //     << M/(D*(1.0 - pow(nu, 2))*patch().deltaCoeffs()) << endl;

    fixedGradientFaPatchField<scalar>::updateCoeffs();

    //Info<< nl << "-------------------------------------------" << nl << endl;
}


void Foam::freeEdgeDisplacementFaPatchScalarField::write
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
        freeEdgeDisplacementFaPatchScalarField
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
