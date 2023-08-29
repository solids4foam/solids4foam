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

#include "fixedVelocityPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedVelocityPressureFvPatchScalarField::
fixedVelocityPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


fixedVelocityPressureFvPatchScalarField::
fixedVelocityPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=(patchInternalField());
}


fixedVelocityPressureFvPatchScalarField::
fixedVelocityPressureFvPatchScalarField
(
    const fixedVelocityPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


#ifndef OPENFOAMFOUNDATION
fixedVelocityPressureFvPatchScalarField::
fixedVelocityPressureFvPatchScalarField
(
    const fixedVelocityPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
{}
#endif


fixedVelocityPressureFvPatchScalarField::
fixedVelocityPressureFvPatchScalarField
(
    const fixedVelocityPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void fixedVelocityPressureFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes
// )
// {
//     if (!this->updated())
//     {
//         this->updateCoeffs();
//     }

//     const fvMesh& mesh = this->patch().boundaryMesh().mesh();

//     if
//     (
//         mesh.found("HU")
//      && mesh.found("rAU")
//      && mesh.found("phi")
//     )
//     {
//         Info << "fixedVelocityPressureFvPatchScalarField::evaluate" << endl;

//         const fvPatchField<vector>& HU =
//             patch().lookupPatchField<volVectorField, vector>("HU");

//         const fvPatchField<scalar>& rAU =
//             patch().lookupPatchField<volScalarField, scalar>("rAU");

//         const fvsPatchField<scalar>& phi =
//             patch().lookupPatchField<surfaceScalarField, scalar>("phi");

//         scalarField magS = this->patch().magSf();
//         vectorField n(patch().nf());
//         scalarField dn = 1.0/this->patch().deltaCoeffs();

//         scalarField nGradP = (n & HU) - phi/rAU/magS;

//         Field<scalar>::operator=
//         (
//             this->patchInternalField() + (dn*nGradP)
//         );

//         fvPatchField<scalar>::evaluate();
//     }
//     else
//     {
//         fixedGradientFvPatchField<scalar>::evaluate();
//     }
// }


void fixedVelocityPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//     const fvPatchField<vector>& U =
//         patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& laplacianU =
        patch().lookupPatchField<volVectorField, vector>("laplacian(U)");

//     const fvPatchField<tensor>& gradU =
//         patch().lookupPatchField<volTensorField, tensor>("grad(U)");

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar nu(transportProperties.lookup("nu"));

    const scalarField dn(1.0/this->patch().deltaCoeffs());

    const vectorField n(patch().nf());

    gradient() = (n & laplacianU.patchInternalField());

    Info<< patch().name() << ", nGradP, max: " << max(gradient())
        << ", avg: " << average(gradient())
        << ", min: " << min(gradient()) << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fixedVelocityPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
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
    fixedVelocityPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
