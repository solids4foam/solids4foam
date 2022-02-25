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

#include "outflowPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

outflowPressureFvPatchScalarField::outflowPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


outflowPressureFvPatchScalarField::outflowPressureFvPatchScalarField
(
    const outflowPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


outflowPressureFvPatchScalarField::outflowPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
//     fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

// #ifdef OPENFOAMESIORFOUNDATION
//     const fvMesh& mesh = internalField().mesh();
//     const pointField& points = mesh.points();
// #else
//     const fvMesh& mesh = dimensionedInternalField().mesh();
//     const pointField& points = mesh.allPoints();
// #endif

//     forAll(Fc_, i)
//     {
//         Fc_[i] = patch().patch()[i].centre(points);
//     }

//     oldFc_ = Fc_;
//     oldoldFc_ = Fc_;
}


#ifndef OPENFOAMFOUNDATION
outflowPressureFvPatchScalarField::outflowPressureFvPatchScalarField
(
    const outflowPressureFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf)
{}
#endif


outflowPressureFvPatchScalarField::outflowPressureFvPatchScalarField
(
    const outflowPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void outflowPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField pp(this->patch().size(), 0);

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    IOdictionary transportProperties
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
    dimensionedScalar nu(transportProperties.lookup("nu"));

    const fvPatchField<vector>& pU =
        patch().lookupPatchField<volVectorField, vector>("U");

    vectorField n(patch().nf());

    pp = nu.value()*(n & pU.snGrad());

    scalarField::operator=(pp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void outflowPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
//     writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    outflowPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
