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

#include "freeSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    curTimeIndex_(-1)
{
    fvPatchScalarField::operator=(patchInternalField());
}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    curTimeIndex_(-1)
{}


Foam::freeSurfacePressureFvPatchScalarField::
freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (true)
//     if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        const fvPatchField<vector>& U =
            patch().lookupPatchField<volVectorField, vector>
            (
                "U"
            );

        vectorField n = this->patch().nf();

        scalarField nGradUn = (n & U.snGrad());

        const fvMesh& mesh =
            this->patch().boundaryMesh().mesh();

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

//         scalarField newFsPressure = -2*nu.value()*nGradUn;
//         scalarField oldFsPressure = *this;
//         scalar relaxationFactor = 0.3;

//         fvPatchScalarField::operator=
//         (
//             oldFsPressure
//           + relaxationFactor*(newFsPressure - oldFsPressure)
//         );

        fvPatchScalarField::operator==(2*nu.value()*nGradUn);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::freeSurfacePressureFvPatchScalarField::snGrad() const
{
    word pName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradp =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + pName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    bool secondOrder_ = true;
    if (secondOrder_)
    {
        scalarField dpP = (k&gradp.patchInternalField());
        scalarField nGradpP = (n&gradp.patchInternalField());

        return
            2
           *(
                *this
              - (this->patchInternalField() + dpP)
            )*this->patch().deltaCoeffs()
          - nGradpP;
    }

    return
    (
        *this
      - (patchInternalField() + (k&gradp.patchInternalField()))
    )*this->patch().deltaCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::freeSurfacePressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    word pName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradp =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + pName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    bool secondOrder_ = true;
    if (secondOrder_)
    {
        scalarField dpP = (k&gradp.patchInternalField());
        scalarField nGradpP = (n&gradp.patchInternalField());

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dpP)
              - this->patchInternalField()
            )
          - nGradpP;
    }

    return this->patch().deltaCoeffs()
       *(*this - (k&gradp.patchInternalField()));
}


void Foam::freeSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
//     writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freeSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
