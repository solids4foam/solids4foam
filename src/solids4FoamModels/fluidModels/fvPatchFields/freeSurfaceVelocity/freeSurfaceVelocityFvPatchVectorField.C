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

#include "freeSurfaceVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    fvPatchVectorField::operator=(patchInternalField());
}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const freeSurfaceVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::freeSurfaceVelocityFvPatchVectorField::
freeSurfaceVelocityFvPatchVectorField
(
    const freeSurfaceVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(fcvpvf, iF),
    phiName_(fcvpvf.phiName_),
    rhoName_(fcvpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfaceVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    {
        const vectorField n(patch().nf());

#ifdef OPENFOAMESIORFOUNDATION
        word UName = this->internalField().name();
#else
        word UName = this->dimensionedInternalField().name();
#endif
        const fvPatchField<tensor>& gradU =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + UName + ")"
            );

        const tensorField gradUP(gradU.patchInternalField());

        const vectorField sGradUn(((I-sqr(n)) & gradUP) & n);
        const vectorField nGradUn(-tr(((I-sqr(n)) & gradUP) & (I-sqr(n)))*n);

        gradient() = nGradUn - sGradUn;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::freeSurfaceVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

//     fixedGradientFvPatchVectorField::evaluate();

    // Evaluate tangential component of velocity
    // using second order discretisation
    if (true)
    {
//         Info << "Evaluate free surface velocity "
//             << "using second order correction"
//             << endl;

        const vectorField n(patch().nf());
        const vectorField delta(patch().delta());
        const vectorField k((I - sqr(n)) & delta);

#ifdef OPENFOAMESIORFOUNDATION
        word UName = this->internalField().name();
#else
        word UName = this->dimensionedInternalField().name();
#endif

        const fvPatchField<tensor>& gradU =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + UName + ")"
            );

        const vectorField dUP(k & gradU.patchInternalField());

        if (true)
        {
            const vectorField nGradUP(n & gradU.patchInternalField());

            Field<vector>::operator=
            (
                this->patchInternalField() + dUP
              + 0.5*(gradient() + nGradUP)
               /this->patch().deltaCoeffs()
            );
        }
        else
        {
            Field<vector>::operator=
            (
                this->patchInternalField() + dUP
              + gradient()/this->patch().deltaCoeffs()
            );
        }
    }

    if (!this->db().objectRegistry::found(phiName_))
    {
        // Flux not available, do not update
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    const vectorField n(patch().nf());
    const Field<scalar>& magS = patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(*this - n*(n & *this) + n*phip/magS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        if (!this->db().objectRegistry::found(rhoName_))
        {
            // Rho not available, do not update
            return;
        }

        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(*this - n*(n & *this) + n*phip/(rhop*magS));
    }
    else
    {
        FatalErrorIn
        (
            "freeSurfaceVelocityFvPatchVectorField::evaluate()"
        )
            << "dimensions of phi are incorrect\n"
            << "    on patch " << this->patch().name()
#ifdef OPENFOAMESIORFOUNDATION
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
#else
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
#endif
            << exit(FatalError);
    }

    fvPatchField<vector>::evaluate();
}


void Foam::freeSurfaceVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);

    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        freeSurfaceVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
