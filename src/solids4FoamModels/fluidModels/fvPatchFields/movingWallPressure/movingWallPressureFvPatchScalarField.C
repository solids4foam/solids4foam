/*---------------------------------------------------------------------------*\
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

#include "movingWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "newMovingWallVelocityFvPatchVectorField.H"
#include "wallNormalMomentum.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<vectorField> movingWallPressureFvPatchScalarField::wallAcceleration() const
{
    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    if (isA<newMovingWallVelocityFvPatchVectorField>(Up))
    {
        return tmp<vectorField>
        (
            new vectorField
            (
                refCast<const newMovingWallVelocityFvPatchVectorField>
                (
                    Up
                ).acceleration()
            )
        );
    }

    // Time derivative of the boundary velocity with the ddt scheme of U
    const volVectorField& U = db().lookupObject<volVectorField>("U");

    return tmp<vectorField>
    (
        new vectorField(fvc::ddt(U)().boundaryField()[patch().index()])
    );
}


tmp<scalarField> movingWallPressureFvPatchScalarField::patchDensity() const
{
    if (db().foundObject<volScalarField>("rho"))
    {
        return tmp<scalarField>
        (
            new scalarField
            (
                patch().lookupPatchField<volScalarField, scalar>("rho")
            )
        );
    }

    // The transport properties of the fluid model, or else those on disk
    autoPtr<IOdictionary> readDictPtr;

    if (!db().foundObject<IOdictionary>("transportProperties"))
    {
        readDictPtr.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "transportProperties",
                    db().time().constant(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
    }

    const dictionary& transportProperties =
        readDictPtr.valid()
      ? readDictPtr()
      : db().lookupObject<IOdictionary>("transportProperties");

    const dimensionedScalar rho(transportProperties.lookup("rho"));

    return tmp<scalarField>(new scalarField(patch().size(), rho.value()));
}


tmp<scalarField> movingWallPressureFvPatchScalarField::momentumGradient
(
    const volVectorField* viscousPtr
) const
{
    const vectorField n(patch().nf());

    vectorField momentum(-wallAcceleration());

    if (viscousTerm_)
    {
        momentum +=
            viscousPtr->boundaryField()[patch().index()].patchInternalField();
    }

    tmp<scalarField> tgrad(new scalarField(n & momentum));

#ifdef OPENFOAM_NOT_EXTEND
    const word& fieldName = internalField().name();
#else
    const word& fieldName = dimensionedInternalField().name();
#endif

    const volScalarField& p = db().lookupObject<volScalarField>(fieldName);

    if (p.dimensions() == dimPressure)
    {
        tmpRef(tgrad) *= patchDensity();
    }

    return tgrad;
}


tmp<scalarField> movingWallPressureFvPatchScalarField::meanWeights() const
{
    tmp<scalarField> tw(new scalarField(patch().magSf()));

    if (db().foundObject<volScalarField>("rAU"))
    {
        tmpRef(tw) *= patch().lookupPatchField<volScalarField, scalar>("rAU");
    }

    return tw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    viscousTerm_(true)
{}


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    viscousTerm_(dict.lookupOrDefault<Switch>("viscousTerm", true))
{
    fvPatchField<scalar>::operator=(patchInternalField());
}


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const movingWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    viscousTerm_(ptf.viscousTerm_)
{}


#ifndef OPENFOAM_ORG
movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const movingWallPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    viscousTerm_(ptf.viscousTerm_)
{}
#endif


movingWallPressureFvPatchScalarField::
movingWallPressureFvPatchScalarField
(
    const movingWallPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    viscousTerm_(ptf.viscousTerm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingWallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#ifdef OPENFOAM_NOT_EXTEND
    const word& fieldName = internalField().name();
#else
    const word& fieldName = dimensionedInternalField().name();
#endif
    const volScalarField& p = db().lookupObject<volScalarField>(fieldName);

    // Enclosed domain: the mean is taken over all movingWallPressure patches
    // of the field, so all of them need their gradient here. The decisions
    // depend on dictionary entries and reduced values only, so all
    // processors take the same ones and call the same collectives
    const bool enclosed = p.needReference();

    bool needViscous = viscousTerm_;

    if (enclosed)
    {
        forAll(p.boundaryField(), patchI)
        {
            if
            (
                isA<movingWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
            )
            {
                needViscous =
                    needViscous
                 || refCast<const movingWallPressureFvPatchScalarField>
                    (
                        p.boundaryField()[patchI]
                    ).viscousTerm_;
            }
        }
    }

    autoPtr<volVectorField> viscousFieldPtr;

    if (needViscous)
    {
        viscousFieldPtr.reset
        (
            wallNormalMomentum::field(db(), "U", patch().name()).ptr()
        );
    }

    const volVectorField* viscousPtr =
        needViscous ? &viscousFieldPtr() : nullptr;

    gradient() = momentumGradient(viscousPtr);

    if (enclosed)
    {
        scalar sumWeightedGrad = 0;
        scalar sumWeights = 0;

        forAll(p.boundaryField(), patchI)
        {
            if
            (
                isA<movingWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
            )
            {
                const movingWallPressureFvPatchScalarField& pp =
                    refCast<const movingWallPressureFvPatchScalarField>
                    (
                        p.boundaryField()[patchI]
                    );

                const scalarField w(pp.meanWeights());
                const scalarField g
                (
                    patchI == patch().index()
                  ? gradient()
                  : pp.momentumGradient(viscousPtr)()
                );

                sumWeightedGrad += sum(w*g);
                sumWeights += sum(w);
            }
        }

        reduce(sumWeightedGrad, sumOp<scalar>());
        reduce(sumWeights, sumOp<scalar>());

        if (sumWeights > VSMALL)
        {
            gradient() -= sumWeightedGrad/sumWeights;
        }
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void movingWallPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("viscousTerm")
        << viscousTerm_ << token::END_STATEMENT << nl;
#ifdef OPENFOAM_ORG
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    movingWallPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
