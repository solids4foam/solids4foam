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

#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elasticSlipWallVelocityFvPatchVectorField::
elasticSlipWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
#ifdef OPENFOAMESIORFOUNDATION
    myTimeIndex_(internalField().mesh().time().timeIndex()),
#else
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero),
    movingWallVelocity_(p.patch().size(),vector::zero)
{}


elasticSlipWallVelocityFvPatchVectorField::
elasticSlipWallVelocityFvPatchVectorField
(
    const elasticSlipWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    myTimeIndex_(ptf.myTimeIndex_),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero),
    movingWallVelocity_(p.patch().size(),vector::zero)
{}


elasticSlipWallVelocityFvPatchVectorField::
elasticSlipWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
#ifdef OPENFOAMESIORFOUNDATION
    myTimeIndex_(internalField().mesh().time().timeIndex()),
#else
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero),
    movingWallVelocity_(p.patch().size(),vector::zero)
{
//     fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = internalField().mesh();
    const pointField& points = mesh.points();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const pointField& points = mesh.allPoints();
#endif

    forAll(Fc_, i)
    {
        Fc_[i] = patch().patch()[i].centre(points);
    }

    oldFc_ = Fc_;
    oldoldFc_ = Fc_;

    //
    if (dict.found("refValue"))
    {
        refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        refValue() = vector::zero;
    }

    if (dict.found("refGradient"))
    {
        refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        refGrad() = vector::zero;
    }

    // Patch normal
    vectorField n(patch().nf());

    valueFraction() = sqr(n);

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        const Field<vector> normalValue(transform(valueFraction(), refValue()));

        const Field<vector> gradValue
        (
            patchInternalField() + refGrad()/patch().deltaCoeffs()
        );

        const Field<vector> transformGradValue
        (
            transform(I - valueFraction(), gradValue)
        );

        Field<vector>::operator=(normalValue + transformGradValue);
    }
}

#ifndef OPENFOAMFOUNDATION
elasticSlipWallVelocityFvPatchVectorField::
elasticSlipWallVelocityFvPatchVectorField
(
    const elasticSlipWallVelocityFvPatchVectorField& pivpvf
)
:
    solidDirectionMixedFvPatchVectorField(pivpvf),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.Fc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_),
    movingWallVelocity_(pivpvf.movingWallVelocity_)
{}
#endif


elasticSlipWallVelocityFvPatchVectorField::
elasticSlipWallVelocityFvPatchVectorField
(
    const elasticSlipWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(pivpvf, iF),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.oldFc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_),
    movingWallVelocity_(pivpvf.movingWallVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticSlipWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    Info << "elasticSlipWallVelocityFvPatchVectorField::updateCoeffs" << endl;

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = internalField().mesh();
    const pointField& points = mesh.points();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const pointField& points = mesh.allPoints();
#endif
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();

    vectorField Up(p.size(), vector::zero);

    //const pointField& oldPoints = mesh.oldPoints();
    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
#ifdef OPENFOAMESIORFOUNDATION
            internalField().name()
#else
            dimensionedInternalField().name()
#endif
        );

    const word ddtScheme
    (
#ifdef OPENFOAMESIORFOUNDATION
        mesh.ddtScheme("ddt(" + U.name() +')')
#else
        mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
#endif
    );

    if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
    {
        if (myTimeIndex_ < mesh.time().timeIndex())
        {
            oldoldFc_ = oldFc_;
            oldFc_ = Fc_;

            forAll(Fc_, i)
            {
                Fc_[i] = pp[i].centre(points);
            }

            myTimeIndex_ = mesh.time().timeIndex();
        }

        const scalar deltaT = mesh.time().deltaT().value();
        scalar deltaT0 = mesh.time().deltaT0().value();

        if
        (
            U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex()
         || U.oldTime().oldTime().timeIndex() < 0
        )
        {
            deltaT0 = GREAT;
        }

        //Set coefficients based on deltaT and deltaT0
        const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));

        Up = coefft*(Fc_ - oldFc_)/deltaT
           - coefft00*(oldFc_ - oldoldFc_)/deltaT;
    }
    else // Euler
    {
        if (myTimeIndex_ < mesh.time().timeIndex())
        {
            oldFc_ = Fc_;

            forAll(Fc_, i)
            {
                Fc_[i] = pp[i].centre(points);
            }

            myTimeIndex_ = mesh.time().timeIndex();
        }

        Up = (Fc_ - oldFc_)/mesh.time().deltaT().value();
    }

    const scalarField phip
    (
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
    );

    const vectorField n(p.nf());
    const scalarField& magSf = p.magSf();
    scalarField Un(phip/(magSf + VSMALL));

    Up += n*(Un - (n & Up));

    movingWallVelocity_ = Up;

    valueFraction() = sqr(n);

    if (mesh.foundObject<surfaceScalarField>("phi"))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>("phi");

        const scalarField phip(phi.boundaryField()[patch().index()]);

        Un = phip/(magSf + VSMALL);
    }
    else
    {
        const volScalarField& pressure =
            mesh.lookupObject<volScalarField>("p");
        const scalarField nGradP
        (
            pressure.boundaryField()[patch().index()].snGrad()
        );

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

        const dimensionedScalar rho
        (
            transportProperties.lookup("rho")
        );

        const vectorField UOld
        (
            U.oldTime().boundaryField()[patch().index()]
        );

        Un = (n & UOld) - nGradP*mesh.time().deltaT().value()/rho.value();
    }

    refValue() = n*Un;

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}

void elasticSlipWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    elasticSlipWallVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
