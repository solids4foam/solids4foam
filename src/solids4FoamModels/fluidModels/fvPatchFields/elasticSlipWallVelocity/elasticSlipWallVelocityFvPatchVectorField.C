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
        this->refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = vector::zero;
    }

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

    word ddtScheme
    (
#ifdef OPENFOAMESIORFOUNDATION
        mesh.ddtScheme("ddt(" + U.name() +')')
#else
        mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
#endif
    );

    if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
    {
        if(myTimeIndex_ < mesh.time().timeIndex())
        {
            oldoldFc_ = oldFc_;
            oldFc_ = Fc_;
//             Fc_ = pp.faceCentres();

            forAll(Fc_, i)
            {
                Fc_[i] = pp[i].centre(points);
            }

            myTimeIndex_ = mesh.time().timeIndex();
        }

        scalar deltaT = mesh.time().deltaT().value();
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
        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
//         scalar coefft0  = coefft + coefft00;

//         Up = (coefft*Fc_ - coefft0*oldFc_ + coefft00*oldoldFc_)
//            /mesh.time().deltaT().value();


        Up = coefft*(Fc_ - oldFc_)/deltaT
          - coefft00*(oldFc_ - oldoldFc_)/deltaT;

//         Info << max(mag(Up)) << endl;
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

//             Fc_ = pp.faceCentres();
            myTimeIndex_ = mesh.time().timeIndex();
        }

        Up = (Fc_ - oldFc_)/mesh.time().deltaT().value();
    }

    scalarField phip =
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

    vectorField n = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    Up += n*(Un - (n & Up));


//     vectorField::operator=(Up);
    movingWallVelocity_ = Up;

    this->valueFraction() = sqr(n);


//     this->refValue() = movingWallVelocity_;
//     Info << max(Un) << endl;

    if (mesh.foundObject<surfaceScalarField>("phi"))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>("phi");

        scalarField phip = phi.boundaryField()[this->patch().index()];
        Un = phip/(magSf + VSMALL);
    }
    else
    {
        const volScalarField& pressure =
            mesh.lookupObject<volScalarField>("p");
        scalarField nGradP =
            pressure.boundaryField()[this->patch().index()].snGrad();

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

        dimensionedScalar rho
        (
            transportProperties.lookup("rho")
        );

        vectorField UOld = U.oldTime().boundaryField()[this->patch().index()];
        Un = (n & UOld) - nGradP*mesh.time().deltaT().value()/rho.value();
    }

    this->refValue() = n*Un;

//     Info << max(mag(nGradP)) << endl;

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
