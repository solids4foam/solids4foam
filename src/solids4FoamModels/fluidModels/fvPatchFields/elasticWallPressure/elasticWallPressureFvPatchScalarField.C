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

#include "elasticWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidSolidInterface.H"
#include "backwardDdtScheme.H"
#include "EulerDdtScheme.H"
#include "fvcMeshPhi.H"
// #include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(p, iF),
#ifdef OPENFOAMESIORFOUNDATION
    timeIndex_(internalField().mesh().time().timeIndex()),
#else
    timeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldOldFc_(p.patch().size(),vector::zero),
    U_(p.patch().size(),vector::zero),
    oldU_(p.patch().size(),vector::zero),
    oldOldU_(p.patch().size(),vector::zero)
{}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    robinFvPatchScalarField(ptf, p, iF, mapper),
    timeIndex_(ptf.timeIndex_),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldOldFc_(p.patch().size(),vector::zero),
    U_(p.patch().size(),vector::zero),
    oldU_(p.patch().size(),vector::zero),
    oldOldU_(p.patch().size(),vector::zero)
{}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    robinFvPatchScalarField(p, iF),
#ifdef OPENFOAMESIORFOUNDATION
    timeIndex_(internalField().mesh().time().timeIndex()),
#else
    timeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldOldFc_(p.patch().size(),vector::zero),
    U_(p.patch().size(),vector::zero),
    oldU_(p.patch().size(),vector::zero),
    oldOldU_(p.patch().size(),vector::zero)
{
    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }

    this->coeff0() = 1.0;
    this->coeff1() = 1.0;

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
    oldOldFc_ = Fc_;
}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf
)
:
    robinFvPatchScalarField(pivpvf),
    timeIndex_(pivpvf.timeIndex_),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_),
    Fc_(pivpvf.Fc_),
    oldFc_(pivpvf.oldFc_),
    oldOldFc_(pivpvf.oldOldFc_),
    U_(pivpvf.U_),
    oldU_(pivpvf.oldU_),
    oldOldU_(pivpvf.oldOldU_)
{}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(pivpvf, iF),
    timeIndex_(pivpvf.timeIndex_),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_),
    Fc_(pivpvf.oldFc_),
    oldFc_(pivpvf.oldFc_),
    oldOldFc_(pivpvf.oldOldFc_),
    U_(pivpvf.oldU_),
    oldU_(pivpvf.oldU_),
    oldOldU_(pivpvf.oldOldU_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticWallPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
}


void elasticWallPressureFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);
}

void elasticWallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    // Looking up fsi solver
    const fluidSolidInterface& fsi =
        mesh.objectRegistry::parent().lookupObject<fluidSolidInterface>
        (
            "fsiProperties"
        );

    // Solid properties
    // PC: hmnn what if the solidModel does not use mu and lambda...
    // It seems that ap is the speed of sound so we just need the stiffness, we
    // can lookup impK
    // Also, what happends if mu/lambda are varying...
    const volScalarField rho = fsi.solid().mechanical().rho()();

    if (rho.size() == 0)
    {
        return;
    }

    const scalar rhoSolid = rho[0];
    // scalar mu =
    //     fsi.solid().mechanical().mu()()[0];
    // scalar lambda =
    //     fsi.solid().mechanical().lambda()()[0];
    const scalar impK =
        fsi.solid().mechanical().impK()()[0];

    //const scalar ap = sqrt((lambda+2*mu)/rhoSolid);
    const scalar ap = sqrt(impK/rhoSolid);
    const scalar hs = ap*mesh.time().deltaT().value();

    // Fluid properties
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

    const dimensionedScalar rhoFluid
    (
        transportProperties.lookup("rho")
    );

    Info<< "rhoSolid = " << rhoSolid << ", hs = " << hs
        << ", rhoFluid = " << rhoFluid.value() << endl;

    // Update velocity and acceleration

    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();
#ifdef OPENFOAMESIORFOUNDATION
    const pointField& points = mesh.points();
#else
    const pointField& points = mesh.allPoints();
#endif

    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            "U"
        );

    scalarField phip =
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));
    vectorField n = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    const word ddtScheme
    (
#ifdef OPENFOAMESIORFOUNDATION
        mesh.ddtScheme("ddt(" + U.name() +')')
#else
        mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
#endif
    );

    if (ddtScheme == fv::EulerDdtScheme<vector>::typeName)
    {
        if (timeIndex_ < mesh.time().timeIndex())
        {
            oldOldFc_ = oldFc_;
            oldFc_ = Fc_;

            oldOldU_ = oldU_;
            oldU_ = U_;

            timeIndex_ = mesh.time().timeIndex();
        }

        forAll(Fc_, i)
        {
            Fc_[i] = pp[i].centre(points);
        }

        U_ = (Fc_ - oldFc_)/mesh.time().deltaT().value();
        U_ += n*(Un - (n & U_));

//         prevAcceleration_ =
//             (U_ - oldU_)/mesh.time().deltaT().value();
    }
    else
    {
        Info<< "elasticWallPressureFvPatchScalarField::updateCoeffs()"
            << endl;
    }

    // Info << "ddtUn, max: " << max(n&prevAcceleration_)
    //     << ", avg: " << average(n&prevAcceleration_)
    //     << ", min: " << min(n&prevAcceleration_) << endl;

    // Info << "p, max: " << max(prevPressure_/rhoFluid.value())
    //     << ", avg: " << average(prevPressure_/rhoFluid.value())
    //     << ", min: " << min(prevPressure_/rhoFluid.value()) << endl;




//     Info << rhoSolid_ << ", " << h_ << endl;

//     if(timeIndex_ < mesh.time().timeIndex())
//     {
//         timeIndex_ = mesh.time().timeIndex();
//     }

//     Info << ap << endl;

#ifdef OPENFOAMESIORFOUNDATION
    const word fieldName = internalField().name();
#else
    const word fieldName = dimensionedInternalField().name();
#endif

    const volScalarField& pressure =
        mesh.lookupObject<volScalarField>(fieldName);

    // const volVectorField& ddtU =
    //     mesh.lookupObject<volVectorField>("ddt(U)");

    // scalarField ddtUn =
    //     (n & ddtU.boundaryField()[patch().index()]);

    // Info << "ddtUn2, max: " << max(ddtUn)
    //     << ", avg: " << average(ddtUn)
    //     << ", min: " << min(ddtUn) << endl;

//     const volVectorField& U =
//         mesh.lookupObject<volVectorField>("U");

//     vectorField n = this->patch().nf();

//     scalarField prevDdtUn =
//         (n & fvc::ddt(U)().boundaryField()[patch().index()]);

    scalarField prevDdtUn = (n & prevAcceleration_);

    if (pressure.dimensions() == dimPressure/dimDensity)
    {
        // p/rho
        this->coeff0() = 1.0;
        this->coeff1() = rhoSolid*hs/rhoFluid.value();
        this->rhs() =
            prevPressure_/rhoFluid.value()
          - rhoSolid*hs*prevDdtUn/rhoFluid.value();
    }
    else
    {
        // p
        this->coeff0() = 1.0;
        this->coeff1() = rhoSolid*hs/rhoFluid.value();
        this->rhs() = prevPressure_ - rhoSolid*hs*prevDdtUn;
    }

//     if (weak_)
//     {
//         this->rhs() =
//             p.oldTime().boundaryField()[patch().index()]
//           + rhoSolid_*h_
//            *(n & ddtU.oldTime().boundaryField()[patch().index()]);
//     }
//     else
//     {
//        this->rhs() =
//            p.prevIter().boundaryField()[patch().index()]
//          - rhoSolid_*h_
//           *(n & ddtU.prevIter().boundaryField()[patch().index()]);
//     }

//     scalarField dn = 1.0/this->patch().deltaCoeffs();

//     Info << "pcoeff " << max(this->coeff()+dn)
//         << ", " << average(this->coeff()+dn) << endl;

    robinFvPatchField<scalar>::updateCoeffs();
}


void elasticWallPressureFvPatchScalarField::patchFlux
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<scalar>& matrix
) const
{
    const label patchI = this->patch().index();

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    const fvsPatchField<scalar>& rAU =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "rAUf"
        );

#ifdef OPENFOAMESIORFOUNDATION
    flux.boundaryFieldRef()[patchI] =
#else
    flux.boundaryField()[patchI] =
#endif
        rAU*this->snGrad()
       *mesh.magSf().boundaryField()[patchI];
}


void elasticWallPressureFvPatchScalarField::write(Ostream& os) const
{
    robinFvPatchScalarField::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    elasticWallPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
