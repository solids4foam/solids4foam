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

#include "elasticWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidSolidInterface.H"

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
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero)
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
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero)
{}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    robinFvPatchScalarField(p, iF),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero)
{
    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }

    this->coeff0() = 1.0;
    this->coeff1() = 1.0;
}


#ifndef OPENFOAM_ORG
elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf
)
:
    robinFvPatchScalarField(pivpvf),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_)
{}
#endif


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(pivpvf, iF),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_)
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

#ifdef OPENFOAM_NOT_EXTEND
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

    // Bug-fix, Mike Tree, see:
    // https://bitbucket.org/philip_cardiff/solids4foam-release/issues/27/elasticwallpressurefvpatchscalarfieldc
    // label patchID = this->patch().index(); // this is the fluid patch ID!

    // Find the solid patch ID corresponding to the current fluid patch
    label patchID = -1;
    forAll(fsi.fluidPatchIndices(), interfaceI)
    {
        if (fsi.fluidPatchIndices()[interfaceI] == patch().index())
        {
            // Take the corresponding solid patch ID
            patchID = fsi.solidPatchIndices()[interfaceI];
            break;
        }
    }

    if (patchID == -1)
    {
        FatalErrorIn
        (
            "void elasticWallPressureFvPatchScalarField::updateCoeffs()"
        )   << "Are you sure this patch is an FSI interface?"
            << abort(FatalError);
    }

    // Solid properties
    // PC: hmnn what if the solidModel does not use mu and lambda...
    // It seems that ap is the speed of sound so we just need the stiffness, we
    // can lookup impK
    // Also, what happends if mu/lambda are varying...

    // Get solid density
    const scalarField rhoSolid =
        fsi.solid().mechanical().rho()().boundaryField()[patchID];

    // Get solid stiffness (impK for generality)
    scalarField impK =
        fsi.solid().mechanical().impK()().boundaryField()[patchID];

    // p-wave propagation speed, ap
    const scalarField ap(sqrt(impK/rhoSolid));

    // Solid "virtual thickness"
    const scalarField hs(ap*mesh.time().deltaT().value());

    // Fluid properties
    const IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // Do not register
        )
    );

    const dimensionedScalar rhoFluid
    (
        transportProperties.lookup("rho")
    );

    if (debug)
    {
        Info<< "rhoSolid = " << max(rhoSolid)
            << ", hs = " << max(hs)
            << ", rhoFluid = " << rhoFluid.value()
            << endl;
    }

    // Update velocity and acceleration

    const fvPatch& p = patch();
    const vectorField n(p.nf());

#ifdef OPENFOAM_NOT_EXTEND
    const word fieldName = internalField().name();
#else
    const word fieldName = dimensionedInternalField().name();
#endif

    const volScalarField& pressure =
        mesh.lookupObject<volScalarField>(fieldName);

    // The previous acceleration is updated at the end of each
    // time step in the fluidSolidInterface
    const scalarField prevDdtUn(n & prevAcceleration_);

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

    robinFvPatchField<scalar>::updateCoeffs();
}


void elasticWallPressureFvPatchScalarField::patchFlux
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<scalar>& matrix
) const
{
    scalarField rAU(patch().size(), 0.0);
    if (db().foundObject<volScalarField>("rAU"))
    {
        Info<< "Found rAU" << endl;
        rAU = patch().lookupPatchField<volScalarField, scalar>("rAU");
    }
    else
    {
        rAU = patch().lookupPatchField<surfaceScalarField, scalar>("rAUf");
    }

#ifdef OPENFOAM_NOT_EXTEND
    flux.boundaryFieldRef()[patch().index()] = rAU*snGrad()*patch().magSf();
#else
    flux.boundaryField()[patch().index()] = rAU*snGrad()*patch().magSf();
#endif
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
