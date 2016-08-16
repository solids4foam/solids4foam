/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "thermoLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IOdictionary.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, thermoLinearElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::thermoLinearElastic::thermoLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    lambda_
    (
        planeStress()
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_(lambda_ + (2.0/3.0)*mu_),
    alpha_(dict.lookup("alpha"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoLinearElastic::~thermoLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::thermoLinearElastic::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::thermoLinearElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            2.0*mu_ + lambda_
        )
    );
}


const Foam::dimensionedScalar& Foam::thermoLinearElastic::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar& Foam::thermoLinearElastic::lambda() const
{
    return lambda_;
}


void Foam::thermoLinearElastic::correct(volSymmTensorField& sigma)
{
    // Lookup the strain tensor from the solver
    const volSymmTensorField epsilon =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<symmTensor>("epsilon", mesh());

    // Lookup the temperature field from the solver
    const volScalarField T =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T", mesh());

    // Lookup the stress-free reference temperature field from the solver
    const volScalarField T0 =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T0", mesh());

    // Calculate stress based on Hooke's law
    sigma = 2.0*mu_*epsilon + lambda_*tr(epsilon)*I - 3.0*K_*alpha_*(T - T0)*I;
}


void Foam::thermoLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Lookup the strain tensor from the solver
    const surfaceSymmTensorField epsilon =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshSurfaceField<symmTensor>("epsilonf", mesh());

    // Lookup the temperature field from the solver
    const volScalarField T =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T", mesh());

    const surfaceScalarField& Tf = fvc::interpolate(T);

    // Lookup the stress-free reference temperature field from the solver
    const volScalarField T0 =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T0", mesh());

    const surfaceScalarField& T0f = fvc::interpolate(T0);

    // Calculate stress based on Hooke's law
    sigma =
        2.0*mu_*epsilon + lambda_*tr(epsilon)*I - 3.0*K_*alpha_*(Tf - T0f)*I;
}


// ************************************************************************* //
