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

#include "neoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElastic::neoHookeanElastic
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
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        mesh.lookupObject<IOdictionary>
        (
            "mechanicalProperties"
        ).lookup("planeStress")
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElastic::~neoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::impK() const
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
            (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
        )
    );
}


void Foam::neoHookeanElastic::correct(volSymmTensorField& sigma)
{
    const fvMesh& mesh = this->mesh();

    // Lookup the total deformation gradient from the solver
    const volTensorField& F = mesh.lookupObject<volTensorField>("F");

    // Lookup the Jacobian of the deformation gradient from the solver
    const volScalarField& J = mesh.lookupObject<volScalarField>("J");

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    volSymmTensorField bEbar = pow(J, -2.0/3.0)*symm(F & F.T());

    // Calculate deviatoric stress
    volSymmTensorField s = mu_*dev(bEbar);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*0.5*K_*(pow(J, 2) - 1)*I + s;
}


void Foam::neoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    const fvMesh& mesh = this->mesh();

    // Lookup the total deformation gradient from the solver
    const surfaceTensorField& F = mesh.lookupObject<surfaceTensorField>("Ff");

    // Lookup the Jacobian of the deformation gradient from the solver
    const surfaceScalarField& J = mesh.lookupObject<surfaceScalarField>("Jf");

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    surfaceSymmTensorField bEbar = pow(J, -2.0/3.0)*symm(F & F.T());

    // Calculate deviatoric stress
    surfaceSymmTensorField s = mu_*dev(bEbar);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*0.5*K_*(pow(J, 2) - 1)*I + s;
}


void Foam::neoHookeanElastic::setMaterialIndex(label curMatIndex)
{
    // Set current material index
    curMaterialIndex() = curMatIndex;

    // Rename fields to avoid conflicts with other mechanical laws
    curMaterial().rename(curMaterial().name() + '_' + Foam::name(curMatIndex));
}


// ************************************************************************* //
