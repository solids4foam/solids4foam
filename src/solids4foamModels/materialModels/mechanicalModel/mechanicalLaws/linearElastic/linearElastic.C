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

#include "linearElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElastic::linearElastic
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
    lambda_("lambda", dimPressure, 0.0),
    mu_(E_/(2.0*(1.0 + nu_))),
    k_("k", dimPressure, 0.0)
{
    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorIn
        (
            "Foam::linearElastic::linearElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unphysical Poisson's ratio: nu should be >= -1.0 and <= 0.5"
            << abort(FatalError);
    }

    // Check for incompressibility
    if (nu_.value() == 0.5)
    {
        Info<< "Material is incompressible: make sure to use a hybrid"
            << " approach solid model" << endl;

        // Set lambda and k to GREAT
        lambda_.value() = GREAT;
        k_.value() = GREAT;
    }
    else
    {
        if (planeStress())
        {
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
        }
        else
        {
            lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        }

        k_ = lambda_ + (2.0/3.0)*mu_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElastic::~linearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::linearElastic::impK() const
{
    if (nu_.value() == 0.5 || mesh().foundObject<volScalarField>("p"))
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
                2.0*mu_
            )
        );
    }
    else
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

}


Foam::tmp<Foam::volScalarField> Foam::linearElastic::K() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            k_
        )
    );
}


const Foam::dimensionedScalar& Foam::linearElastic::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar& Foam::linearElastic::lambda() const
{
    return lambda_;
}


void Foam::linearElastic::correct(volSymmTensorField& sigma)
{
    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate stress based on incremental form of Hooke's law
        sigma = sigma.oldTime() + mu_*twoSymm(gradDD) + lambda_*tr(gradDD)*I;
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Check if a hybrid approach is being used
        // Currently, this will only work when there is no material interface
        if (mesh().foundObject<volScalarField>("p"))
        {
            const volScalarField& p = mesh().lookupObject<volScalarField>("p");

            // Calculate stress using Hooke's law in uncoupled deviatoric
            // volumetric form
            sigma = mu_*dev(twoSymm(gradD)) - p*I;
        }
        else
        {
            // Calculate stress based on Hooke's law
            sigma = mu_*twoSymm(gradD) + lambda_*tr(gradD)*I;
        }
    }
}


void Foam::linearElastic::correct(surfaceSymmTensorField& sigma)
{
    if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Calculate stress based on incremental form of Hooke's law
        sigma = sigma.oldTime() + mu_*twoSymm(gradDD) + lambda_*tr(gradDD)*I;
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Calculate stress based on Hooke's law
        sigma = mu_*twoSymm(gradD) + lambda_*tr(gradD)*I;
    }
}


// ************************************************************************* //
