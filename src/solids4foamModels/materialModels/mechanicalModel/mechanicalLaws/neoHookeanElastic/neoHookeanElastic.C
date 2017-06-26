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
#include "solidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::neoHookeanElastic::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElastic::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ =
        new volTensorField
        (
            IOobject
            (
                "F",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::neoHookeanElastic::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


void Foam::neoHookeanElastic::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElastic::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "Ff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::surfaceTensorField& Foam::neoHookeanElastic::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElastic::neoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    FPtr_(NULL),
    FfPtr_(NULL)
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
    // Check if the solidModel is enforcing linearity for convergence
    // If it is then we will calculate the stress using Hooke's law
    const Switch& enforceLinear =
        mesh().lookupObject<solidModel>("solidProperties").enforceLinear();

    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the total deformation gradient
        F() = (I + gradDD.T()) & F().oldTime();

        if (enforceLinear)
        {
            WarningIn
            (
                "void Foam::neoHookeanElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu_*symm(gradDD) + (K_ - (2.0/3.0)*mu_)*tr(gradDD)*I;

            return;
        }
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Update the total deformation gradient
        F() = I + gradD.T();

        if (enforceLinear)
        {
            WarningIn
            (
                "void Foam::neoHookeanElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma = 2.0*mu_*symm(gradD) + (K_ - (2.0/3.0)*mu_)*tr(gradD)*I;

            return;
        }
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F());

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField bEbar = pow(J, -2.0/3.0)*symm(F() & F().T());

    // Calculate the deviatoric stress
    const volSymmTensorField s = mu_*dev(bEbar);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


void Foam::neoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Check if the solidModel is enforcing linearity for convergence
    // If it is then we will calculate the stress using Hooke's law
    const Switch& enforceLinear =
        mesh().lookupObject<solidModel>("solidProperties").enforceLinear();

    if (mesh().foundObject<volTensorField>("grad(DD)f"))
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the total deformation gradient
        Ff() = (I + gradDD.T()) & Ff().oldTime();

            if (enforceLinear)
        {
            WarningIn
            (
                "void Foam::neoHookeanElastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu_*symm(gradDD) + (K_ - (2.0/3.0)*mu_)*tr(gradDD)*I;

            return;
        }
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Update the total deformation gradient
        Ff() = I + gradD.T();

        if (enforceLinear)
        {
            WarningIn
            (
                "void Foam::neoHookeanElastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma = 2.0*mu_*symm(gradD) + (K_ - (2.0/3.0)*mu_)*tr(gradD)*I;

            return;
        }
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J = det(Ff());

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField bEbar = pow(J, -2.0/3.0)*symm(Ff() & Ff().T());

    // Calculate deviatoric stress
    const surfaceSymmTensorField s = mu_*dev(bEbar);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


// ************************************************************************* //
