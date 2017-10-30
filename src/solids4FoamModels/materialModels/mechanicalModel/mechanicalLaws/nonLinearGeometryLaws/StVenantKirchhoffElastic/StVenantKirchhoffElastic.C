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

#include "StVenantKirchhoffElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(StVenantKirchhoffElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, StVenantKirchhoffElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::StVenantKirchhoffElastic::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::StVenantKirchhoffElastic::makeF()")
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


Foam::volTensorField& Foam::StVenantKirchhoffElastic::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


void Foam::StVenantKirchhoffElastic::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::StVenantKirchhoffElastic::makeFf()")
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


Foam::surfaceTensorField& Foam::StVenantKirchhoffElastic::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::StVenantKirchhoffElastic::StVenantKirchhoffElastic
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
    lambda_
    (
        planeStress()
      ? (nu_*E_)/((1.0 + nu_)*(1.0 - nu_))
      : (nu_*E_)/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    mu_(E_/(2.0*(1.0 + nu_))),
    FPtr_(NULL),
    FfPtr_(NULL)
{
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::StVenantKirchhoffElastic::~StVenantKirchhoffElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::StVenantKirchhoffElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::StVenantKirchhoffElastic::impK() const
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


void Foam::StVenantKirchhoffElastic::correct(volSymmTensorField& sigma)
{
    // Update deformation gradient F
    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the total deformation gradient
        F() = (I + gradDD.T()) & F().oldTime();

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::StVenantKirchhoffElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu_*symm(gradDD) + lambda_*tr(gradDD)*I;

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

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::StVenantKirchhoffElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma = 2.0*mu_*symm(gradD) + lambda_*tr(gradD)*I;

            return;
        }
    }

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField c = symm(F().T() & F());

    // Calculate the Green strain tensor
    const volSymmTensorField E = 0.5*(c - I);

    // Calculate the 2nd Piola Kirchhoff stress
    const volSymmTensorField S = 2.0*mu_*E + lambda_*tr(E)*I;

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F());

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(F() & S & F().T());
    sigma = (1.0/J)*transform(F(), S);
}


void Foam::StVenantKirchhoffElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update deformation gradient F
    if (mesh().foundObject<volTensorField>("grad(DD)f"))
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the total deformation gradient
        Ff() = (I + gradDD.T()) & Ff().oldTime();

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::StVenantKirchhoffElastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu_*symm(gradDD) + lambda_*tr(gradDD)*I;

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

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::StVenantKirchhoffElastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma = 2.0*mu_*symm(gradD) + lambda_*tr(gradD)*I;

            return;
        }
    }

    // Calculate the right Cauchy–Green deformation tensor
    const surfaceSymmTensorField c = symm(Ff().T() & Ff());

    // Calculate the Green strain tensor
    const surfaceSymmTensorField E = 0.5*(c - I);

    // Calculate the 2nd Piola Kirchhoff stress
    const surfaceSymmTensorField S = 2.0*mu_*E + lambda_*tr(E)*I;

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J = det(Ff());

    // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
    // sigma = (1.0/J)*symm(Ff() & S & Ff().T());
    sigma = (1.0/J)*transform(Ff(), S);
}


// ************************************************************************* //
