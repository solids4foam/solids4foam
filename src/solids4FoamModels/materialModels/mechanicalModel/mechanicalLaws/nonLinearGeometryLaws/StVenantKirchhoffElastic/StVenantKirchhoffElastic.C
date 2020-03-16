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
    lambda_("lambda", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    FPtr_(NULL),
    FfPtr_(NULL)
{
    // Read mechanical properties
    dimensionedScalar E("E", dimPressure, 0.0);
    dimensionedScalar nu("nu", dimless, 0.0);
    if
    (
        dict.found("E") && dict.found("nu")
     && !dict.found("mu") && !dict.found("K")
    )
    {
        E = dimensionedScalar(dict.lookup("E"));
        nu = dimensionedScalar(dict.lookup("nu"));

        mu_ = (E/(2.0*(1.0 + nu)));
    }
    else if
    (
        dict.found("mu") && dict.found("K")
     && !dict.found("E") && !dict.found("nu")
    )
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        const dimensionedScalar K = dimensionedScalar(dict.lookup("K"));

        E = 9*K*mu_/(3*K + mu_);
        nu = (3*K - 2*mu_)/(2*(3*K + mu_));
    }
    else
    {
        FatalErrorIn(type())
            << "Either E and nu or mu and K should be specified"
            << abort(FatalError);
    }

    if (planeStress())
    {
        lambda_ = nu*E/((1 + nu)*(1 - nu));
    }
    else
    {
        lambda_ = nu*E/((1 + nu)*(1 - 2.0*nu));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::StVenantKirchhoffElastic::~StVenantKirchhoffElastic()
{
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


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
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the relative deformation gradient: not needed
        const volTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF & F().oldTime();

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
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for incremental total Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Update the total deformation gradient
        F() = I + gradD.T();

        // Update the relative deformation gradient: not needed
        //relF() = F() & inv(F().oldTime());

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
    else
    {
        FatalErrorIn
        (
            "void Foam::StVenantKirchhoffElastic::"
            "correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
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
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(surfaceSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient: not needed
        const surfaceTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        Ff() = relF & Ff().oldTime();

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
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            FatalErrorIn(type() + "::correct(surfaceSymmTensorField& sigma)")
                << "Not implemented for incremental total Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Update the total deformation gradient
        Ff() = I + gradD.T();

        // Update the relative deformation gradient: not needed
        //relF() = F() & inv(F().oldTime());

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
    else
    {
        FatalErrorIn
        (
            "void Foam::StVenantKirchhoffElastic::"
            "correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
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
