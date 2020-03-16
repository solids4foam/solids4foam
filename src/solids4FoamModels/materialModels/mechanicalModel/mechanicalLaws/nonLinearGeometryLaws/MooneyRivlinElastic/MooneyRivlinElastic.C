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

#include "MooneyRivlinElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MooneyRivlinElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, MooneyRivlinElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::MooneyRivlinElastic::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::MooneyRivlinElastic::makeF()")
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


Foam::volTensorField& Foam::MooneyRivlinElastic::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


void Foam::MooneyRivlinElastic::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::MooneyRivlinElastic::makeFf()")
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


Foam::surfaceTensorField& Foam::MooneyRivlinElastic::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::MooneyRivlinElastic::MooneyRivlinElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    mu1_(dict.lookup("mu1")),
    mu2_(dict.lookup("mu2")),
    K_(dict.lookup("K")),
    FPtr_(NULL),
    FfPtr_(NULL)
{
    if (planeStress())
    {
        notImplemented
        (
            type() + " mechanical law is not implemented for planeStress"
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MooneyRivlinElastic::~MooneyRivlinElastic()
{
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::MooneyRivlinElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::MooneyRivlinElastic::impK() const
{
    const dimensionedScalar mu = mu1_ + mu2_;

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
            (4.0/3.0)*mu + K_ // == 2*mu + lambda
        )
    );
}


void Foam::MooneyRivlinElastic::correct(volSymmTensorField& sigma)
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

        // Calculate the relative deformation gradient
        const volTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF & F().oldTime();

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::MooneyRivlinElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            const dimensionedScalar mu = mu1_ + mu2_;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu*symm(gradDD) + (K_ - (2.0/3.0)*mu)*tr(gradDD)*I;

            return;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            F() = F().oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            //relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::MooneyRivlinElastic::"
                    "correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                const dimensionedScalar mu = mu1_ + mu2_;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K_*tr(gradDD)*I;

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

            // Update the relative deformation gradient: not needed
            //relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::MooneyRivlinElastic::"
                    "correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                const dimensionedScalar mu = mu1_ + mu2_;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K_*tr(gradD)*I;

                return;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::MooneyRivlinElastic::"
            "correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F());

    // Calculate the left Cauchy Green strain
    const volSymmTensorField b = symm(F() & F().T());

    // Calculate hydrostatic pressure term
    const volScalarField p = -K_*(J - 1.0);

    // Calculate deviatoric stress term
    const volSymmTensorField s =
        (mu1_/pow(J, 5.0/3.0))*dev(b)
      + (mu2_/pow(J, 7.0/3.0))
       *symm
        (
           tr(b)*b - pow(tr(b), 2)*I/3.0 - (b & b) + (b && b)*I/3
        );

    // Calculate the Cauchy stress
    sigma = s - p*I;
}


void Foam::MooneyRivlinElastic::correct(surfaceSymmTensorField& sigma)
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
                "void Foam::MooneyRivlinElastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            const dimensionedScalar mu = mu1_ + mu2_;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu*dev(symm(gradDD)) + K_*tr(gradDD)*I;

            return;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const surfaceTensorField& gradDD =
                mesh().lookupObject<surfaceTensorField>("grad(DD)f");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            Ff() = Ff().oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            //relFf() = Ff() & inv(Ff().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::MooneyRivlinElastic::"
                    "correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                const dimensionedScalar mu = mu1_ + mu2_;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K_*tr(gradDD)*I;

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

            // Update the relative deformation gradient: not needed
            //relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::MooneyRivlinElastic::"
                    "correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                const dimensionedScalar mu = mu1_ + mu2_;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K_*tr(gradD)*I;

                return;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::MooneyRivlinElastic::"
            "correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J = det(Ff());

    // Calculate the left Cauchy Green strain
    const surfaceSymmTensorField b = symm(Ff() & Ff().T());

    // Calculate hydrostatic pressure term
    const surfaceScalarField p = -K_*(J - 1.0);

    // Calculate deviatoric stress term
    const surfaceSymmTensorField s =
        (mu1_/pow(J, 5.0/3.0))*dev(b)
      + (mu2_/pow(J, 7.0/3.0))
       *symm
        (
           tr(b)*b - pow(tr(b), 2)*I/3.0 - (b & b) + (b && b)*I/3
        );

    // Calculate the Cauchy stress
    sigma = s - p*I;
}


// ************************************************************************* //
