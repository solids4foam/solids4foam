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

InClass
    Foam::mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "lookupSolidModel.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalLaw, 0);
    defineRunTimeSelectionTable(mechanicalLaw, linGeomMechLaw);
    defineRunTimeSelectionTable(mechanicalLaw, nonLinGeomMechLaw);
}

// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::mechanicalLaw::makeF()
{
    if (FPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_.set
    (
        new volTensorField
        (
            IOobject
            (
                "F+" + type(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeFf()
{
    if (FfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_.set
    (
        new surfaceTensorField
        (
            IOobject
            (
                "Ff_" + type(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


// * * * * * * * * * * * * * * Protected Members * * * * * * * * * * * * * * //

bool Foam::mechanicalLaw::planeStress() const
{
    if (mesh_.foundObject<IOdictionary>("mechanicalProperties"))
    {
        return
            Switch
            (
                mesh_.lookupObject<IOdictionary>
                (
                    "mechanicalProperties"
                ).lookup("planeStress")
            );
    }
    else
    {
        // It is not straight-forward to lookup the mechanicalProperties from
        // here as we only have access to a subMesh fvMesh objectRegistry
        // We will read it here again; this switch only gets called at the start
        // of a simulation so it is not a problem
        IOdictionary mechProp
        (
            IOobject
            (
                "mechanicalProperties",
                "constant",
                mesh_.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        return Switch(mechProp.lookup("planeStress"));
    }
}


Foam::volTensorField& Foam::mechanicalLaw::F()
{
    if (FPtr_.empty())
    {
        makeF();
    }

    return FPtr_();
}


Foam::surfaceTensorField& Foam::mechanicalLaw::Ff()
{
    if (FfPtr_.empty())
    {
        makeFf();
    }

    return FfPtr_();
}


bool Foam::mechanicalLaw::updateF
(
    volSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
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
                "void Foam::MooneyRivlinThreeParametersElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu*symm(gradDD) + (K - (2.0/3.0)*mu)*tr(gradDD)*I;

            return true;
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
                    "void " + type() + "::correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                return true;
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
                    "void " + type() + "::correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;

                return true;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "::correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // linearised elasticity was not enforced
    return false;
}


bool Foam::mechanicalLaw::updateFf
(
    surfaceSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
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
                "void " + type() + "::correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

            return true;
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
                    "void " + type()
                  + "::correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                return true;
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
                    "void " + type()
                  + "::correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;

                return true;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "::correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // linearised elasticity was not enforced
    return false;
}


const Foam::Switch& Foam::mechanicalLaw::enforceLinear() const
{
    // Lookup the solideModel
    const solidModel& solMod = lookupSolidModel(mesh(), baseMeshRegionName_);

    return solMod.enforceLinear();
}


bool Foam::mechanicalLaw::incremental() const
{
    // Lookup the solideModel
    const solidModel& solMod = lookupSolidModel(mesh(), baseMeshRegionName_);

    return solMod.incremental();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalLaw::mechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    baseMeshRegionName_(),
    nonLinGeom_(nonLinGeom),
    FPtr_(),
    FfPtr_()
{
    // Set the base mesh region name
    // For an FSI case, the region will be called solid, else it will be called
    // region0.
    if (mesh.time().foundObject<fvMesh>("solid"))
    {
        baseMeshRegionName_ = "solid";
    }
    else if (mesh.time().foundObject<fvMesh>("region0"))
    {
        baseMeshRegionName_ = "region0";
    }
    else
    {
        FatalErrorIn
        (
            "Foam::mechanicalLaw::mechanicalLaw\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        ) << "solid region name not found" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalLaw::impKf() const
{
    return fvc::interpolate(impK());
}


void Foam::mechanicalLaw::correct(surfaceSymmTensorField&)
{
    notImplemented
    (
        type() + "::correct(surfaceSymmTensorField&)\n"
        "The correct(surfaceSymmTensorField&) function is not implemented\n"
        " for the " + type() + " mechanical law"
    );
}


Foam::scalar Foam::mechanicalLaw::residual()
{
    // Default to zero; this can be overwritten by any derived mechanical law
    return 0.0;
}


Foam::scalar Foam::mechanicalLaw::newDeltaT()
{
    // Default to a large number
    return mesh_.time().endTime().value();
}


// ************************************************************************* //
