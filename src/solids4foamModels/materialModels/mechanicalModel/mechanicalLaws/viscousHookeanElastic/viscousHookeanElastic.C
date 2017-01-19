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

#include "viscousHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscousHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, viscousHookeanElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::viscousHookeanElastic::viscousHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    rho_(dict.lookup("rho")),
    EInf_(dict.lookup("EInfinity")),
    E_(dict.lookup("E")),
    tau_(dict.lookup("relaxationTimes")),
    gammaInf_(0.0),
    gamma_(E_.size(), 0.0),
    nu_(dict.lookup("nu")),
    lambda_("lambda", dimPressure, 0.0),
    mu_(EInf_/(2.0*(1.0 + nu_))),
    k_("k", dimPressure, 0.0),
    h_(),
    hf_(),
    s_
    (
        IOobject
        (
            "s",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    ),
    sf_
    (
        IOobject
        (
            "sf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    )
{
    // Check for physical Poisson's ratio
    if (nu_.value() < -1.0 || nu_.value() > 0.5)
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
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
            lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - nu_));
        }
        else
        {
            lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - 2.0*nu_));
        }

        k_ = lambda_ + (2.0/3.0)*mu_;
    }

    // Check E_ and tau_ are the same length
    if (E_.size() != tau_.size())
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The E and tau lists should have the same length!"
            << abort(FatalError);
    }

    // Calculate relative modulii

    const scalar E0 = EInf_.value() + sum(E_);

    gammaInf_ = EInf_.value()/E0;

    forAll(gamma_, i)
    {
        gamma_[i] = E_[i]/E0;
    }

    // Check all the relaxation times are positive
    if (min(tau_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "All relaxation times should be positive!"
            << abort(FatalError);
    }

    // Check all the E values are positive
    if (min(E_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::viscousHookeanElastic::viscousHookeanElastic\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "All values of stiffness E should be positive!"
            << abort(FatalError);
    }

    // Print out the relative module

    Info<< "Relative modulii" << nl
        << "    gammaInfinity: " << gammaInf_ << nl;

    forAll(gamma_, i)
    {
        Info<< "    gamma[" << i << "] : " << gamma_[i] << nl;
    }

    Info<< endl;

    // Create the internal stress variables for each Maxwell model

    h_.setSize(gamma_.size());
    hf_.setSize(gamma_.size());

    forAll(h_, MaxwellModelI)
    {
        h_.set
        (
            MaxwellModelI,
            new volSymmTensorField
            (
                IOobject
                (
                    "h" + name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );

        hf_.set
        (
            MaxwellModelI,
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "hf" + name(MaxwellModelI),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );

        // We need to store the old time field
        h_[MaxwellModelI].storeOldTime();
        hf_[MaxwellModelI].storeOldTime();
    }

    // Store the old time s field
    s_.storeOldTime();
    sf_.storeOldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::viscousHookeanElastic::~viscousHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscousHookeanElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::viscousHookeanElastic::impK() const
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


Foam::tmp<Foam::volScalarField> Foam::viscousHookeanElastic::K() const
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


void Foam::viscousHookeanElastic::correct(volSymmTensorField& sigma)
{
    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate deviatoric component of the strain increment
        const volSymmTensorField De = dev(symm(gradDD));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        s_ = s_.oldTime() + 2.0*mu_*De;

        // Set the volumetric component of the total stress
        sigma = tr(sigma.oldTime())*symmTensor(I) + k_*tr(gradDD)*I;
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate deviatoric component of total strain
        const volSymmTensorField e = dev(symm(gradD));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        s_ = 2.0*mu_*e;

        // Set the volumetric component of the total stress
        sigma = k_*tr(gradD)*symmTensor(I);
    }

    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model

    const scalar deltaT = mesh().time().deltaTValue();

    forAll(h_, MaxwellModelI)
    {
        h_[MaxwellModelI] =
            Foam::exp(-deltaT/tau_[MaxwellModelI])*h_[MaxwellModelI].oldTime()
          + Foam::exp(-deltaT/(2.0*tau_[MaxwellModelI]))*(s_ - s_.oldTime());
    }

    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic
    sigma += gammaInf_*s_;

    forAll(h_, MaxwellModelI)
    {
        sigma += gamma_[MaxwellModelI]*h_[MaxwellModelI];
    }
}


void Foam::viscousHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Calculate deviatoric component of the strain increment
        const surfaceSymmTensorField De = dev(symm(gradDD));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        sf_ = sf_.oldTime() + 2.0*mu_*De;

        // Set the volumetric component of the total stress
        sigma = tr(sigma.oldTime())*symmTensor(I) + k_*tr(gradDD)*I;
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Calculate deviatoric component of total strain
        const surfaceSymmTensorField e = dev(symm(gradD));

        // Calculate deviatoric component of the initial stress, based on
        // Hooke's law
        sf_ = 2.0*mu_*e;

        // Set the volumetric component of the total stress
        sigma = k_*tr(gradD)*symmTensor(I);
    }

    // Update internal stress variables, representing stress relaxations for
    // each Maxwell model

    const scalar deltaT = mesh().time().deltaTValue();

    forAll(hf_, MaxwellModelI)
    {
        hf_[MaxwellModelI] =
            Foam::exp(-deltaT/tau_[MaxwellModelI])*hf_[MaxwellModelI].oldTime()
          + Foam::exp(-deltaT/(2.0*tau_[MaxwellModelI]))*(sf_ - sf_.oldTime());
    }

    // Calculate the current total stress, where the volumetric term is
    // elastic and the deviatoric term is viscoelastic
    sigma += gammaInf_*sf_;

    forAll(hf_, MaxwellModelI)
    {
        sigma += gamma_[MaxwellModelI]*hf_[MaxwellModelI];
    }
}


// ************************************************************************* //
