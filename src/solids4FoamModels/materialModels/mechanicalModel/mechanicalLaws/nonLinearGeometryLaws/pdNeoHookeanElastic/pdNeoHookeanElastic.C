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

#include "pdNeoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pdNeoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, pdNeoHookeanElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::pdNeoHookeanElastic::pdNeoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    nu_("nu", dimless, 0.5),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0)
{
    // Read mechanical properties
    if
    (
        dict.found("E") && dict.found("nu")
     && !dict.found("mu")
    )
    {
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));
        nu_ = dimensionedScalar(dict.lookup("nu"));

        mu_ = (E/(2.0*(1.0 + nu_)));

        // Set the bulk modulus
        if (nu_.value() < 0.5-SMALL)
        {
            K_ = (nu_*E/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_.value() = GREAT;
        }
    }
    // else if
    // (
    //     dict.found("mu")
    //  && !dict.found("E") && !dict.found("nu")
    // )
    // {
    //     mu_ = dimensionedScalar(dict.lookup("mu"));
    //     K_.value() = GREAT;
    // }
    else
    {
        FatalErrorIn(type())
            << "One should specify E and nu"
            << abort(FatalError);
    }

    // Read F and Ff if present (need for restart)
    F();
    Ff();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdNeoHookeanElastic::~pdNeoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::pdNeoHookeanElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::pdNeoHookeanElastic::impK() const
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
            mu_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::pdNeoHookeanElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulusLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            K_
        )
    );
}


void Foam::pdNeoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    dimensionedScalar K_("K", mu_.dimensions(), 0);
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    const volTensorField& F = this->F();
    const volTensorField FT(F.T());

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F);

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField b = symm(F & FT);

    // Calculate the deviatoric stress
    const volSymmTensorField s = mu_*(b-I)/J;

    // Pressure term coeffs
    dimensionedScalar lambdaByK = 3*nu_/(1+nu_);

    // Lookup pressure field
    const volScalarField& p =
        mesh().lookupObject<volScalarField>("p");

    // Calculate the Cauchy stress
    sigma = -lambdaByK*p*I + s;
}


void Foam::pdNeoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    dimensionedScalar K_("K", mu_.dimensions(), 0);
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    const surfaceTensorField& Ff = this->Ff();
    const surfaceTensorField FfT(Ff.T());

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField Jf = det(Ff);

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField b = symm(Ff & FfT);

    // Calculate deviatoric stress
    const surfaceSymmTensorField s = mu_*(b-I)/Jf;

    // Lookup pressure field
    const surfaceScalarField& pf =
        mesh().lookupObject<surfaceScalarField>("pf");

    // Pressure term coeffs
    dimensionedScalar lambdaByK = 3*nu_/(1+nu_);

    // Calculate the Cauchy stress
    sigma = -lambdaByK*pf*I + s;
}


// ************************************************************************* //
