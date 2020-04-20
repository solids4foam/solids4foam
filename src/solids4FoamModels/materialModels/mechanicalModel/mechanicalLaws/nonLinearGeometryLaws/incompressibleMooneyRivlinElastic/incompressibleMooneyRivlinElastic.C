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

#include "incompressibleMooneyRivlinElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleMooneyRivlinElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, incompressibleMooneyRivlinElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::incompressibleMooneyRivlinElastic::calculateHydrostaticStress
(
    volScalarField& sigmaHyd,
    const volScalarField& Jacobian
)
{
    if (solvePressureEqn_)
    {
        // Store previous iteration to allow relaxation, if needed
        sigmaHyd.storePrevIter();

        // Lookup the momentum equation inverse diagonal field
        const volScalarField* ADPtr = NULL;
        if (mesh().foundObject<volScalarField>("DEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DEqnA");
        }
        else if (mesh().foundObject<volScalarField>("DDEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DDEqnA");
        }
        else
        {
            FatalErrorIn
            (
                "void Foam::incompressibleMooneyRivlinElastic::"
                "calculateHydrostaticStress\n"
                "(\n"
                "    volScalarField& sigmaHyd,\n"
                "    const volScalarField& trEpsilon\n"
                ")"
            )   << "Cannot find the DEqnA or DDEqnA field: this should be "
                << "stored in the solidModel" << abort(FatalError);
        }
        const volScalarField& AD = *ADPtr;

        // Pressure diffusivity field multiple by (4.0/3.0)*mu + K, which is
        // equivalent to (2*mu + lambda)
        // Note: we can scale this coefficient by pressureSmoothingCoeff to
        // provide greater smoothing
        const surfaceScalarField rDAf
        (
            "rDAf",
            pressureSmoothingCoeff_
           *fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        const dimensionedScalar one("one", dimless, 1.0);

        // Solve pressure laplacian
        // Note: the the laplacian and div terms combine to produce a
        // third-order smoothing/dissipation term
        // If we only used the laplacian term then the smoothing/dissipation
        // would be second-order.
        // It would be interesting to see how this compares to the JST 4th
        // order dissipation term
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(one, sigmaHyd)
          - fvm::laplacian(rDAf, sigmaHyd, "laplacian(DA,sigmaHyd)")
         ==
            0.5*K_*(pow(Jacobian, 2.0) - 1.0)
          - fvc::div(rDAf*fvc::interpolate(gradSigmaHyd_) & mesh().Sf())
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the pressure field
        sigmaHyd.relax();
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = 0.5*K_*(pow(Jacobian, 2.0) - 1.0);
    }

    // Update the gradient
    gradSigmaHyd_ = fvc::grad(sigmaHyd_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::incompressibleMooneyRivlinElastic::incompressibleMooneyRivlinElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    c10_(dict.lookup("c10")),
    c01_(dict.lookup("c01")),
    c11_(dict.lookup("c11")),
    mu_(2*(c10_ + c01_)),
    K_(dict.lookup("K")),
    solvePressureEqn_
    (
        dict.lookupOrDefault<Switch>("solvePressureEqn", false)
    ),
    pressureSmoothingCoeff_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingCoeff", 1.0)
    ),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    gradSigmaHyd_
    (
        IOobject
        (
            "grad(" + sigmaHyd_.name() + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimPressure/dimLength, vector::zero),
        calculatedFvPatchScalarField::typeName
    )
{
    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingCoeff: " << pressureSmoothingCoeff_
            << endl;
    }

    if (planeStress())
    {
        notImplemented
        (
            type() + " mechanical law is not implemented for planeStress"
        );
    }

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleMooneyRivlinElastic::~incompressibleMooneyRivlinElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::incompressibleMooneyRivlinElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::incompressibleMooneyRivlinElastic::impK() const
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
            (4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::incompressibleMooneyRivlinElastic::correct
(
    volSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F());

    // Calculate the hydrostatic stress
    calculateHydrostaticStress(sigmaHyd_, J);

    // Calculate the left Cauchy-Green deformation tensor
    const volSymmTensorField B = symm(F() & F().T());

    // Compute invariants fields
    const volScalarField I1 = tr(B);
    const volScalarField I2 = 0.5*(pow(tr(B), 2.0) - tr(B & B));

    // Calculate the deviatoric part of the Cauchy stress
    const volSymmTensorField s =
        2*
        (   c10_
          + I1*c01_
          + c11_*(I1*(I1 - 3) + I2 - 3)
        )*B
      - 2*
        (
            c01_
          + c11_*(I1 - 3)
        )*symm(B & B);

    // Calculate the Cauchy stress
    // for the Mooney-Rivlin model
    // The last RHS term is the initial hydrostatic pressure field
    // This term is important to assure the underformed configuration
    // to be stress-free
    sigma = sigmaHyd_*I + s - 2.0*(c10_ + 2.0*c01_)*I;
}


void Foam::incompressibleMooneyRivlinElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J = det(Ff());

    // To-do: interpolate from vol field
    if (solvePressureEqn_)
    {
        notImplemented
        (
            "solvePressureEqn cannot currently be used with 'uns' solid models"
        );
    }

    // Calculate pressure field with bulk modulus to approximate
    // incompressibility
    const surfaceScalarField sigmaHydf = 0.5*K_*(pow(J, 2.0) - 1.0);

    // Calculate the left Cauchy-Green deformation tensor
    const surfaceSymmTensorField Bf = symm(Ff() & Ff().T());

    // Compute invariants fields
    const surfaceScalarField I1 = tr(Bf);
    const surfaceScalarField I2 = 0.5*(pow(tr(Bf), 2.0) - tr(Bf & Bf));

    // Calculate the deviatoric part of the Cauchy stress
    const surfaceSymmTensorField s =
        2*
        (   c10_
          + I1*c01_
          + c11_*(I1*(I1 - 3) + I2 - 3)
        )*Bf
      - 2*
        (
            c01_
          + c11_*(I1 - 3)
        )*symm(Bf & Bf);

    // Calculate the Cauchy stress
    // for the Mooney-Rivlin model
    // The last RHS term is the initial hydrostatic pressure field
    // This term is important to assure the underformed configuration
    // to be stress-free
    sigma = sigmaHydf*I + s - 2.0*(c10_ + 2.0*c01_)*I;
}


// ************************************************************************* //
