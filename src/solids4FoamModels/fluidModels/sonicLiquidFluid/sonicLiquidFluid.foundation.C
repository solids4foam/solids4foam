/*---------------------------------------------------------------------------*\ =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sonicLiquidFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "adjustPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "CorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sonicLiquidFluid, 0);
addToRunTimeSelectionTable(fluidModel, sonicLiquidFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void sonicLiquidFluid::compressibleContinuityErrs()
{
        scalar sumLocalContErr =
        (sum(mag(rho_ - rho0_ - psi_*(p() - p0_)))/sum(rho_)).value();

    scalar globalContErr =
                (
                        sum(rho_ - rho0_ - psi_*(p() - p0_))/sum(rho_)
                ).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = "
                << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_ << endl;
}


void sonicLiquidFluid::solveRhoEqn()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_)
      + fvc::div(phi())
     ==
        models().source(rho_)
    );

    constraints().constrain(rhoEqn);

    rhoEqn.solve();

    constraints().constrain(rho_);
}


void sonicLiquidFluid::CorrectFlux()
{
    Info<< "Correcting flux for moving mesh" << endl;

    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    autoPtr<volScalarField> divrhoU;

    divrhoU.set
    (
        new volScalarField
        (
                "divrhoU",
                fvc::div(phi())
        )
    );

    // Store momentum to set rhoUf for introduced faces.
    autoPtr<volVectorField> rhoU;

    // Calculate absolute flux
    // from the mapped surface velocity
    phi() = mesh().Sf() & rhoUf_();

    // Define volScalarField to hold phi
    // to pass to compressible CorrectPhi function
    const volScalarField psi
    (
         IOobject
         (
                "psi",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
         ),
         mesh(),
         psi_
    );

    CorrectPhi
    (
        phi(),
        p(),
        rho_,
        psi,
        dimensionedScalar("rAUf", dimTime, 1),
        divrhoU(),
        pimple()
        //,
        //true
    );
}

void sonicLiquidFluid::compressibleCourantNo()
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    surfaceScalarField SfUfbyDelta
        (
            mesh().surfaceInterpolation::deltaCoeffs()*mag(phi())
            /fvc::interpolate(rho_)
        );

    CoNum = max
            (
                SfUfbyDelta/mesh().magSf()
            ).value()*runTime().deltaT().value();

    meanCoNum = (
                    sum(SfUfbyDelta)/sum(mesh().magSf())
                ).value()*runTime().deltaT().value();

    Info<< "Region: " << mesh().name()
        << " Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sonicLiquidFluid::sonicLiquidFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mu_(transportProperties_.lookup("mu")),
    thermodynamicProperties_
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rho0_(thermodynamicProperties_.lookup("rho0")),
    p0_(thermodynamicProperties_.lookup("p0")),
    // pis (compressibility) can be read directly or calculated from the bulk
    // modulus
    psi_
    (
        bool(thermodynamicProperties_.found("psi"))
      ? dimensionedScalar(thermodynamicProperties_.lookup("psi"))
      : rho0_/dimensionedScalar(thermodynamicProperties_.lookup("K"))
    ),
    rhoO_
    (
        IOobject
        (
            "rhoO",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        rho0_ - psi_*p0_
    ),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoO_ + psi_*p()
    ),
    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        runTime.deltaT()/rho0_,
        zeroGradientFvPatchScalarField::typeName
    ),
    cumulativeContErr_(0)
{
    UisRequired();
    pisRequired();

    // Reset phi dimensions: compressible
    Info<< "Resetting the dimensions of phi" << endl;
    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);

    phi() = linearInterpolate(rho_*U()) & mesh().Sf();

    mesh().setFluxRequired(p().name());

    if (mesh().dynamic())
    {
        Info<< "Constructing face momentum rhoUf" << endl;

        rhoUf_.set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "rhoUf",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rho_*U())
            )
        );
    }

    // Check that either psi or K are specified, not both
    if
    (
        thermodynamicProperties_.found("psi")
     && thermodynamicProperties_.found("K")
    )
    {
        FatalErrorIn(type() + "::" + type() + "()")
            << "Either psi OR K should be specified in "
            << thermodynamicProperties_.name() << ", not both!"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> sonicLiquidFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF.ref() = mu_.value()*U().boundaryField()[patchID].snGrad();

    return tvF;
}


tmp<scalarField> sonicLiquidFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    // Pressure here is already in Pa
    tpF.ref() = p().boundaryField()[patchID];

    return tpF;
}


bool sonicLiquidFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    dynamicFvMesh& mesh = this->mesh();

    bool meshChanged = false;

    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());
    }

    if (meshChanged)
    {
        const Time& runTime = fluidModel::runTime();
#       include "volContinuity.H"
    }

    bool correctPhi
    (
        pimple().dict().lookupOrDefault("correctPhi", false)
    );

    if (correctPhi && meshChanged)
    {
        CorrectFlux();
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), rho_, U());

    // Calculate CourantNo
    compressibleCourantNo();

    // Pressure-velocity corrector
    while (pimple().loop())
    {
        if (!pimple().simpleRho())
        {
            solveRhoEqn();
        }

        // --- Momentum equation
        // Time-derivative matrix
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho_, U())
          + fvm::div(phi(), U())
          - fvm::laplacian(mu_, U())
         ==
            models().source(rho_, U())
        );

        UEqn.relax();

        constraints().constrain(UEqn);

        if (pimple().momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p()));
        }
        else
        {
            // Explicit update
            U() = (UEqn.H() - fvc::grad(p()))/UEqn.A();
            U().correctBoundaryConditions();
        }

        constraints().constrain(U());

        // --- Pressure corrector loop
        while (pimple().correct())
        {
            rAU_ = 1.0/UEqn.A();

            const surfaceScalarField rhorAUf
            (
                "rhorAUf",
                fvc::interpolate(rho_*rAU_)
            );

            U() = rAU_*UEqn.H();

            // Compute advective transport coeff. to pressure
            surfaceScalarField phid
            (
                "phid",
                psi_
               *(
                    fvc::flux(U())
                  + rhorAUf*fvc::ddtCorr(rho_, U(), phi())/fvc::interpolate(rho_)
                )
            );

            // Make flux relative to mesh motion
            fvc::makeRelative(phid, psi_, U());

            // Recompute flux for pressure equation
            phi() = fvc::interpolate(rhoO_)*phid/psi_;

            // Pressure equation
            // essential to account for ddt(rho0) term for mesh motion)
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi_, p())
              + fvc::div(phi())
              + fvm::div(phid, p())
              - fvm::laplacian(rhorAUf, p())
              + fvc::ddt(rhoO_)
            );

            pEqn.solve();

            phi() += pEqn.flux();

            solveRhoEqn();

            // State equation errors
            compressibleContinuityErrs();

            // Correct velocity
            U() -= rAU_*fvc::grad(p());
            U().correctBoundaryConditions();
        }

        gradU() = fvc::grad(U());
    }

    // Update rho based on equation of state
    rho_ = rhoO_ + psi_*p();

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi(), rho_, U());

    // Print variables for inspection
    // Density variation
    scalar deltaRho = max(rho_).value() - min(rho_).value();
    scalar refDeltaRho = 0.01*rho0_.value();

    // Info variable values
    Info<< nl << "Density: min " << min(rho_).value()
        << " max " << max(rho_).value() << endl;

    Info<< "Density variation: " << deltaRho
        << " ref: " << refDeltaRho << endl;

    Info<< "Pressure: min " << min(p()).value()
        << " max " << max(p()).value() << endl;

    Info<< "Velocity: min " << min(mag(U())).value()
        << " max " << max(mag(U())).value() << nl << nl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
