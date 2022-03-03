/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAMESIORFOUNDATION

#include "buoyantPimpleFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "CorrectPhi.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(buoyantPimpleFluid, 0);
addToRunTimeSelectionTable(physicsModel, buoyantPimpleFluid, fluid);
addToRunTimeSelectionTable(fluidModel, buoyantPimpleFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void buoyantPimpleFluid::compressibleContinuityErrs()
{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho_);

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho_ - thermo_.rho()))/totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rho_ - thermo_.rho())/totalMass).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void buoyantPimpleFluid::solveRhoEqn()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_)
      + fvc::div(phi())
#ifdef OPENFOAMFOUNDATION
     ==
        models().source(rho_)
#elif OPENFOAMESI
     ==
        options()(rho_)
#endif
    );

#ifdef OPENFOAMFOUNDATION
        constraints().constrain(rhoEqn);
#elif OPENFOAMESI
        options().constrain(rhoEqn);
#endif

    rhoEqn.solve();

#ifdef OPENFOAMFOUNDATION
        constraints().constrain(rho_);
#elif OPENFOAMESI
        options().correct(rho_);
#endif
}


tmp<fvVectorMatrix> buoyantPimpleFluid::solveUEqn()
{
    // Solve the Momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho_, U())
      + fvm::div(phi(), U())
#ifdef OPENFOAMFOUNDATION
      + turbulence_->divDevTau(U())
     ==
        models().source(rho_, U())
#else
      + turbulence_->divDevRhoReff(U())
     ==
        fvOptions_(rho_, U())
#endif
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
#ifdef OPENFOAMFOUNDATION
    constraints().constrain(UEqn);
#elif OPENFOAMESI
    options().constrain(UEqn);
#endif

    if (pimple().momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf_*fvc::snGrad(rho_)
                  - fvc::snGrad(p_rgh_)
                )*mesh().magSf()
            )
        );
#ifdef OPENFOAMFOUNDATION
        constraints().constrain(U());
#elif OPENFOAMESI
        options().correct(U());
#endif
        K_ = 0.5*magSqr(U());
    }

    return tUEqn;
}


void buoyantPimpleFluid::solveEEqn()
{
    volScalarField& he = thermo_.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, he) + fvm::div(phi(), he)
      + fvc::ddt(rho_, K_) + fvc::div(phi(), K_)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi()/fvc::interpolate(rho_), U()),
                p(),
                "div(phiv,p)"
            )
          : -dpdt_
        )
#ifdef OPENFOAMFOUNDATION
      + thermophysicalTransport_->divq(he)
#else
      - fvm::laplacian(turbulence_->alphaEff(), he)
#endif
     ==
        rho_*(U() & g())
#ifdef OPENFOAMFOUNDATION
      + radiation_->Sh(thermo_, he)
      + models().source(rho_, he)
#elif
      + options()(rho_, he)
#endif
    );

    EEqn.relax();

#ifdef OPENFOAMFOUNDATION
    constraints().constrain(EEqn);
#elif OPENFOAMESI
    options().constrain(EEqn);
#endif

    EEqn.solve();

#ifdef OPENFOAMFOUNDATION
    constraints().constrain(he);
#elif OPENFOAMESI
    options().correct(he);
#endif

    thermo_.correct();
#ifdef OPENFOAMFOUNDATION
    radiation_->correct();
#endif
}


void buoyantPimpleFluid::solvePEqn
(
    const fvVectorMatrix& UEqn
)
{
#ifdef OPENFOAMFOUNDATION
    if (!pimple().simpleRho())
    {
        rho_ = thermo_.rho();
    }
#else
    rho_ = thermo_.rho();
#endif

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi_*p());

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p_rgh_));

    surfaceScalarField phig(-rhorAUf*ghf_*fvc::snGrad(rho_)*mesh().magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::interpolate(rho_)*fvc::flux(HbyA)
          + rhorAUf*fvc::ddtCorr(rho_, U(), phi(), rhoUf_)
        )
      + phig
    );

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh_, rho_, U(), phiHbyA, rhorAUf);

    fvc::makeRelative(phiHbyA, rho_, U());

    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho_) + psi_*correction(fvm::ddt(p_rgh_))
      + fvc::div(phiHbyA)
      ==
#ifdef OPENFOAMFOUNDATION
        models().source(psi_, p_rgh_, rho_.name())
#elif OPENFOAMESI
        options()(psi_, p_rgh_, rho_.name())
#endif
    );

    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh_)
        );

        p_rghEqn.setReference
        (
            pressureControl_.refCell(),
            pressureControl_.refValue()
        );

        p_rghEqn.solve();

        gradp() = fvc::grad(p());

        if (pimple().finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi() = phiHbyA + p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh_.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U() = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
            U().correctBoundaryConditions();
#ifdef OPENFOAMFOUNDATION
            constraints().constrain(U());
#elif OPENFOAMESI
            options().correct(U());
#endif
            K_ = 0.5*magSqr(U());
        }
    }

    p() = p_rgh_ + rho_*gh_;

#ifdef OPENFOAMFOUNDATION
    bool limitedp = !thermo_.incompressible();
#else
    bool limitedp = pressureControl_.limit(p());
#endif

    if (limitedp)
    {
        p_rgh_ = p() - rho_*gh_;
    }

    // Thermodynamic density update
    thermo_.correctRho(psi_*p() - psip0);

    if (limitedp)
    {
        rho_ = thermo_.rho();
    }

    solveRhoEqn();

#ifdef OPENFOAMFOUNDATION
    if (pimple().simpleRho())
    {
        rho_ = thermo_.rho();
    }
#endif

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf_, rho_, U(), phi());

    if (thermo_.dpdt())
    {
        dpdt_ = fvc::ddt(p());

        if (mesh().moving())
        {
            dpdt_ -= fvc::div(fvc::meshPhi(rho_, U()), p());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

buoyantPimpleFluid::buoyantPimpleFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    pThermo_(fluidThermo::New(mesh())),
    thermo_(pThermo_()),
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
        thermo_.rho()
    ),
    psi_(thermo_.psi()),
    dpdt_
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh()
        ),
        mesh(),
        dimensionedScalar("zero", p().dimensions()/dimTime, 0)
    ),
    K_
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.5*magSqr(U())
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("zero", dimLength, 0)
    ),
    ghRef_(-mag(g())*hRef_),
    gh_("gh", (g() & mesh().C()) - ghRef_),
    ghf_("ghf", (g() & mesh().Cf()) - ghRef_),
#ifdef OPENFOAMFOUNDATION
    turbulence_
    (
        compressible::momentumTransportModel::New
        (
            rho_, U(), phi(), thermo_
        )
    ),
    thermophysicalTransport_
    (
        fluidThermophysicalTransportModel::New(turbulence_(), thermo_)
    ),
    radiation_(radiationModel::New(thermo_.T())),
    pressureControl_
    (
        p(),
        p_rgh_,
        pimple().dict(),
        thermo_.incompressible()
    ),
#else
    turbulence_
    (
        compressible::turbulenceModel::New
        (
            rho_, U(), phi(), thermo_
        )
    ),
    pressureControl_(p(), rho_, pimple().dict(), false),
#endif
    cumulativeContErr_(0)
{
    UisRequired();
    pisRequired();

    thermo_.validate("solids4Foam", "h", "e");

    // Reset phi dimensions: compressible
    Info<< "Resetting the dimensions of phi" << endl;
    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);
    phi() = linearInterpolate(rho_*U()) & mesh().Sf();

    mesh().setFluxRequired(p_rgh_.name());

    // Force p_rgh to be consistent with p
    p_rgh_ = p() - rho_*gh_;

    turbulence_->validate();

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
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> buoyantPimpleFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF.ref() =
        (
            mesh().boundary()[patchID].nf()
#ifdef OPENFOAMFOUNDATION
          & (-turbulence_->devTau()().boundaryField()[patchID])
#else
          & (-turbulence_->devRhoReff()().boundaryField()[patchID])
#endif
        );

    return tvF;
}


tmp<scalarField> buoyantPimpleFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF.ref() = p().boundaryField()[patchID];

    return tpF;
}


bool buoyantPimpleFluid::evolve()
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
        pimple().dict().lookupOrDefault("correctPhi", mesh.dynamic())
    );

    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    autoPtr<volScalarField> divrhoU;
    if (correctPhi)
    {
        divrhoU.set
        (
            new volScalarField
            (
                "divrhoU",
                fvc::div(phi())
            )
        );
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), rho_, U());

    // Calculate CourantNo
    fluidModel::CourantNo();

    // Pressure-velocity corrector
    while (pimple().loop())
    {
#ifdef OPENFOAMESI
        if (pimple().firstIter())
#else
        if (pimple().firstPimpleIter())
#endif
        {
            // Store momentum to set rhoUf for introduced faces.
            autoPtr<volVectorField> rhoU;
            if (rhoUf_.valid())
            {
                rhoU.set(new volVectorField("rhoU", rho_*U()));
            }

            if (meshChanged)
            {
                gh_ = (g() & mesh.C()) - ghRef_;
                ghf_ = (g() & mesh.Cf()) - ghRef_;

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi() = mesh.Sf() & rhoUf_();

                    CorrectPhi
                    (
#ifndef OPENFOAMFOUNDATION
                        U(),
#endif
                        phi(),
                        p_rgh_,
                        rho_,
                        psi_,
                        dimensionedScalar("rAUf", dimTime, 1),
                        divrhoU(),
                        pimple()
#ifndef OPENFOAMESIORFOUNDATION
                        ,
                        true
#endif
                    );

                    // Make the fluxes relative to the mesh-motion
                    fvc::makeRelative(phi(), rho_, U());
                }
            }

#ifdef OPENFOAMFOUNDATION
            if (!pimple().simpleRho())
            {
                solveRhoEqn();
            }
#else
            solveRhoEqn();
#endif
        }

        // Momentum equation
        tmp<fvVectorMatrix> tUEqn = solveUEqn();

        // Energy equation
        solveEEqn();

        // --- Pressure corrector loop
        while (pimple().correct())
        {
            solvePEqn(tUEqn());
        }

        compressibleContinuityErrs();

        tUEqn.clear();

        gradU() = fvc::grad(U());

        if (pimple().turbCorr())
        {
            turbulence_->correct();
        }
    }

    rho_ = thermo_.rho();

    Info<< "Fluid temperature min/max(T) = " << min(thermo_.T()).value()
	<< ", " << max(thermo_.T()).value() << " [K]" << endl;

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi(), rho_, U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

#endif // OPENFOAMESIORFOUNDATION

// ************************************************************************* //
