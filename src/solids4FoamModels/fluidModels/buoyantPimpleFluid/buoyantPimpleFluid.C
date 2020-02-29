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
      ==
        fvOptions_(rho_)
    );

    fvOptions_.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions_.correct(rho_);
}


tmp<fvVectorMatrix> buoyantPimpleFluid::solveUEqn()
{
    // Solve the Momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho_, U())
      + fvm::div(phi(), U())
      + turbulence_->divDevRhoReff(U())
     ==
        fvOptions_(rho_, U())
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions_.constrain(UEqn);

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

        fvOptions_.correct(U());
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
      - fvm::laplacian(turbulence_->alphaEff(), he)
     ==
        rho_*(U() & g_)
#ifdef OPENFOAMFOUNDATION
      + radiation_->Sh(thermo_, he)
#endif
      + fvOptions_(rho_, he)
    );

    EEqn.relax();

    fvOptions_.constrain(EEqn);

    EEqn.solve();

    fvOptions_.correct(he);

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
        fvOptions_(psi_, p_rgh_, rho_.name())
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
            fvOptions_.correct(U());
            K_ = 0.5*magSqr(U());
        }
    }

    p() = p_rgh_ + rho_*gh_;

    bool limitedp = pressureControl_.limit(p());

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
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
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
    ghRef_(-mag(g_)*hRef_),
    gh_("gh", (g_ & mesh().C()) - ghRef_),
    ghf_("ghf", (g_ & mesh().Cf()) - ghRef_),
    fvOptions_(fv::options::New(mesh())),
    turbulence_
    (
        compressible::turbulenceModel::New
        (
            rho_, U(), phi(), thermo_
        )
    ),
#ifdef OPENFOAMFOUNDATION
    radiation_(radiationModel::New(thermo_.T())),
#endif
    pressureControl_(p(), rho_, pimple().dict(), false),
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

        rhoUf_ = new surfaceVectorField
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
        );
    }

    // Check if any finite volume option is present
    if (!fvOptions_.optionList::size())
    {
        Info << "No finite volume options present\n" << endl;
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
          & (-turbulence_->devRhoReff()().boundaryField()[patchID])
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
        divrhoU = new volScalarField
        (
            "divrhoU",
            fvc::div(phi())
        );
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), rho_, U());

    // Calculate CourantNo
    fluidModel::CourantNo();

    // Pressure-velocity corrector
    while (pimple().loop())
    {
        if (pimple().firstPimpleIter())
        {
            // Store momentum to set rhoUf for introduced faces.
            autoPtr<volVectorField> rhoU;
            if (rhoUf_.valid())
            {
                rhoU = new volVectorField("rhoU", rho_*U());
            }

            if (meshChanged)
            {
                gh_ = (g_ & mesh.C()) - ghRef_;
                ghf_ = (g_ & mesh.Cf()) - ghRef_;

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi() = mesh.Sf() & rhoUf_();

                    CorrectPhi
                    (
                        U(),
                        phi(),
                        p_rgh_,
                        rho_,
                        psi_,
                        dimensionedScalar("rAUf", dimTime, 1),
                        divrhoU(),
                        pimple(),
                        true
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
