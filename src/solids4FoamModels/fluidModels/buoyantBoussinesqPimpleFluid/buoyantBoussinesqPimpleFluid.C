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

#include "buoyantBoussinesqPimpleFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "constrainHbyA.H"
    #include "constrainPressure.H"
#else
    #include "fvc.H"
    #include "fvm.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(buoyantBoussinesqPimpleFluid, 0);
addToRunTimeSelectionTable(physicsModel, buoyantBoussinesqPimpleFluid, fluid);
addToRunTimeSelectionTable(fluidModel, buoyantBoussinesqPimpleFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> buoyantBoussinesqPimpleFluid::solveUEqn()
{
    // Solve the momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U())
      + fvm::div(phi(), U())
#ifdef OPENFOAMFOUNDATION
      + turbulence_->divDevSigma(U())
     ==
        models().source(U())
#elif OPENFOAMESI
      + turbulence_->divDevReff(U()))
     ==
        options()(U())
#else
      + turbulence_->divDevReff()
#endif
    );
#ifdef OPENFOAMESIORFOUNDATION
    fvVectorMatrix& UEqn = tUEqn.ref();
#else
    fvVectorMatrix& UEqn = tUEqn();
#endif

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
                  - ghf_*fvc::snGrad(rhok_)
                  - fvc::snGrad(p_rgh_)
                )*mesh().magSf()
            )
        );

#ifdef OPENFOAMFOUNDATION
    constraints().constrain(U());
#elif OPENFOAMESI
    options().correct(U());
#endif
    }

    return tUEqn;
}


void buoyantBoussinesqPimpleFluid::solveTEqn()
{
    alphaEff_ = turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_;
    alphaEff_.correctBoundaryConditions();

    fvScalarMatrix TEqn
    (
        fvm::ddt(T_)
      + fvm::div(phi(), T_)
      - fvm::laplacian(alphaEff_, T_)
#ifdef OPENFOAMFOUNDATION
     ==
        radiation_->ST(rhoCpRef_, T_)
      + models().source(T_)
#elif OPENFOAMESI
     ==
        options()(T_)
#endif
    );

    TEqn.relax();

#ifdef OPENFOAMFOUNDATION
    constraints().constrain(TEqn);
#elif OPENFOAMESI
    options().constrain(TEqn);
#endif


    TEqn.solve();

#ifdef OPENFOAMFOUNDATION
    radiation_->correct();
    constraints().constrain(T_);
#elif OPENFOAMESI
    options().correct(T_);
#endif
    rhok_ = 1.0 - beta_*(T_ - TRef_);
}


void buoyantBoussinesqPimpleFluid::solvePEqn
(
    const fvVectorMatrix& UEqn
)
{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

#ifdef OPENFOAMESIORFOUNDATION
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p_rgh_));
#else
    volVectorField HbyA("HbyA", U());
    HbyA = rAU*UEqn.H();
#endif

    surfaceScalarField phig(-rAUf*ghf_*fvc::snGrad(rhok_)*mesh().magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
#ifdef OPENFOAMESIORFOUNDATION
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U(), phi())
#else
        (fvc::interpolate(HbyA) & mesh().Sf())
      + fvc::ddtPhiCorr(rAU, U(), phi())
#endif
      + phig
    );

#ifdef OPENFOAMESIORFOUNDATION
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh_, U(), phiHbyA, rAUf);
#endif

    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh_) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell_, getRefCellValue(p_rgh_, pRefCell_));

        p_rghEqn.solve
        (
#ifndef OPENFOAMESIORFOUNDATION
            mesh().solutionDict().solver
            (
                p_rgh_.select(pimple().finalInnerIter())
            )
#endif
        );

        gradp() = fvc::grad(p());

        if (pimple().finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi() = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh_.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U() = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U().correctBoundaryConditions();
#ifdef OPENFOAMFOUNDATION
            constraints().constrain(U());
#elif OPENFOAMESI
            options().correct(U());
#endif
        }
    }

    fluidModel::continuityErrs();

    p() = p_rgh_ + rhok_*gh_;

    if (p_rgh_.needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pRefCell_)
        );
        p_rgh_ = p() - rhok_*gh_;
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

buoyantBoussinesqPimpleFluid::buoyantBoussinesqPimpleFluid
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
    T_
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    laminarTransport_(U(), phi()),
    rho_(laminarTransport_.lookup("rho")),
    beta_(laminarTransport_.lookup("beta")),
    TRef_(laminarTransport_.lookup("TRef")),
    Pr_(laminarTransport_.lookup("Pr")),
    Prt_(laminarTransport_.lookup("Prt")),
#ifdef OPENFOAMFOUNDATION
    turbulence_
    (
        incompressible::momentumTransportModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
#else
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
#endif
    alphaEff_
    (
        IOobject
        (
            "alphaEff",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_
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
    rhok_
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh()
        ),
        1.0 - beta_*(T_ - TRef_)
    ),
#ifdef OPENFOAMFOUNDATION
    radiation_(radiationModel::New(T_)),
    rhoCpRef_
    (
        "rhoCpRef",
        dimDensity*dimEnergy/dimMass/dimTemperature,
        1.0
    ),
#endif
    pRefCell_(0),
    pRefValue_(0)
{
    UisRequired();

    // Reset p dimensions
    Info<< "Resetting the dimensions of p and its gradient" << endl;
    p().dimensions().reset(dimPressure/dimDensity);
    gradp().dimensions().reset(dimPressure/dimDensity/dimLength);
    p() = p_rgh_ + rhok_*gh_;

#ifdef OPENFOAMESIORFOUNDATION
    turbulence_->validate();

    setRefCell(p(), p_rgh_, pimple().dict(), pRefCell_, pRefValue_);
    mesh().setFluxRequired(p_rgh_.name());

    #ifdef OPENFOAMFOUNDATION
    if (!isType<radiationModels::noRadiation>(radiation_()))
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false // Do not register
            )
        );

        dimensionedScalar rhoRef
        (
            "rhoRef",
            dimDensity,
            transportProperties
        );

        dimensionedScalar CpRef
        (
            "CpRef",
            dimSpecificHeatCapacity,
            transportProperties
        );

        rhoCpRef_ = rhoRef*CpRef;
    }
    #endif
#else
    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p_rgh_.name());
#endif

    if (p_rgh_.needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pRefCell_)
        );
        p_rgh_ = p() - rhok_*gh_;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> buoyantBoussinesqPimpleFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

#ifdef OPENFOAMESIORFOUNDATION
    tvF.ref() =
#else
    tvF() =
#endif
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
#ifdef OPENFOAMFOUNDATION
          & (-turbulence_->devSigma()().boundaryField()[patchID])
#else
          & (-turbulence_->devReff()().boundaryField()[patchID])
#endif
        );

    return tvF;
}


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAMESIORFOUNDATION
    tpF.ref() =
#else
    tpF() =
#endif
        rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


bool buoyantBoussinesqPimpleFluid::evolve()
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

    // Update gh fields as the mesh may have moved
    gh_ = (g() & mesh.C()) - ghRef_;
    ghf_ = (g() & mesh.Cf()) - ghRef_;

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // Calculate CourantNo
    CourantNo();

    // Pressure-velocity corrector
    while (pimple().loop())
    {
        // Momentum equation
        tmp<fvVectorMatrix> tUEqn = solveUEqn();

        // Temperature equation
        solveTEqn();

        // --- Pressure corrector loop
        while (pimple().correct())
        {
#ifdef OPENFOAMESIORFOUNDATION
            solvePEqn(tUEqn.ref());
#else
            solvePEqn(tUEqn());
#endif
        }

        tUEqn.clear();

        gradU() = fvc::grad(U());

#ifdef OPENFOAMESIORFOUNDATION
        laminarTransport_.correct();
#endif
        turbulence_->correct();
    }

    Info<< "Fluid temperature min/max(T) = " << min(T_).value()
	<< ", " << max(T_).value() << " [K]" << endl;

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
