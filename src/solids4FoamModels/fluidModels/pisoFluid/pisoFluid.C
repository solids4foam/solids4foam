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

#include "pisoFluid.H"
#include "fvm.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "constrainHbyA.H"
    #include "constrainPressure.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pisoFluid, 0);
addToRunTimeSelectionTable(physicsModel, pisoFluid, fluid);
addToRunTimeSelectionTable(fluidModel, pisoFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pisoFluid::pisoFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    laminarTransport_(U(), phi()),
    turbulence_
    (
#ifdef OPENFOAMFOUNDATION
        incompressible::momentumTransportModel::New
#else
        incompressible::turbulenceModel::New
#endif
        (
            U(), phi(), laminarTransport_
        )
    ),
    rho_
    (
        laminarTransport_.lookup("rho")
    ),
#ifdef OPENFOAMESI
    fvOptions_(fv::options::New(mesh())),
#endif
    pRefCell_(0),
    pRefValue_(0)
{
    UisRequired();
    pisRequired();
    
    setRefCell(p(), piso().dict(), pRefCell_, pRefValue_);

#ifdef OPENFOAMESIORFOUNDATION
    turbulence_->validate();

    mesh().setFluxRequired(p().name());
#else
    mesh().schemesDict().setFluxRequired(p().name());
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> pisoFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

#ifdef FOAMEXTEND
    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );
#elif OPENFOAMESI
    tvF.ref() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );
#else
    tvF.ref() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devSigma()().boundaryField()[patchID])
        );
#endif

    return tvF;
}


tmp<scalarField> pisoFluid::patchPressureForce(const label patchID) const
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


bool pisoFluid::evolve()
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
    
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // CourantNo
    fluidModel::CourantNo();

    // Time-derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U()));

    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi(), U())
#ifdef FOAMEXTEND
      + turbulence_->divDevReff()
#elif OPENFOAMESI
      + turbulence_->divDevReff(U())
     ==
        fvOptions_(U())
#else
      + turbulence_->divDevSigma(U())
#endif
    );

#if (defined(OPENFOAMESIORFOUNDATION) || FOAMEXTEND < 41)
    fvVectorMatrix UEqn(ddtUEqn + HUEqn);
#endif

    UEqn.relax();

#ifdef OPENFOAMESI
    fvOptions_.constrain(UEqn);
#endif

    if (piso().momentumPredictor())
    {
#if (defined(OPENFOAMESIORFOUNDATION) || FOAMEXTEND < 41)
        solve(UEqn == -fvc::grad(p()));
#else
        solve(ddtUEqn + HUEqn == -fvc::grad(p()));
#endif

#ifdef OPENFOAMESI
        fvOptions_.correct(U());
#endif
    }

    // --- PISO loop
    while (piso().correct())
    {
#if FOAMEXTEND > 40
        // Prepare clean 1/a_p without time derivative contribution
        volScalarField rAU(1.0/HUEqn.A());
#else
        volScalarField rAU(1.0/UEqn.A());
#ifndef OPENFOAMESIORFOUNDATION
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
#endif
#endif

#if FOAMEXTEND > 40
        // Calculate U from convection-diffusion matrix
        U() = rAU*HUEqn.H();

        // Consistently calculate flux
        piso().calcTransientConsistentFlux(phi(), U(), rAU, ddtUEqn);
#else
#ifdef OPENFOAMESIORFOUNDATION
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p()));

        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + fvc::interpolate(rAU)*fvc::ddtCorr(U(), phi())
        );
#else
        // Calculate U from convection-diffusion matrix
        U() = rAU*UEqn.H();

        phi() = (fvc::interpolate(U()) & mesh.Sf());
#endif
#endif

#ifdef OPENFOAMESIORFOUNDATION
        adjustPhi(phiHbyA, U(), p());

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p(), U(), phiHbyA, rAU);
#else
        adjustPhi(phi(), U(), p());
#endif

        // Non-orthogonal pressure corrector loop
        while (piso().correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
#if FOAMEXTEND > 40
                    fvc::interpolate(rAU)/piso().aCoeff(U().name()),
                    p(),
                    "laplacian(rAU," + p().name() + ')'
#else
#ifdef OPENFOAMESIORFOUNDATION
                    rAU, p()
#else
                    rAUf, p(), "laplacian((1|A(U)),p)"
#endif
#endif
                )
#ifdef OPENFOAMESIORFOUNDATION
             == fvc::div(phiHbyA)
#else
             == fvc::div(phi())
#endif
            );

            pEqn.setReference(pRefCell_, pRefValue_);

#ifdef OPENFOAMESIORFOUNDATION
            pEqn.solve();
#else
            pEqn.solve
            (
                mesh.solutionDict().solver(p().select(piso().finalInnerIter()))
            );
#endif

            gradp() = fvc::grad(p());

            if (piso().finalNonOrthogonalIter())
            {
#ifdef OPENFOAMESIORFOUNDATION
                phi() = phiHbyA - pEqn.flux();
#else
                phi() -= pEqn.flux();
#endif
            }
        }

        fluidModel::continuityErrs();

#if FOAMEXTEND > 40
        // Consistently reconstruct velocity after pressure equation.
        // Note: flux is made relative inside the function
        piso().reconstructTransientVelocity(U(), phi(), ddtUEqn, rAU, p());
#else
        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi(), U());

    #ifdef OPENFOAMESIORFOUNDATION
        U() = HbyA - rAU*fvc::grad(p());
        U().correctBoundaryConditions();
    #else
        U() -= rAU*gradp();
        U().correctBoundaryConditions();
    #endif

    #ifdef OPENFOAMESI
        fvOptions_.correct(U());
    #endif

#endif
        gradU() = fvc::grad(U());
    }

#ifdef OPENFOAMESIORFOUNDATION
    laminarTransport_.correct();
#endif
    turbulence_->correct();

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
