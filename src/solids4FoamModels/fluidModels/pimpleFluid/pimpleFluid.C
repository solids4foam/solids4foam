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

#include "pimpleFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "CorrectPhi.H"
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

defineTypeNameAndDebug(pimpleFluid, 0);
addToRunTimeSelectionTable(physicsModel, pimpleFluid, fluid);
addToRunTimeSelectionTable(fluidModel, pimpleFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

#ifdef OPENFOAMESIORFOUNDATION
void pimpleFluid::correctPhi()
{
    CorrectPhi
    (
        U(),
        phi(),
        p(),
        dimensionedScalar("rAUf", dimTime, 1),
        geometricZeroField(),
        pimple(),
        true
    );

    fluidModel::continuityErrs();
}
#endif

#ifdef OPENFOAMESIORFOUNDATION
tmp<fvVectorMatrix> pimpleFluid::solveUEqn(tmp<fvVectorMatrix> totalUEqn)
{
    tmp<fvVectorMatrix> tUEqn(totalUEqn);
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions_.constrain(UEqn);

    if (pimple().momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p()));

        fvOptions_.correct(U());
    }

    return tUEqn;
}
#endif

#ifdef OPENFOAMESIORFOUNDATION
void pimpleFluid::solvePEqn
(
    const fvVectorMatrix& UEqn
)
{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p()));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rAU)*fvc::ddtCorr(U(), phi(), Uf_)
    );

    if (p().needReference())
    {
        fvc::makeRelative(phiHbyA, U());
        adjustPhi(phiHbyA, U(), p());
        fvc::makeAbsolute(phiHbyA, U());
    }

    tmp<volScalarField> rAtU(rAU);

    if (pimple().consistent())
    {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p())*mesh().magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p());
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p(), U(), phiHbyA, rAtU());

    // Non-orthogonal pressure corrector loop
    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p()) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell_, pRefValue_);

        pEqn.solve();

        gradp() = fvc::grad(p());

        if (pimple().finalNonOrthogonalIter())
        {
            phi() = phiHbyA - pEqn.flux();
        }
    }

    fluidModel::continuityErrs();

    // Explicitly relax pressure for momentum corrector
    p().relax();

    U() = HbyA - rAtU*fvc::grad(p());
    U().correctBoundaryConditions();
    fvOptions_.correct(U());

    gradU() = fvc::grad(U());

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf_, U(), phi());

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());
}
#endif

#ifndef OPENFOAMESIORFOUNDATION
void pimpleFluid::solveUEqn
(
    const scalar& UUrf,
    const fvVectorMatrix& ddtUEqn,
    const fvVectorMatrix& HUEqn
)
{
    if (pimple().momentumPredictor())
    {
#if FOAMEXTEND > 40
        solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p()));
#else
        solve
        (
            ddtUEqn
          + relax(HUEqn, UUrf)
         ==
          - fvc::grad(p()),
            mesh().solutionDict().solver((U().select(pimple().finalIter())))
        );
#endif
    }
    else
    {
        // Explicit update
        U() =
            (ddtUEqn.H() + HUEqn.H() - fvc::grad(p()))/(HUEqn.A() + ddtUEqn.A());
        U().correctBoundaryConditions();
    }
}
#endif

#ifndef OPENFOAMESIORFOUNDATION
void pimpleFluid::solvePEqn
(
    const scalar& UUrf,
    const fvVectorMatrix& ddtUEqn,
    const fvVectorMatrix& HUEqn
)
{
    p().boundaryField().updateCoeffs();

#if FOAMEXTEND > 40
    // Prepare clean 1/a_p without time derivative and under-relaxation
    // contribution
    rAU_ = 1.0/HUEqn.A();

    // Calculate U from convection-diffusion matrix
    U() = rAU_*HUEqn.H();

    // Consistently calculate flux
    pimple().calcTransientConsistentFlux(phi(), U(), rAU_, ddtUEqn);
#else
    // Prepare clean Ap without time derivative contribution and
    // without contribution from under-relaxation
    // HJ, 26/Oct/2015
    aU_ = HUEqn.A();

    // Store velocity under-relaxation point before using U for the flux
    // precursor
    U().storePrevIter();

    U() = HUEqn.H()/aU_;
    phi() = (fvc::interpolate(U()) & mesh().Sf());
#endif

    adjustPhi(phi(), U(), p());

    // Non-orthogonal pressure corrector loop
    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
#if FOAMEXTEND > 40
                fvc::interpolate(rAU_)/pimple().aCoeff(U().name()),
                p(),
                "laplacian(rAU," + p().name() + ')'
#else
                1/aU_, p(), "laplacian((1|A(U)),p)"
#endif
            )
         == fvc::div(phi())
        );

        pEqn.setReference(pRefCell_, pRefValue_);
        pEqn.solve
        (
            mesh().solutionDict().solver
            (
                p().select(pimple().finalInnerIter())
            )
        );

        gradp() = fvc::grad(p());

        if (pimple().finalNonOrthogonalIter())
        {
            phi() -= pEqn.flux();
        }
    }

    fluidModel::continuityErrs();

    // Explicitly relax pressure for momentum corrector
    // except for last corrector
    if (!pimple().finalIter())
    {
        p().relax();
    }

#if FOAMEXTEND > 40
    // Consistently reconstruct velocity after pressure equation.
    // Note: flux is made relative inside the function
    pimple().reconstructTransientVelocity(U(), phi(), ddtUEqn, rAU_, p());
#else
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    U() = UUrf*
        (
            1.0/(aU_ + ddtUEqn.A())*
            (
                U()*aU_ - fvc::grad(p()) + ddtUEqn.H()
            )
        )
      + (1 - UUrf)*U().prevIter();
    U().correctBoundaryConditions();
#endif

    gradU() = fvc::grad(U());
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pimpleFluid::pimpleFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    laminarTransport_(U(), phi()),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
    rho_(laminarTransport_.lookup("rho")),
#ifdef OPENFOAMESIORFOUNDATION
    fvOptions_(fv::options::New(mesh())),
#else
#if FOAMEXTEND > 40
    rAU_
#else
    aU_
#endif
    (
        IOobject
        (
#if FOAMEXTEND > 40
            "rAU",
#else
            "aU",
#endif
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
#if FOAMEXTEND > 40
        runTime.deltaT(),
#else
        1/runTime.deltaT(),
#endif
        zeroGradientFvPatchScalarField::typeName
    ),
#endif
    pRefCell_(0),
    pRefValue_(0)
{
    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);

#ifdef OPENFOAMESIORFOUNDATION
    mesh().setFluxRequired(p().name());

    turbulence_->validate();

    if (mesh().dynamic())
    {
        Info<< "Constructing face velocity Uf\n" << endl;

        Uf_ = new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U())
        );
    }
#else
    mesh().schemesDict().setFluxRequired(p().name());
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> pimpleFluid::patchViscousForce(const label patchID) const
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
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> pimpleFluid::patchPressureForce(const label patchID) const
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


bool pimpleFluid::evolve()
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

#ifdef OPENFOAMESIORFOUNDATION
    if (meshChanged)
    {
        // Calculate absolute flux
        // from the mapped surface velocity
        phi() = mesh.Sf() & Uf_();

        correctPhi();
    }
#endif
    
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // CourantNo
    fluidModel::CourantNo();

    // --- PIMPLE loop
    while (pimple().loop())
    {
        // Time-derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(U()));

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(phi(), U())
#ifndef OPENFOAMESIORFOUNDATION
          + turbulence_->divDevReff()
#else
          + turbulence_->divDevReff(U())
         ==
            fvOptions_(U())
#endif
        );

#ifdef OPENFOAMESIORFOUNDATION
        tmp<fvVectorMatrix> tUEqn = solveUEqn(ddtUEqn + HUEqn);
#else
#if FOAMEXTEND < 41
        // Get under-relaxation factor
        scalar UUrf =
            mesh.solutionDict().equationRelaxationFactor
            (
                U().select(pimple().finalIter())
            );

        solveUEqn(UUrf, ddtUEqn, HUEqn);
#else
        solveUEqn(scalar(1), ddtUEqn, HUEqn);
#endif
#endif

        // --- PISO loop
        while (pimple().correct())
        {
#ifdef OPENFOAMESIORFOUNDATION
            solvePEqn(tUEqn.ref());
#else
#if FOAMEXTEND < 41
            solvePEqn(UUrf, ddtUEqn, HUEqn);
#else
            solvePEqn(scalar(1), ddtUEqn, HUEqn);
#endif
#endif
        }

#ifdef OPENFOAMESIORFOUNDATION
        laminarTransport_.correct();
#endif
        turbulence_->correct();
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
