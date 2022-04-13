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

#include "icoFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(icoFluid, 0);
addToRunTimeSelectionTable(fluidModel, icoFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

icoFluid::icoFluid
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
    nu_(transportProperties_.lookup("nu")),
    rho_(transportProperties_.lookup("rho"))
{
    UisRequired();
    pisRequired();
    
#ifdef OPENFOAMESIORFOUNDATION
    mesh().setFluxRequired(p().name());
#else
    mesh().schemesDict().setFluxRequired(p().name());
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> icoFluid::patchViscousForce(const label patchID) const
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
        rho_.value()*nu_.value()*U().boundaryField()[patchID].snGrad();

    return tvF;
}


tmp<scalarField> icoFluid::patchPressureForce(const label patchID) const
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


bool icoFluid::evolve()
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

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), piso().dict(), pRefCell, pRefValue);
    
    // Time-derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U()));

    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi(), U())
      - fvm::laplacian(nu_, U())
    );

    if (piso().momentumPredictor())
    {
        solve(ddtUEqn + HUEqn == -fvc::grad(p()));
    }

    // --- PISO loop

#if FOAMEXTEND > 40
    // Prepare clean 1/a_p without time derivative contribution
    volScalarField rAU(1.0/HUEqn.A());
#else
    volScalarField rAU(1.0/(HUEqn.A() + ddtUEqn.A()));
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
#endif

    while (piso().correct())
    {
#if FOAMEXTEND > 40
        // Calculate U from convection-diffusion matrix
        U() = rAU*HUEqn.H();

        // Consistently calculate flux
        piso().calcTransientConsistentFlux(phi(), U(), rAU, ddtUEqn);
#else
        // Calculate U from convection-diffusion matrix
        U() = rAU*(HUEqn.H() + ddtUEqn.H());

        phi() = (fvc::interpolate(U()) & mesh.Sf());
#endif

        adjustPhi(phi(), U(), p());

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
                    rAUf, p(), "laplacian((1|A(U)),p)"
#endif
                )
             == fvc::div(phi())
            );

            pEqn.setReference(pRefCell, pRefValue);
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
                phi() -= pEqn.flux();
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

        U() -= rAU*gradp();
        U().correctBoundaryConditions();
#endif

        gradU() = fvc::grad(U());
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());
    
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
