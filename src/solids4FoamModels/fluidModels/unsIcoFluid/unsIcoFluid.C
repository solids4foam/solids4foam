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

#include "unsIcoFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

#include "wedgeFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "EulerDdtScheme.H"
// #include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsIcoFluid, 0);
addToRunTimeSelectionTable(physicsModel, unsIcoFluid, fluid);
addToRunTimeSelectionTable(fluidModel, unsIcoFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void unsIcoFluid::updateRobinFsiInterfaceFlux()
{
    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
             && isA<elasticSlipWallVelocityFvPatchVectorField>
                (
                    U().boundaryField()[patchI]
                )
            )
         || (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
             && isA<elasticWallVelocityFvPatchVectorField>
                (
                    U().boundaryField()[patchI]
                )
            )
        )
        {
            Info<< "Mesh changed: updating flux on Robin interface.";
#ifdef OPENFOAMESIORFOUNDATION
            phi().boundaryFieldRef()[patchI] = 0;
#else
            phi().boundaryField()[patchI] = 0;
#endif
        }
    }
}

void unsIcoFluid::updateRobinFsiInterface()
{
    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
             && isA<elasticSlipWallVelocityFvPatchVectorField>
                (
                    U().boundaryField()[patchI]
                )
            )
         || (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
             && isA<elasticWallVelocityFvPatchVectorField>
                (
                    U().boundaryField()[patchI]
                )
            )
        )
        {
#ifdef OPENFOAMESIORFOUNDATION
            word ddtScheme
            (
                mesh().ddtScheme("ddt(" + U().name() +')')
            );
#else
            word ddtScheme
            (
                mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
            );
#endif

            if
            (
                ddtScheme
             == fv::EulerDdtScheme<vector>::typeName
            )
            {
#ifdef OPENFOAMESIORFOUNDATION
                phi().boundaryFieldRef()[patchI] =
                    phi().oldTime().boundaryField()[patchI];

                rAUf_.boundaryFieldRef()[patchI] =
                    runTime().deltaT().value();
#else
                phi().boundaryField()[patchI] =
                    phi().oldTime().boundaryField()[patchI];

                rAUf_.boundaryField()[patchI] =
                    runTime().deltaT().value();
#endif
        }
            else if
            (
                ddtScheme
             == fv::backwardDdtScheme<vector>::typeName
            )
            {
                if(runTime().timeIndex() == 1)
                {
#ifdef OPENFOAMESIORFOUNDATION
                    phi().boundaryFieldRef()[patchI] =
                        phi().oldTime().boundaryField()[patchI];

                    rAUf_.boundaryFieldRef()[patchI] =
                        runTime().deltaT().value();
#else
                    phi().boundaryField()[patchI] =
                        phi().oldTime().boundaryField()[patchI];

                    rAUf_.boundaryField()[patchI] =
                        runTime().deltaT().value();
#endif

                    phi().oldTime().oldTime();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

#ifdef OPENFOAMESIORFOUNDATION
                    phi().boundaryFieldRef()[patchI] =
#else
                    phi().boundaryField()[patchI] =
#endif
                        (Co/Cn)*phi().oldTime().boundaryField()
                        [
                            patchI
                        ]
                      - (Coo/Cn)*phi().oldTime().oldTime().boundaryField()
                        [
                            patchI
                        ];

#ifdef OPENFOAMESIORFOUNDATION
                    rAUf_.boundaryFieldRef()[patchI] = deltaT/Cn;
#else
                    rAUf_.boundaryField()[patchI] = deltaT/Cn;
#endif
                }
            }
        }
    }
}

void unsIcoFluid::correctRobinFsiInterfaceFlux()
{
    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
             && isA<elasticSlipWallVelocityFvPatchVectorField>
                (
                    U().boundaryField()[patchI]
                )
            )
         || (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p().boundaryField()[patchI]
                )
             && isA<elasticWallVelocityFvPatchVectorField>
                (
                    U().boundaryField()[patchI]
                )
            )
        )
        {
#ifdef OPENFOAMESIORFOUNDATION
            word ddtScheme
            (
                mesh().ddtScheme("ddt(" + U().name() +')')
            );
#else
            word ddtScheme
            (
                mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
            );
#endif

            if
            (
                ddtScheme
             == fv::EulerDdtScheme<vector>::typeName
            )
            {
#ifdef OPENFOAMESIORFOUNDATION
                phi().boundaryFieldRef()[patchI] =
#else
                phi().boundaryField()[patchI] =
#endif
                    phi().oldTime().boundaryField()[patchI]
                  - p().boundaryField()
                    [
                        patchI
                    ].snGrad()*runTime().deltaT().value()
                   *mesh().magSf().boundaryField()
                    [
                        patchI
                    ];
            }
            else if
            (
                ddtScheme
             == fv::backwardDdtScheme<vector>::typeName
            )
            {
                if(runTime().timeIndex() == 1)
                {
#ifdef OPENFOAMESIORFOUNDATION
                    phi().boundaryFieldRef()[patchI] =
                        phi().oldTime().boundaryField()[patchI];

                    rAUf_.boundaryFieldRef()[patchI] =
                        runTime().deltaT().value();
#else
                    phi().boundaryField()[patchI] =
                        phi().oldTime().boundaryField()[patchI];

                    rAUf_.boundaryField()[patchI] =
                        runTime().deltaT().value();
#endif

                    phi().oldTime().oldTime();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

#ifdef OPENFOAMESIORFOUNDATION
                    phi().boundaryFieldRef()[patchI] =
#else
                    phi().boundaryField()[patchI] =
#endif
                        (Co/Cn)*phi().oldTime().boundaryField()
                        [
                            patchI
                        ]
                      - (Coo/Cn)*phi().oldTime().oldTime().boundaryField()
                        [
                            patchI
                        ];

#ifdef OPENFOAMESIORFOUNDATION
                    rAUf_.boundaryFieldRef()[patchI] = deltaT/Cn;
#else
                    rAUf_.boundaryField()[patchI] = deltaT/Cn;
#endif
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsIcoFluid::unsIcoFluid
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
    rho_(transportProperties_.lookup("rho")),
    ddtU_(fvc::ddt(U())),
    Uf_
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
    ),
    rAUf_
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        calculatedFvPatchScalarField::typeName
    )
{
    UisRequired();
    pisRequired();

    phi().oldTime();
    Uf_.oldTime();

#ifdef OPENFOAMESIORFOUNDATION
    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + U().name() +')')
    );
#else
    word ddtScheme
    (
        mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
    );
#endif

    if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        Uf_.oldTime().oldTime();
    }

    ddtU_.checkIn();

#ifdef OPENFOAMESIORFOUNDATION
    mesh().setFluxRequired(p().name());
#else
    mesh().schemesDict().setFluxRequired(p().name());
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> unsIcoFluid::patchViscousForce(const label patchID) const
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


tmp<scalarField> unsIcoFluid::patchPressureForce(const label patchID) const
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

void unsIcoFluid::evolveConsistent()
{
    Info<< "Evolving fluid model with consistent strategy: "
                << this->type() << endl;

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

    if (meshChanged)
    {
        updateRobinFsiInterfaceFlux();
    }

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
    //
        volScalarField AU("AU", HUEqn.A() + ddtUEqn.A());

#   include "correctCentralCoeffs.H"

#if FOAMEXTEND > 40
    // Prepare clean 1/a_p without time derivative contribution
    const volScalarField rAU = 1.0/HUEqn.A();

    rAUf_ = fvc::interpolate(rAU);
#else
    const volScalarField rAU(1.0/AU);

    rAUf_ = 1.0/fvc::interpolate(AUcorr, "interpolate(U)");
#endif

    while (piso().correct())
    {
#if FOAMEXTEND > 40
        // Calculate U from convection-diffusion matrix
        U() = rAU*HUEqn.H();

        // Consistently calculate flux
        piso().calcTransientConsistentFlux(phi(), U(), rAU, ddtUEqn);
#else
#               include "calcPhiConsistent.H"

        // Calculate U from convection-diffusion matrix
        U() = rAU*(HUEqn.H() + ddtUEqn.H());
#endif

        updateRobinFsiInterface();

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
                    rAUf_, p(), "laplacian((1|A(U)),p)"
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
                        if (Pstream::parRun())
                        {
                                // Needed for extrapolated BC
                                p().correctBoundaryConditions();
                        }

            gradp() = fvc::grad(p());

            if (piso().finalNonOrthogonalIter())
            {
                phi() -= pEqn.flux();
            }
        }

        correctRobinFsiInterfaceFlux();

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
                ddtU_ = fvc::ddt(U());
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

        // Update velocity on faces
    Uf_ = fvc::interpolate(U(), "interpolate(U)");

        // Get tangential velocity
    Uf_ -= (mesh.Sf() & Uf_)*mesh.Sf()/magSqr(mesh.Sf());

        // Update with normal from flux
    Uf_ += phi()*mesh.Sf()/magSqr(mesh.Sf());
}

void unsIcoFluid::evolveInconsistent()
{
    Info<< "Evolving fluid model with inconsistent strategy: "
                << this->type() << endl;

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

    if (meshChanged)
    {
        updateRobinFsiInterfaceFlux();
    }

    // CourantNo
    fluidModel::CourantNo();

    // Prepare for the pressure solution
    label  pRefCell  = 0;
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
#endif

    rAUf_ = fvc::interpolate(rAU);

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

        updateRobinFsiInterface();

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
                    rAUf_, p(), "laplacian((1|A(U)),p)"
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
                        if (Pstream::parRun())
                        {
                                // Needed for extrapolated BC
                                p().correctBoundaryConditions();
                        }

            gradp() = fvc::grad(p());

            if (piso().finalNonOrthogonalIter())
            {
                phi() -= pEqn.flux();
            }
        }

        correctRobinFsiInterfaceFlux();

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
}

bool unsIcoFluid::evolve()
{
    Switch consistent = true;
    if (fluidProperties().found("consistent"))
    {
        consistent = Switch(fluidProperties().lookup("consistent"));
    }

    if (consistent)
    {
        evolveConsistent();
    }
    else
    {
        evolveInconsistent();
    }

    return 0;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
