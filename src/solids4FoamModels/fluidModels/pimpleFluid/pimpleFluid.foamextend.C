/*---------------------------------------------------------------------------*\
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

#include "pimpleFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "EulerDdtScheme.H"
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

defineTypeNameAndDebug(pimpleFluid, 0);
addToRunTimeSelectionTable(fluidModel, pimpleFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void pimpleFluid::updateRobinFsiInterfaceFlux()
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
            phi().boundaryField()[patchI] = 0;
        }
    }
}

void pimpleFluid::updateRobinFsiInterface()
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
            word ddtScheme
            (
                mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
            );

            if
            (
                ddtScheme
             == fv::EulerDdtScheme<vector>::typeName
            )
            {
                phi().boundaryField()[patchI] =
                    phi().oldTime().boundaryField()[patchI];

                rAUf_.boundaryField()[patchI] =
                    runTime().deltaT().value();
            }
            else if
            (
                ddtScheme
             == fv::backwardDdtScheme<vector>::typeName
            )
            {
                if(runTime().timeIndex() == 1)
                {
                    phi().boundaryField()[patchI] =
                        phi().oldTime().boundaryField()[patchI];

                    rAUf_.boundaryField()[patchI] =
                        runTime().deltaT().value();

                    phi().oldTime().oldTime();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

                    phi().boundaryField()[patchI] =
                        (Co/Cn)*phi().oldTime().boundaryField()
                        [
                            patchI
                        ]
                      - (Coo/Cn)*phi().oldTime().oldTime().boundaryField()
                        [
                            patchI
                        ];

                    rAUf_.boundaryField()[patchI] = deltaT/Cn;
                }
            }
        }
    }
}

void pimpleFluid::correctRobinFsiInterfaceFlux()
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
            word ddtScheme
            (
                mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
            );

            if
            (
                ddtScheme
             == fv::EulerDdtScheme<vector>::typeName
            )
            {
                phi().boundaryField()[patchI] =
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
                    phi().boundaryField()[patchI] =
                        phi().oldTime().boundaryField()[patchI];

                    rAUf_.boundaryField()[patchI] =
                        runTime().deltaT().value();

                    phi().oldTime().oldTime();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

                    phi().boundaryField()[patchI] =
                        (Co/Cn)*phi().oldTime().boundaryField()
                        [
                            patchI
                        ]
                      - (Coo/Cn)*phi().oldTime().oldTime().boundaryField()
                        [
                            patchI
                        ];

                    rAUf_.boundaryField()[patchI] = deltaT/Cn;
                }
            }
        }
    }
}

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
    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
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
    ),
    pRefCell_(0),
    pRefValue_(0),
    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    sumLocalContErr_(0),
    globalContErr_(0),
    cumulativeContErr_(0)
{
    UisRequired();
    pisRequired();

    phi().oldTime();
    Uf_.oldTime();

    word ddtScheme
    (
        mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
    );

    if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        Uf_.oldTime().oldTime();
    }

    ddtU_.checkIn();

    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> pimpleFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() = rho_.value()
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

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


bool pimpleFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    pimpleControl& pimple = this->pimple();
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    surfaceScalarField& phi = this->phi();
    volScalarField& rAU = rAU_;
    const label pRefCell = pRefCell_;
    const scalar& pRefValue = pRefValue_;
    scalar& sumLocalContErr = sumLocalContErr_;
    scalar& globalContErr = globalContErr_;
    scalar& cumulativeContErr = cumulativeContErr_;

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
#       include "volContinuity.H"
    }

    if (correctPhi_ && meshChanged)
    {
        // Fluxes will be corrected to absolute velocity
        // HJ, 6/Feb/2009
#       include "correctPhi.foamextend.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    if (meshChanged)
    {
        // Update flux in FSI interface with Robin BC
        updateRobinFsiInterfaceFlux();
    }

    // CourantNo
    if (mesh.moving() && checkMeshCourantNo_)
    {
#       include "meshCourantNo.H"
    }

    if (meshChanged)
    {
#       include "CourantNo.H"
    }

    // --- PIMPLE loop
    while (pimple.loop())
    {
        // UEqn.H
        // Time-derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(U));

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(phi, U)
          + turbulence_->divDevReff()
        );

        if (pimple.momentumPredictor())
        {
            // Solve momentum predictor
            solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
        }

        // --- PISO loop
        while (pimple.correct())
        {
            p.boundaryField().updateCoeffs();

            // Prepare clean 1/a_p without time derivative and under-relaxation
            // contribution
            rAU = 1.0/HUEqn.A();
            rAUf_ = fvc::interpolate(rAU);

            // Calculate U from convection-diffusion matrix
            U = rAU*HUEqn.H();

            // Consistently calculate flux
            pimple.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

            // Update flux on the FSI interface for Robin BCs
            updateRobinFsiInterface();

            // Global flux balance
            adjustPhi(phi, U, p);

            // Non-orthogonal pressure corrector loop
            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        fvc::interpolate(rAU)/pimple.aCoeff(U.name()),
                        p,
                        "laplacian(rAU," + p.name() + ')'
                    )
                 ==
                    fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve
                (
                    mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
                );

                if (Pstream::parRun())
                {
                    // Needed for extrapolated BC
                    p.correctBoundaryConditions();
                }

                gradp() = fvc::grad(p);

                if (pimple.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

            // Correct flux due for Robin BCs according to temporal
            // discretization
            correctRobinFsiInterfaceFlux();

            // Explicitly relax pressure for momentum corrector except for last
            // corrector
            if (!pimple.finalIter())
            {
                p.relax();
            }

#           include "movingMeshContinuityErrs.H"

            // Consistently reconstruct velocity after pressure equation.
            // Note: flux is made relative inside the function
            pimple.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);

            gradU() = fvc::grad(U);
            ddtU_ = fvc::ddt(U);
        }

        turbulence_->correct();
    }

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi, U);

    // Update velocity on faces to account for modifications in the flux
    // when using the Robin BCs
    Uf_ = fvc::interpolate(U, "interpolate(U)");

    // Get tangential velocity
    Uf_ -= (mesh.Sf() & Uf_)*mesh.Sf()/magSqr(mesh.Sf());

    // Update with normal from flux
    Uf_ += phi*mesh.Sf()/magSqr(mesh.Sf());

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
