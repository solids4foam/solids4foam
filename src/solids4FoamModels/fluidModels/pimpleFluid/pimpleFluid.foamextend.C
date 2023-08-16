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
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "wedgeFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "extrapolatedFvPatchFields.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "flowRateOutletPressureFvPatchScalarField.H"
#include "thermalRobinFvPatchScalarField.H"

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
         // || (
         //        isA<elasticWallPressureFvPatchScalarField>
         //        (
         //            p().boundaryField()[patchI]
         //        )
         //     && isA<robinElasticWallVelocityFvPatchVectorField>
         //        (
         //            U().boundaryField()[patchI]
         //        )
         //    )
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
         // || (
         //        isA<elasticWallPressureFvPatchScalarField>
         //        (
         //            p().boundaryField()[patchI]
         //        )
         //     && isA<robinElasticWallVelocityFvPatchVectorField>
         //        (
         //            U().boundaryField()[patchI]
         //        )
         //    )
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
                (
                    Uf_.oldTime().boundaryField()[patchI] &
                    mesh().Sf().boundaryField()[patchI]
                );
                // phi().oldTime().boundaryField()[patchI];

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
                    (
                        Uf_.oldTime().boundaryField()[patchI] &
                        mesh().Sf().boundaryField()[patchI]
                    );

                    rAUf_.boundaryField()[patchI] =
                        runTime().deltaT().value();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

                    phi().boundaryField()[patchI] =
                        (Co/Cn)*
                        (
                            Uf_.oldTime().boundaryField()[patchI] &
                            mesh().Sf().boundaryField()[patchI]
                        )
                      - (Coo/Cn)*
                        (
                            Uf_.oldTime().oldTime().boundaryField()[patchI] &
                            mesh().Sf().boundaryField()[patchI]
                        );

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
         // || (
         //        isA<elasticWallPressureFvPatchScalarField>
         //        (
         //            p().boundaryField()[patchI]
         //        )
         //     && isA<robinElasticWallVelocityFvPatchVectorField>
         //        (
         //            U().boundaryField()[patchI]
         //        )
         //    )
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
    ddtU_
    (
        IOobject
        (
            "ddt(U)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimVelocity/dimTime, vector::zero)
    ),
    laplNuU_
    (
        IOobject
        (
            "laplacian(nuEff,U)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimVelocity/dimTime, vector::zero),
        extrapolatedFvPatchVectorField::typeName
    ),
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
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        calculatedFvPatchScalarField::typeName
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
    refPressureSeries_(),
    ddtCorr_
    (
        fluidProperties().lookupOrDefault<Switch>
        (
            "ddtCorr",
            false
        )
    ),
    adjustPhi_(true),
    solveEnergyEq_
    (
        fluidProperties().lookupOrDefault<Switch>
        (
            "solveEnergyEq",
            false
        )
    ),
    TPtr_(),
    lambdaEffPtr_()
{
    UisRequired();
    pisRequired();

    phi().oldTime();
    Uf_.oldTime();

    word ddtScheme
    (
        mesh().schemesDict().ddtScheme("ddt(" + U().name() +')')
    );

    // Check if flux should be adjusted
    forAll(p().boundaryField(), patchI)
    {
        if
        (
            isA<flowRateOutletPressureFvPatchScalarField>
            (
                p().boundaryField()[patchI]
            )
        )
        {
            adjustPhi_ = false;
        }
    }

    if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        Uf_.oldTime().oldTime();
    }

    ddtU_.checkIn();

#ifdef OPENFOAM_NOT_EXTEND
    mesh().setFluxRequired(p().name());
#else
    mesh().schemesDict().setFluxRequired(p().name());
#endif

    // Read reference pressure series if exist
    if (pimple().dict().found("refPressureSeries"))
    {
        // if (debug)
        {
            Info<< "    reference pressure is time-varying" << endl;
        }

        refPressureSeries_ =
            interpolationTable<scalar>
            (
                pimple().dict().subDict("refPressureSeries")
            );
    }

    // Create temperature field if necessary
    if (solveEnergyEq_)
    {
        TPtr_.set
        (
            new volScalarField
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
            )
        );

        // Heat capacity [J/kgK]
        dimensionedScalar Cp(laminarTransport_.lookup("Cp"));

        // Thermal conductivity [W/(m K)]
        dimensionedScalar lambda(laminarTransport_.lookup("lambda"));

        // Turbulent Prandtl number
        dimensionedScalar Prt(laminarTransport_.lookup("Prt"));

        lambdaEffPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "lambdaEff",
                    runTime.timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                lambda + rho_*Cp*(turbulence_->nut()/Prt)
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> pimpleFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tvF.ref() =
#else
    tvF() =
#endif
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );
        // rho_.value()*nu_.value()*U().boundaryField()[patchID].snGrad();

    return tvF;
}


tmp<scalarField> pimpleFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tpF.ref() =
#else
    tpF() =
#endif
        rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


bool pimpleFluid::evolve()
{
    Info << "Evolving fluid model: " << this->type() << endl;

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

    phi().oldTime().oldTime();
    Uf_.oldTime().oldTime();

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // CourantNo
    fluidModel::CourantNo();

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), pimple().dict(), pRefCell, pRefValue);

    if (refPressureSeries_.size())
    {
        pRefValue = refPressureSeries_(runTime().timeOutputValue());
    }

    while (pimple().loop())
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(U())
          + fvm::div(phi(), U())
          + turbulence_->divDevReff()
          - boussinesqMomentumSource()
        );

        UEqn.relax();

        if (pimple().momentumPredictor())
        {
            solve(UEqn == -gradp());
        }

        // --- PISO loop
        //
        volScalarField AU("AU", UEqn.A());

        while (pimple().correct())
        {
            rAU_ = 1.0/UEqn.A();
            rAUf_ = fvc::interpolate(rAU_);

            U() = rAU_*UEqn.H();
            phi() = (mesh.Sf() & fvc::interpolate(U()));

            // Correct phi for temporal consistency
            if (ddtCorr_)
            {
                phi() += fvc::ddtConsistentPhiCorr(Uf_, U(), rAUf_);
            }

            updateRobinFsiInterface();

            if (adjustPhi_)
            {
                adjustPhi(phi(), U(), p());
            }

            // Non-orthogonal pressure corrector loop
            while (pimple().correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        rAUf_, p(), "laplacian((1|A(U)),p)"
                    )
                 == fvc::div(phi())
                );

                pEqn.setReference(pRefCell, pRefValue);

#ifdef OPENFOAM_NOT_EXTEND
                pEqn.solve();
#else
                pEqn.solve
                (
                    mesh.solutionDict()
                   .solver(p().select(pimple().finalInnerIter()))
                );
#endif

                if (Pstream::parRun())
                {
                    // Needed for extrapolated BC
                    p().correctBoundaryConditions();
                }

                if (pimple().finalNonOrthogonalIter())
                {
                    phi() -= pEqn.flux();
                }
            }

            fluidModel::continuityErrs();

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi(), U());

            p().relax();

            gradp() = fvc::grad(p());

            U() -= rAU_*gradp();
            U().correctBoundaryConditions();

            gradU() = fvc::grad(U());
            ddtU_ = fvc::ddt(U());
        }

        // Needed for wallPressure bc
        laplNuU_ = fvc::laplacian(turbulence_->nuEff(), U());
        laplNuU_.correctBoundaryConditions();

        // Solve energy equation if required
        solveEnergyEq();
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    // Update velocity on faces
    Uf_ = fvc::interpolate(U(), "interpolate(U)");

    // Get tangential velocity
    Uf_ -= (mesh.Sf() & Uf_)*mesh.Sf()/magSqr(mesh.Sf());

    // Update with normal from flux
    Uf_ += phi()*mesh.Sf()/magSqr(mesh.Sf());

    return true;
}


void pimpleFluid::writeFields(const Time& runTime)
{
    fluidModel::writeFields(runTime);
}


void pimpleFluid::solveEnergyEq()
{
    if (solveEnergyEq_)
    {
        // Heat capacity [J/kgK]
        dimensionedScalar Cp(laminarTransport_.lookup("Cp"));

        // Thermal conductivity [W/(m K)]
        dimensionedScalar lambda(laminarTransport_.lookup("lambda"));

        // Turbulent Prandtl number
        dimensionedScalar Prt(laminarTransport_.lookup("Prt"));

        // Update effective thermal conductivity
        lambdaEffPtr_() =
            lambda + rho_*Cp*(turbulence_->nut()/Prt);

        // Store fields for under-relaxation and residual calculation
        TPtr_().storePrevIter();

        fvScalarMatrix TEqn
        (
            rho_*Cp*
            (
                fvm::ddt(TPtr_())
              + fvm::div(phi(), TPtr_())
            )
          - fvm::laplacian(lambdaEffPtr_(), TPtr_())
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        TEqn.solve();

        // Under-relax the field
        TPtr_().relax();
    }
}


tmp<volVectorField> pimpleFluid::boussinesqMomentumSource() const
{
    tmp<volVectorField> tSource
    (
        new volVectorField
        (
            IOobject
            (
                "boussinesqMomentumSource",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector("0", dimVelocity/dimTime, vector::zero)
        )
    );

    if (solveEnergyEq_)
    {
        // Thermal expansion coefficient [1/K]
        dimensionedScalar beta(laminarTransport_.lookup("beta"));

        // Reference temperature [K]
        dimensionedScalar TRef(laminarTransport_.lookup("TRef"));

        tSource() = -beta*(TPtr_() - TRef)*g();
    }

    return tSource;
}


void pimpleFluid::setTemperatureAndHeatFlux
(
    fvPatchScalarField& temperaturePatch,
    const scalarField& temperature,
    const scalarField& heatFlux
)
{
    if (temperaturePatch.type() == thermalRobinFvPatchScalarField::typeName)
    {
        thermalRobinFvPatchScalarField& patchT =
            refCast<thermalRobinFvPatchScalarField>(temperaturePatch);

        patchT.temperature() = temperature;
        patchT.heatFlux() = heatFlux;
    }
    else
    {
        FatalErrorIn
        (
            "void pimpleFluid::setTemperature\n"
            "(\n"
            "    fvPatchScalarField& temperaturePatch,\n"
            "    const scalarField& temperature,\n"
            "    const scalarField& heatFlux\n"
            ")"
        )   << "Boundary condition "
            << temperaturePatch.type()
            << " for patch " << temperaturePatch.patch().name()
            << " should instead be type "
            << thermalRobinFvPatchScalarField::typeName
            << abort(FatalError);
    }
}


void pimpleFluid::setTemperatureAndHeatFlux
(
    const label interfaceI,
    const label patchID,
    const scalarField& faceZoneTemperature,
    const scalarField& faceZoneHeatFlux
)
{
    const scalarField patchTemperature
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTemperature)
    );

    const scalarField patchHeatFlux
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneHeatFlux)
    );

#ifdef OPENFOAM_NOT_EXTEND
    setTemperatureAndHeatFlux
    (
        TPtr_().boundaryFieldRef()[patchID],
        patchTemperature,
        patchHeatFlux
    );
#else
    setTemperatureAndHeatFlux
    (
        TPtr_().boundaryField()[patchID],
        patchTemperature,
        patchHeatFlux
    );
#endif
}


void pimpleFluid::setEqInterHeatTransferCoeff
(
    fvPatchScalarField& temperaturePatch,
    const scalarField& heatTransferCoeff
)
{
    if (temperaturePatch.type() == thermalRobinFvPatchScalarField::typeName)
    {
        thermalRobinFvPatchScalarField& patchT =
            refCast<thermalRobinFvPatchScalarField>(temperaturePatch);

        patchT.eqInterHeatTransferCoeff() = heatTransferCoeff;
    }
    else
    {
        FatalErrorIn
        (
            "void pimpleFluid::setEqInterHeatTransferCoeff\n"
            "(\n"
            "    fvPatchScalarField& temperaturePatch,\n"
            "    const scalarField& heatTransferCoeff\n"
            ")"
        )   << "Boundary condition "
            << temperaturePatch.type()
            << " for patch " << temperaturePatch.patch().name()
            << " should instead be type "
            << thermalRobinFvPatchScalarField::typeName
            << abort(FatalError);
    }
}


void pimpleFluid::setEqInterHeatTransferCoeff
(
    const label interfaceI,
    const label patchID,
    const scalarField& faceZoneHTC
)
{
    const scalarField patchHTC
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneHTC)
    );

#ifdef OPENFOAM_NOT_EXTEND
    setEqInterHeatTransferCoeff
    (
        TPtr_().boundaryFieldRef()[patchID],
        patchHTC
    );
#else
    setEqInterHeatTransferCoeff
    (
        TPtr_().boundaryField()[patchID],
        patchHTC
    );
#endif
}


tmp<scalarField> pimpleFluid::patchTemperature(const label patchID) const
{
    tmp<scalarField> tT
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tT.ref() =
#else
    tT() =
#endif
        TPtr_().boundaryField()[patchID];

    return tT;
}

tmp<scalarField> pimpleFluid::patchHeatFlux(const label patchID) const
{
    tmp<scalarField> tHF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tHF.ref() =
#else
    tHF() =
#endif
        lambdaEffPtr_().boundaryField()[patchID]*
        TPtr_().boundaryField()[patchID].snGrad();

    return tHF;
}


tmp<scalarField> pimpleFluid::patchHeatTransferCoeff
(
    const label patchID
) const
{
    tmp<scalarField> tHTC
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tHTC.ref() =
#else
    tHTC() =
#endif
        (1.0/mesh().deltaCoeffs().boundaryField()[patchID])/
        lambdaEffPtr_().boundaryField()[patchID];

    return tHTC;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
