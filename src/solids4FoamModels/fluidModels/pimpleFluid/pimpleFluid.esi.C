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
#include "CorrectPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "findRefCell.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
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

void pimpleFluid::updateRobinFsiInterface
(
    surfaceScalarField& phiHbyA
)
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
                mesh().ddtScheme("ddt(" + U().name() +')')
            );

            if
            (
                ddtScheme
             == fv::EulerDdtScheme<vector>::typeName
            )
            {
                phiHbyA.boundaryFieldRef()[patchI] =
                (
                    Uf_().oldTime().boundaryField()[patchI] &
                    mesh().Sf().boundaryField()[patchI]
                );

                rAU_.boundaryFieldRef()[patchI] =
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
                    phiHbyA.boundaryFieldRef()[patchI] =
                    (
                        Uf_().oldTime().boundaryField()[patchI] &
                        mesh().Sf().boundaryField()[patchI]
                    );

                    rAU_.boundaryFieldRef()[patchI] =
                        runTime().deltaT().value();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

                    phiHbyA.boundaryFieldRef()[patchI] =
                        (Co/Cn)*
                        (
                            Uf_().oldTime().boundaryField()[patchI] &
                            mesh().Sf().boundaryField()[patchI]
                        )
                      - (Coo/Cn)*
                        (
                            Uf_().oldTime().oldTime().boundaryField()[patchI] &
                            mesh().Sf().boundaryField()[patchI]
                        );

                    rAU_.boundaryFieldRef()[patchI] = deltaT/Cn;
                }
            }
        }
        else if
        (
            isA<movingWallPressureFvPatchScalarField>
            (
                p().boundaryField()[patchI]
            )
        )
        {
            phiHbyA.boundaryFieldRef()[patchI] +=
                rAU_.boundaryField()[patchI]
               *p().boundaryField()[patchI].snGrad()
               *mesh().magSf().boundaryField()[patchI];
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
    Uf_(),
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
    pRefCell_(0),
    pRefValue_(0),
    laminarTransport_(U(), phi()),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
    rho_(laminarTransport_.lookup("rho")),
    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    moveMeshOuterCorrectors_
    (
        pimple().dict().lookupOrDefault("moveMeshOuterCorrectors", false)
    ),
    cumulativeContErr_(0),
    solveEnergyEq_
    (
        fluidProperties().lookupOrAddDefault<Switch>
        (
            "solveEnergyEq",
            false
        )
    ),
    TPtr_(),
    lambdaEffPtr_()
{
    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().setFluxRequired(p().name());
    turbulence_->validate();

    if (mesh().dynamic())
    {
        Info<< "Constructing face velocity Uf\n" << endl;

        Uf_.reset
        (
            new surfaceVectorField
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
            )
        );

        Uf_().oldTime();

        if
        (
            word(mesh().ddtScheme("ddt(" + U().name() +')'))
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            Uf_().oldTime().oldTime();
        }
    }

    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& phi = this->phi();
    #include "CourantNo.H"

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
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
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
                    mesh,
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

    tvF.ref() = rho_.value()
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

    tpF.ref() = rho_.value()*p().boundaryField()[patchID];

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
    autoPtr<surfaceVectorField>& Uf = Uf_;
    scalar& cumulativeContErr = cumulativeContErr_;
    const bool correctPhi = correctPhi_;
    const bool checkMeshCourantNo = checkMeshCourantNo_;
    const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            // fvModels not added yet
            // fvModels.preUpdateMesh();

            // Ideally we would not need a specific FSI mesh update function
            // Hopefully we can remove the need for it soon
            if (fluidModel::fsiMeshUpdate())
            {
                // The FSI interface is in charge of calling mesh.update()
                fluidModel::fsiMeshUpdateChanged();
            }
            else
            {
                // Do any mesh changes
                mesh.controlledUpdate();
            }

            if (mesh.changing())
            {
                // MRF not added yet
                // MRF.update();

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh.Sf() & Uf();

                    #include "correctPhi.esi.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }

        // Solve the Momentum equation

        // MRF not implemented yet
        // MRF.correctBoundaryVelocity(U);

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + turbulence_->divDevReff(U)
          - boussinesqMomentumSource()
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();

        // fvOptions not implemented yet
        // fvOptions.constrain(UEqn);

        if (pimple.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
            // fvOptions not implemented yet
            // fvOptions.correct(U);
        }

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            // pEqn.H

            rAU_ = 1.0/UEqn.A();
            volVectorField HbyA(constrainHbyA(rAU_*UEqn.H(), U, p));
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

            // MRF not yet implemented
            // if (pimple.ddtCorr())
            // {
            //     phiHbyA +=
            //     MRF.zeroFilter
            //     (
            //         fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf)
            //     );
            // }
            // else
            // {
            //     phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
            // }
            // MRF.makeRelative(phiHbyA);

            if (pimple.ddtCorr())
            {
                phiHbyA +=
                    fvc::interpolate(rAU_)*fvc::ddtCorr(U, phi, Uf);
            }

            updateRobinFsiInterface(phiHbyA);

            if (p.needReference())
            {
                fvc::makeRelative(phiHbyA, U);
                adjustPhi(phiHbyA, U, p);
                fvc::makeAbsolute(phiHbyA, U);
            }

            tmp<volScalarField> rAtU(rAU_);

            if (pimple.consistent())
            {
                rAtU = 1.0/max(1.0/rAU_ - UEqn.H1(), 0.1/rAU_);
                phiHbyA +=
                    fvc::interpolate(rAtU() - rAU_)*
                    fvc::snGrad(p)*mesh.magSf();
                HbyA -= (rAU_ - rAtU())*fvc::grad(p);
            }

            if (pimple.nCorrPISO() <= 1)
            {
                tUEqn.clear();
            }

            // Update the pressure BCs to ensure flux consistency
            // constrainPressure(p, U, phiHbyA, rAtU(), MRF);
            constrainPressure(p, U, phiHbyA, rAtU());

            // Non-orthogonal pressure corrector loop
            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell_, pRefValue_);

                pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

                if (pimple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            U = HbyA - rAtU*fvc::grad(p);
            U.correctBoundaryConditions();
         // fvOptions.correct(U);
        }

        // Correct Uf if the mesh is moving
        fvc::correctUf(Uf, U, phi);

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (pimple.turbCorr())
        {
            laminarTransport_.correct();
            turbulence_->correct();
        }

        // Solve energy equation if required
        solveEnergyEq();
    }

    return 0;
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

        tSource.ref() = -beta*(TPtr_() - TRef)*g();
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
