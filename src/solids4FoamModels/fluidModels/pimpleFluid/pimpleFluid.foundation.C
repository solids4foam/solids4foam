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
#include "CorrectPhi.H"
#include "fv.H"
#include "fvc.H"
#include "fvm.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "volFields.H"
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

void pimpleFluid::updateRobinFsiInterface(surfaceScalarField& phiHbyA)
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
            const word ddtScheme(mesh().ddtScheme("ddt(" + U().name() +')'));

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
                if (runTime().timeIndex() == 1)
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
    LTS_(fv::localEulerDdt::enabled(mesh())),
    trDeltaT_(),
    Uf_(),
    ddtU_(fvc::ddt(U())),
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
        calculatedFvPatchScalarField::typeName
    ),
    pressureReference_(p(), pimple().dict()),
    laminarTransport_(U(), phi()),
    turbulence_
    (
        incompressible::momentumTransportModel::New
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
    cumulativeContErr_(0)
{
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

    Uf_().oldTime();

    word ddtScheme
    (
        mesh().ddtScheme("ddt(" + U().name() +')')
    );

    if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        Uf_().oldTime().oldTime();
    }

    
    if (LTS_)
    {
        Info<< "Using LTS" << endl;

        trDeltaT_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
    else
    {
        const fvMesh& mesh = this->mesh();
        const surfaceScalarField& phi = this->phi();
        #include "CourantNo.H"
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
          & (-turbulence_->devTau()().boundaryField()[patchID])
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
    pressureReference& pressureReference = pressureReference_;
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
        if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
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
                mesh.update();
            }

            if (mesh.changing())
            {
                // MRF not added yet
                //MRF.update();

                if (correctPhi)
                {
                    #include "correctPhi.foundation.H"
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }

        // fvModels not implemented yet
        //fvModels.correct();

        // UEqn.H

        // Solve the Momentum equation

        // MRF not implemented yet
        // MRF.correctBoundaryVelocity(U);

        tmp<fvVectorMatrix> tUEqn
        (
            fvm::ddt(U) + fvm::div(phi, U)
            // + MRF.DDt(U)
          + turbulence_->divDevSigma(U)
         // ==
         //    fvModels.source(U)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        UEqn.relax();

        // fvConstraints not implemented yet
        //fvConstraints.constrain(UEqn);

        if (pimple.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));

            // fvConstraints not implemented yet
            //fvConstraints.constrain(U);
        }

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            rAU_ = 1.0/UEqn.A();
            // volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU_*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU_)*fvc::ddtCorr(U, phi, Uf)
            //+ MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf))
            );

            //MRF.makeRelative(phiHbyA);

            // Update flux on the FSI interface for Robin BCs
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
                    fvc::interpolate(rAtU() - rAU_)*fvc::snGrad(p)*mesh.magSf();
                HbyA -= (rAU_ - rAtU())*fvc::grad(p);
            }

            if (pimple.nCorrPiso() <= 1)
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

                pEqn.setReference
                (
                    pressureReference.refCell(),
                    pressureReference.refValue()
                );

                pEqn.solve();

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
            // fvConstraints.constrain(U);

            gradU() = fvc::grad(U);
            ddtU_ = fvc::ddt(U);
            
            // Correct Uf if the mesh is moving
            fvc::correctUf(Uf, U, phi);
            
            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);
        }


        if (pimple.turbCorr())
        {
            laminarTransport_.correct();
            turbulence_->correct();
        }
    }

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi, U);


    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
