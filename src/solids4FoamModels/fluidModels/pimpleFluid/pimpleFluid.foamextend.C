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
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pimpleFluid, 0);
addToRunTimeSelectionTable(fluidModel, pimpleFluid, dictionary);

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

#   include "volContinuity.H"

    if (correctPhi_ && meshChanged)
    {
        // Fluxes will be corrected to absolute velocity
        // HJ, 6/Feb/2009
#       include "correctPhi.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

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

        // Time derivative matrix
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

            // Calculate U from convection-diffusion matrix
            U = rAU*HUEqn.H();

            // Consistently calculate flux
            pimple.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

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

                pEqn.setReference(pRefCell_, pRefValue_);
                pEqn.solve
                (
                    mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
                );

                if (pimple.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

            // Explicitly relax pressure for momentum corrector except for last
            // corrector
            if (!pimple.finalIter())
            {
                p.relax();
            }

#           include "movingMeshContinuityErrs.H"

            // Consistently reconstruct velocity after pressure equation. Note: flux is
            // made relative inside the function
            pimple.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);
        }

        turbulence_->correct();
    }

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi, U);

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
