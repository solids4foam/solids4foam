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

#include "transientSimpleFluid.H"
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

defineTypeNameAndDebug(transientSimpleFluid, 0);
addToRunTimeSelectionTable(physicsModel, transientSimpleFluid, fluid);
addToRunTimeSelectionTable(fluidModel, transientSimpleFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transientSimpleFluid::transientSimpleFluid
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
    rho_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("rho")
    ),
    consistencyByJasak_
    (
        fluidProperties().lookupOrDefault<Switch>("consistencyByJasak", false)
    )
{
    UisRequired();
    pisRequired();

    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> transientSimpleFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & turbulence_->devReff()().boundaryField()[patchID]
        );

    // PC: why is this commented
    //vectorField n = mesh().boundary()[patchID].nf();
    //tvF() -= n*(n & tvF());

    return tvF;
}


tmp<scalarField> transientSimpleFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


bool transientSimpleFluid::evolve()
{
    if (consistencyByJasak_)
    {
        return evolveConsistentByJasak();
    }
    else
    {
        return evolveInconsistent();
    }
}


bool transientSimpleFluid::evolveInconsistent()
{
    Info<< "Evolving fluid model" << endl;

    const fvMesh& mesh = fluidModel::mesh();

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    const int nOuterCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));

    scalar convergenceCriterion = 0;
    fluidProperties().readIfPresent("convergence", convergenceCriterion);

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pRefCell, pRefValue);

    phi().oldTime();

    for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
    {
        scalar eqnResidual = 1, maxResidual = 0;
        p().storePrevIter();

        if (mesh.moving())
        {
            // Make the fluxes relative
            phi() -= fvc::meshPhi(U());
        }

        // Calculate CourantNo
        {
            scalar CoNum = 0.0;
            scalar meanCoNum = 0.0;
            scalar velMag = 0.0;
            fluidModel::CourantNo(CoNum, meanCoNum, velMag);
        }

        // Construct momentum equation
        fvVectorMatrix UEqn
        (
            fvm::ddt(U())
          + fvm::div(phi(), U())
          + turbulence_->divDevReff()
        );

        UEqn.relax();

        // Solve momentum equation
        eqnResidual =
            solve(UEqn == -gradp()).initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

        volScalarField aU = UEqn.A();

        U() = UEqn.H()/aU;
        phi() = (fvc::interpolate(U()) & mesh.Sf());

#       include "adjustPhi.H"

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            // Construct pressure equation
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1/aU, p()) == fvc::div(phi())
            );

            // Solve pressure equation

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi() -= pEqn.flux();
            }
        }

        // Calculate Continuity error
        fluidModel::continuityErrs();

        // Explicitly relax pressure for momentum corrector
        p().relax();

        gradp() = fvc::grad(p());

        U() -= gradp()/aU;
        U().correctBoundaryConditions();

        turbulence_->correct();

        if (maxResidual < convergenceCriterion)
        {
            Info<< "reached convergence criterion: "
                << convergenceCriterion << endl;
            Info<< "Number of iterations: " << oCorr << endl;
            break;
        }
    }

    return 0;
}


bool transientSimpleFluid::evolveConsistentByJasak()
{
    Info<< "Evolving fluid model (consistency by Jasak)" << endl;

    const fvMesh& mesh = fluidModel::mesh();

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    const int nOuterCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));

    scalar convergenceCriterion = 0;
    fluidProperties().readIfPresent("convergence", convergenceCriterion);

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pRefCell, pRefValue);

    phi().oldTime();

    for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
    {
        scalar eqnResidual = 1, maxResidual = 0;
        p().storePrevIter();

        if (mesh.moving())
        {
            // Make the fluxes relative
            phi() -= fvc::meshPhi(U());
        }

        // Calculate CourantNo
        {
            scalar CoNum = 0.0;
            scalar meanCoNum = 0.0;
            scalar velMag = 0.0;
            fluidModel::CourantNo(CoNum, meanCoNum, velMag);
        }

        // Construct momentum equation
        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(phi(), U())
          + turbulence_->divDevReff()
        );

        // Time derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(U()));

        tmp<fvVectorMatrix> UEqn(ddtUEqn + HUEqn);

        UEqn().relax();

        // Solve momentum equation
        eqnResidual =
            solve(UEqn() == -gradp()).initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

        UEqn.clear();

        volScalarField aU = HUEqn.A();

        U() = HUEqn.H()/aU;
        phi() = (fvc::interpolate(U()) & mesh.Sf());

#       include "adjustPhi.H"

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            // Construct pressure equation
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1/aU, p()) == fvc::div(phi())
            );

            // Solve pressure equation

            pEqn.setReference(pRefCell, pRefValue);

            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi() -= pEqn.flux();
            }
        }

        // Calculate Continuity error
        fluidModel::continuityErrs();

        // Explicitly relax pressure for momentum corrector
        p().relax();

        gradp() = fvc::grad(p());

        U() = 1.0/(aU + ddtUEqn.A())*
            (
                U()*aU - gradp() + ddtUEqn.H()
            );
        U().correctBoundaryConditions();

        gradU() = fvc::grad(U());

        turbulence_->correct();

        if (maxResidual < convergenceCriterion)
        {
            Info<< "reached convergence criterion: "
                << convergenceCriterion << endl;
            Info<< "Number of iterations: " << oCorr << endl;
            break;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
