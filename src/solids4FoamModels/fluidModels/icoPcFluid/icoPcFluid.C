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

#include "icoPcFluid.H"
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

defineTypeNameAndDebug(icoPcFluid, 0);
addToRunTimeSelectionTable(physicsModel, icoPcFluid, fluid);
addToRunTimeSelectionTable(fluidModel, icoPcFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

icoPcFluid::icoPcFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    laminarTransport_(U(), phi()),
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
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> icoPcFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()*laminarTransport_.nu().boundaryField()[patchID]
       *U().boundaryField()[patchID].snGrad();

    const vectorField n = mesh().boundary()[patchID].nf();

    tvF() -= n*(n & tvF());

    return tvF;
}


tmp<scalarField> icoPcFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> icoPcFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
        rho_.value()*laminarTransport_.nu().boundaryField()[patchID];

    tmp<scalarField> tMuEff
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& muEff = tMuEff();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pMuEff, I)
    {
        muEff[mesh().faceZones()[zoneID].whichFace(patchStart + I)] =
            pMuEff[I];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(muEff, sumOp<scalarField>());

    return tMuEff;
}


bool icoPcFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    const fvMesh& mesh = fluidModel::mesh();

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    const int nOuterCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));

    const scalar convergenceCriterion =
        fluidProperties().lookupOrDefault<int>("convergence", 0);

    const Switch pressureCorr
    (
        fluidProperties().lookup("pressureCorr")
    );

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pRefCell, pRefValue);

    for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
    {
        p().storePrevIter();

        scalar eqnResidual = 1;
        scalar maxResidual = 0;
    
        // Calculate CourantNo
        {
            scalar CoNum = 0.0;
            scalar meanCoNum = 0.0;
            scalar velMag = 0.0;
            fluidModel::CourantNo(CoNum, meanCoNum, velMag);
        }

        fvVectorMatrix UEqn
        (
            fvm::div(phi(), U())
          - fvm::laplacian(laminarTransport_.nu(), U())
        );

        volScalarField rUA = 1.0/UEqn.A();

        UEqn += fvm::ddt(U());
        
        UEqn.relax();

        eqnResidual = solve(UEqn == -gradp()).initialResidual();

        maxResidual = max(eqnResidual, maxResidual);

        if (pressureCorr)
        {
            surfaceVectorField avgU = fvc::interpolate(U());
            surfaceVectorField avgGradp = fvc::interpolate(gradp());

            forAll(p().boundaryField(), patchI)
            {
                if (p().boundaryField()[patchI].fixesValue())
                {
                    avgU.boundaryField()[patchI] = 
                        U().boundaryField()[patchI].patchInternalField();
                    avgGradp.boundaryField()[patchI] = 
                        gradp().boundaryField()[patchI].patchInternalField();
                }
            }

            phi() =
                (avgU & mesh.Sf()) 
              + fvc::interpolate(rUA)*(mesh.Sf() & avgGradp);

            U() += rUA*gradp();
        }
        else
        {
            U() = rUA*UEqn.H();
            phi() = (fvc::interpolate(U()) & mesh.Sf());

            forAll(p().boundaryField(), patchI)
            {
                if 
                (
                    !p().boundaryField()[patchI].fixesValue()
                 && isA<zeroGradientFvPatchVectorField>
                    (
                        U().boundaryField()[patchI]
                    )
                )
                {
                    // If zeroGradient bc for pressure and velocity
                    // at outlet
                    phi().boundaryField()[patchI] =
                    (
                        U().boundaryField()[patchI].patchInternalField()
                      & mesh.Sf().boundaryField()[patchI]
                    );
                }
            }
        }

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rUA, p())
             == fvc::div(phi())
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi() -= pEqn.flux();
            }
        }

        // Continuity error
        fluidModel::continuityErrs();

        p().relax();
        gradp() = fvc::grad(p());

        // Momentum corrector
        U() -= rUA*gradp();
        U().correctBoundaryConditions();

        gradU() = fvc::grad(U());

        if (maxResidual < convergenceCriterion)
        {
            Info<< "Reached convergence criterion: " 
                << convergenceCriterion << endl
                << "Number of iterations: " << oCorr << endl;
            break;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
