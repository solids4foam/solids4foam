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
addToRunTimeSelectionTable(physicsModel, icoFluid, fluid);
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
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> icoFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() = rho_.value()*nu_.value()*U().boundaryField()[patchID].snGrad();

    const vectorField n = mesh().boundary()[patchID].nf();

    tvF() -= n*(n & tvF());

    return tvF;
}


tmp<scalarField> icoFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> icoFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    tmp<scalarField> tMuEff
    (
        new scalarField
        (
            mesh().faceZones()[zoneID].size(),
            rho_.value()*nu_.value()
        )
    );

    return tMuEff;
}


bool icoFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    const fvMesh& mesh = this->mesh();

    const int nCorr = readInt(fluidProperties().lookup("nCorrectors"));

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    const int nOutCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pRefCell, pRefValue);

    for (int oCorr = 0; oCorr < nOutCorr; oCorr++)
    {
        if (mesh.moving())
        {
            // Make the fluxes relative
            phi() -= fvc::meshPhi(U());
        }

        // CourantNo
        {
          scalar CoNum = 0.0;
          scalar meanCoNum = 0.0;
          scalar velMag = 0.0;
          fluidModel::CourantNo(CoNum, meanCoNum, velMag);
        }

        fvVectorMatrix UEqn
        (
            fvm::ddt(U())
          + fvm::div(phi(), U())
          - fvm::laplacian(nu_, U())
        );

        solve(UEqn == -gradp());

        // --- PISO loop

        volScalarField rAU = 1.0/UEqn.A();
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        for (int corr = 0; corr < nCorr; corr++)
        {
            U() = rAU*UEqn.H();
            phi() = (fvc::interpolate(U()) & mesh.Sf());
            //+ fvc::ddtPhiCorr(rUA, U(), phi());

            fluidModel::updateRobinFsiInterface(p(), U(), phi(), rAUf);

            adjustPhi(phi(), U(), p());

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        rAUf, p(), "laplacian((1|A(U)),p)"
                    )
                 == fvc::div(phi())
                 // fvm::laplacian(rAUf, p()) == fvc::div(phi())
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                gradp() = fvc::grad(p());

                if (nonOrth == nNonOrthCorr)
                {
                    phi() -= pEqn.flux();
                }
            }

            fluidModel::continuityErrs();

            U() -= rAU*gradp();
            U().correctBoundaryConditions();

            gradU() = fvc::grad(U());
        }
    }

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
