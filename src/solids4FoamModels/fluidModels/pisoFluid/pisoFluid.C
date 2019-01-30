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

#include "pisoFluid.H"
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

defineTypeNameAndDebug(pisoFluid, 0);
addToRunTimeSelectionTable(physicsModel, pisoFluid, fluid);
addToRunTimeSelectionTable(fluidModel, pisoFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pisoFluid::pisoFluid
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
    )
{
    UisRequired();
    pisRequired();

    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> pisoFluid::patchViscousForce(const label patchID) const
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

    // PC: why is this commented?
    //vectorField n = mesh().boundary()[patchID].nf();
    //tvF() -= n*(n & tvF());

    return tvF;
}


tmp<scalarField> pisoFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> pisoFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
        rho_.value()*turbulence_->nuEff()().boundaryField()[patchID];

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


bool pisoFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    const fvMesh& mesh = fluidModel::mesh();

    const int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

    const int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pRefCell, pRefValue);

    // if (mesh.moving())
    // {
    //     // Make the fluxes relative
    //     phi() -= fvc::meshPhi(U());
    // }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

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

    solve(UEqn == -gradp());

    // --- PISO loop

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rUA));

    for (int corr=0; corr < nCorr; corr++)
    {
        U() = rUA*UEqn.H();
        phi() = (fvc::interpolate(U()) & mesh.Sf());

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf, p()) == fvc::div(phi())
            );

            pEqn.setReference(pRefCell, pRefValue);

            if
            (
                corr == nCorr-1
             && nonOrth == nNonOrthCorr
            )
            {
                pEqn.solve(mesh.solutionDict().solver("pFinal"));
            }
            else
            {
                pEqn.solve();
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi() -= pEqn.flux();
            }
        }

        // Continuity error
        fluidModel::continuityErrs();

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi(), U());

        gradp() = fvc::grad(p());

        U() -= rUA*gradp();
        U().correctBoundaryConditions();

        gradU() = fvc::grad(U());
    }

    turbulence_->correct();

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
