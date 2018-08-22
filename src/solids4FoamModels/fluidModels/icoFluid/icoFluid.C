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
#include "pisoControl.H"

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
    rho_(transportProperties_.lookup("rho")),
    pRefCell_(0),
    pRefValue_(0)
{
    setRefCell(p(), piso().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p().name());
}

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

    fvMesh& mesh = this->mesh();

    bool meshChanged = false;
    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = refCast<dynamicFvMesh>(mesh).update();
        reduce(meshChanged, orOp<bool>());
    }

    if (meshChanged)
    {
        const Time& runTime = fluidModel::runTime();
#       include "volContinuity.H"
    }
    
    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        fluidModel::CourantNo(CoNum, meanCoNum, velMag);
    }

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

    // Prepare clean 1/a_p without time derivative contribution
    volScalarField rAU = 1.0/HUEqn.A();

    while (piso().correct())
    {
        // Calculate U from convection-diffusion matrix
        U() = rAU*HUEqn.H();

        // Consistently calculate flux
        piso().calcTransientConsistentFlux(phi(), U(), rAU, ddtUEqn);

        adjustPhi(phi(), U(), p());

        // Non-orthogonal pressure corrector loop
        while (piso().correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    fvc::interpolate(rAU)/piso().aCoeff(U().name()),
                    p(),
                    "laplacian(rAU," + p().name() + ')'
                )
             == fvc::div(phi())
            );

            pEqn.setReference(pRefCell_, pRefValue_);
            pEqn.solve
            (
                mesh.solutionDict().solver(p().select(piso().finalInnerIter()))
            );

            if (piso().finalNonOrthogonalIter())
            {
                phi() -= pEqn.flux();
            }
        }

        fluidModel::continuityErrs();

        // Consistently reconstruct velocity after pressure equation.
        // Note: flux is made relative inside the function
        piso().reconstructTransientVelocity(U(), phi(), ddtUEqn, rAU, p());
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());
    
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
