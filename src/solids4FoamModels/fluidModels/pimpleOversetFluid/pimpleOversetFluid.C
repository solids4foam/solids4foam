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

#if FOAMEXTEND

#include "pimpleOversetFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"

#include "oversetFvPatchFields.H"
#include "oversetAdjustPhi.H"
#include "globalOversetAdjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pimpleOversetFluid, 0);
addToRunTimeSelectionTable(fluidModel, pimpleOversetFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pimpleOversetFluid::pimpleOversetFluid
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
    pRefCell_(0),
    pRefValue_(0)
{
    osMesh();
    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> pimpleOversetFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> pimpleOversetFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> pimpleOversetFluid::faceZoneMuEff
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


bool pimpleOversetFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    fvMesh& mesh = fluidModel::mesh();

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
    fluidModel::oversetCourantNo();

    // --- PIMPLE loop
    while (pimple().loop())
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(U())
          + fvm::div(phi(), U())
          + turbulence_->divDevReff()
        );

        UEqn.relax
        (
            mesh.solutionDict().equationRelaxationFactor
            (
                U().select(pimple().finalIter())
            )
        );

        if (pimple().momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p()));
        }

        // --- PISO loop

        while (pimple().correct())
        {
            p().boundaryField().updateCoeffs();
            
            volScalarField rAU = 1.0/UEqn.A();
            oversetFvPatchScalarField::oversetInterpolate(rAU);
            surfaceScalarField rAUf = fvc::interpolate(rAU);

            U() = rAU*UEqn.H();
            oversetFvPatchVectorField::oversetInterpolate(U());

            phi() = fvc::interpolate(U()) & mesh.Sf();
            
            // Adjust overset fluxes
            oversetAdjustPhi(phi(), U()); // Fringe flux adjustment
            globalOversetAdjustPhi(phi(), U(), p()); // Global flux adjustment

            // Non-orthogonal pressure corrector loop
            while (pimple().correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAUf, p())
                 == fvc::div(phi())
                );

                // Adjust non-orthogonal fringe fluxes if necessary
                osMesh().correctNonOrthoFluxes(pEqn, U());

                pEqn.setReference(pRefCell_, pRefValue_);
                pEqn.solve
                (
                    mesh.solutionDict().solver
                    (
                        p().select(pimple().finalInnerIter())
                    )
                );

                if (pimple().finalNonOrthogonalIter())
                {
                    phi() -= pEqn.flux();
                }

                // Perform overset interpolation (after flux reconstruction)
                oversetFvPatchScalarField::oversetInterpolate(p());
            }

            fluidModel::oversetContinuityErrs();

            // Explicitly relax pressure for momentum corrector
            // except for last corrector
            if (!pimple().finalIter())
            {
                p().relax();
            }

            U() -= rAU*fvc::grad(p());
            U().correctBoundaryConditions();
            oversetFvPatchVectorField::oversetInterpolate(U());

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi(), U());
        }

        turbulence_->correct();
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

#endif

// ************************************************************************* //
