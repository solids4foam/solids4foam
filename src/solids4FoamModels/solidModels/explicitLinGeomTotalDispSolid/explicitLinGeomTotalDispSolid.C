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

#include "explicitLinGeomTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "linearElastic.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitLinGeomTotalDispSolid, 0);
addToRunTimeSelectionTable(physicsModel, explicitLinGeomTotalDispSolid, solid);
addToRunTimeSelectionTable(solidModel, explicitLinGeomTotalDispSolid, dictionary);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void explicitLinGeomTotalDispSolid::checkEnergies()
{
    // Calculate kinetic energy
    const scalar kineticEnergy =
        gSum(0.5*rho().internalField()*mesh().V()*(U() & U()));

    // Calculate external work energy
    scalar externalWork = 0.0;
    forAll(externalWorkField_.boundaryField(), patchI)
    {
        if (!externalWorkField_.boundaryField()[patchI].coupled())
        {
            externalWorkField_.boundaryField()[patchI] =
                externalWorkField_.oldTime().boundaryField()[patchI]
              + (
                    (
                        mesh().Sf().boundaryField()[patchI]
                      & sigma().boundaryField()[patchI]
                    )
                  & DD().boundaryField()[patchI]
                );

            externalWork += sum(externalWorkField_.boundaryField()[patchI]);
        }
    }

    // Sync in parallel
    reduce(externalWork, sumOp<scalar>());

    // Calculate internal energy
    internalEnergyField_.internalField() =
        internalEnergyField_.oldTime().internalField()
      + (sigma() && symm(gradDD()))*mesh().V();

    const scalar internalEnergy = gSum(internalEnergyField_.internalField());

    // Calculate the net energy
    // This should stay less than 1% of the max energy component
    const scalar netEnergy = externalWork - internalEnergy - kineticEnergy;
    const scalar netEnergyPercent =
        100.0*mag(netEnergy)/max
        (
            SMALL, max(externalWork, max(internalEnergy, kineticEnergy))
        );

    Info<< "External work = " << externalWork << " J" << nl
        << "Internal energy = " << internalEnergy << " J" << nl
        << "Kinetic energy = " << kineticEnergy << " J" << nl
        << "Net energy = " << netEnergy << " J" << nl
        << "Net energy (% of max) = " << netEnergyPercent << " %" << endl;

//     if (netEnergyPercent > 50)
//     {
//         FatalErrorIn(type() + "::checkEnergies()")
//             << "The energy imbalance is greater than 50%" << abort(FatalError);
//     }

    if (netEnergyPercent > 10)
    {
        WarningIn(type() + "::checkEnergies()")
            << "The energy imbalance is greater than 10%" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitLinGeomTotalDispSolid::explicitLinGeomTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    linearBulkViscosityCoeff_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "linearBulkViscosityCoeff", 0.06
        )
    ),
    waveSpeed_
    (
        IOobject
        (
            "waveSpeed",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(Foam::sqrt(impK_/rho()))
    ),
    epsilonVol_
    (
        IOobject
        (
            "epsilonVol",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    externalWorkField_
    (
        IOobject
        (
            "externalWorkField",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimEnergy, 0.0)
    ),
    internalEnergyField_
    (
        IOobject
        (
            "internalEnergyField",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimEnergy, 0.0)
    )
{
    epsilonVol_.oldTime();
    externalWorkField_.oldTime();
    internalEnergyField_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void explicitLinGeomTotalDispSolid::setDeltaT(Time& runTime)
{
    const scalar requiredDeltaT =
        1.0/
        gMin
        (
            mesh().surfaceInterpolation::deltaCoeffs().internalField()
           *waveSpeed_.internalField()
        );

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.7);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    Info<< "Setting deltaT to " << requiredDeltaT << nl << endl;

    runTime.setDeltaT(newDeltaT);
}


bool explicitLinGeomTotalDispSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        blockLduMatrix::debug = 0;

        Info<< "Solving the momentum equation for D" << endl;

        // Stablisation viscosity
        const dimensionedScalar eta =
            solidModelDict().lookupOrDefault<dimensionedScalar>
            (
                "numericalViscosity",
                dimensionedScalar("eta", dimless/dimTime, 0.0)
            );

        // Calculate linear bulk damping e.g. as used by Abaqus/explicit
        epsilonVol_ = tr(gradD())/3.0;
        const surfaceScalarField viscousPressure =
            linearBulkViscosityCoeff_*fvc::interpolate
            (
                rho()*fvc::ddt(epsilonVol_)
            )*waveSpeed_/mesh().deltaCoeffs();

        // Scale factor for Rhie-Chow smoothing term
        const scalar RhieChowScale =
            solidModelDict().lookupOrDefault<scalar>
            (
                "RhieChowScale", 0.0
            );

        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
          + eta*rho()*fvm::ddt(D())
         == fvc::div(sigma(), "div(sigma)")
          + rho()*g()
          + RhieChowScale*mechanical().RhieChowCorrection(D(), gradD())
          + fvc::div(viscousPressure*mesh().Sf())
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Solve the linear system
        DEqn.solve();

        // Update increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());

        // Check energies
        checkEnergies();
    }
    while (mesh().update());

    return true;
}


tmp<vectorField> explicitLinGeomTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
