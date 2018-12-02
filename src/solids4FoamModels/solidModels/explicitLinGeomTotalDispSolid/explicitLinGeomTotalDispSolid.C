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
    rImpK_(1.0/impK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void explicitLinGeomTotalDispSolid::setDeltaT(Time& runTime)
{
    // For now, only the linearElastic material law is allowed
    const PtrList<mechanicalLaw>& mechLaws = mechanical();
    if (mechLaws.size() != 1)
    {
        FatalErrorIn(type() + "::" + type())
            << type() << " can currently only be used with a single material"
            << "\nConsider using one of the other solidModels."
            << abort(FatalError);
    }
    else if (!isA<linearElastic>(mechLaws[0]))
    {
        FatalErrorIn(type() + "::" + type())
            << type() << " can only be used with the linearElastic "
            << "mechanicalLaw" << nl
            << abort(FatalError);
    }

    // Cast the mechanical law to a linearElastic mechanicalLaw
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Calculate the speed of sound
    const scalar waveVelocity  =
        Foam::sqrt
        (
            mech.E()*(1 - mech.nu())
           /(mech.rhoScalar()*(1 + mech.nu())*(1 - 2*mech.nu()))
        ).value();

    // Courant number
    const scalarField Co =
        waveVelocity*runTime.deltaT().value()
        *mesh().surfaceInterpolation::deltaCoeffs().internalField();

    // Calculate required time-step for a Courant number of 1.0
    const scalar requiredDeltaT =
        1.0/
        gMax
        (
            mesh().surfaceInterpolation::deltaCoeffs().internalField()
           *waveVelocity
        );

    //const scalar averageCo = gAverage(Co);
    //const scalar maxCo = gMax(Co);
    //const scalar averageWaveVel = waveVelocity;
    //const scalar maxWaveVel = waveVelocity;

    const scalar newDeltaT =
        readScalar(runTime.controlDict().lookup("maxCo"))*requiredDeltaT;

//     Info<< "Courant Number\n\tmean: " << averageCo
//         << "\n\tmax: " << maxCo << nl
//         << "Wave velocity magnitude\n\tmean " << averageWaveVel
//         << "\n\tmax: " << maxWaveVel << nl
//         << "deltaT for a maximum Courant number of 1.0: "
//         << requiredDeltaT << nl;
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

        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
          + eta*rho()*fvm::ddt(D())
         == fvc::div(sigma(), "div(sigma)")
          + rho()*g()
          + mechanical().RhieChowCorrection(D(), gradD())
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Solve the linear system
        DEqn.solve();

        // Print out damping
        Info<< "Max inertia: "
            << gMax(mag(rho()*fvc::d2dt2(D())().internalField()))
            << nl
            << "Max viscous damping: "
            << gMax(mag(rho()*eta*fvc::ddt(D())().internalField()))
            << endl;

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
