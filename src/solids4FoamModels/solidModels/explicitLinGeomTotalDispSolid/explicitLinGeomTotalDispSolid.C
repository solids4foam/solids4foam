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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitLinGeomTotalDispSolid, 0);
addToRunTimeSelectionTable(physicsModel, explicitLinGeomTotalDispSolid, solid);
addToRunTimeSelectionTable
(
    solidModel, explicitLinGeomTotalDispSolid, dictionary
);

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
    RhieChowScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "RhieChowScale", 0.0
        )
    ),
    eta_
    (
        solidModelDict().lookupOrDefault<dimensionedScalar>
        (
            "numericalViscosity",
            dimensionedScalar("eta", dimless/dimTime, 0.0)
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
    energies_(mesh(), solidModelDict())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitLinGeomTotalDispSolid::setDeltaT(Time& runTime)
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value
    const scalar requiredDeltaT =
        1.0/
        gMax
        (
            mesh().surfaceInterpolation::deltaCoeffs().internalField()
           *waveSpeed_.internalField()
        );

    Info<< "Min delta: "
        << 1.0/
        gMax
        (
            mesh().surfaceInterpolation::deltaCoeffs().internalField()
        ) << endl;

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.7071);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    Info<< "maxCo = " << maxCo << nl
        << "deltaT = " << requiredDeltaT << nl << endl;

    runTime.setDeltaT(newDeltaT);
}


bool explicitLinGeomTotalDispSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Do not write default linear solver residuals
    blockLduMatrix::debug = 0;

    // Mesh update loop
    do
    {
        Info<< "Solving the momentum equation for D" << endl;

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
          - eta_*rho()*fvm::ddt(D())
         == fvc::div(sigma(), "div(sigma)")
          + rho()*g()
          + RhieChowScaleFactor_*mechanical().RhieChowCorrection(D(), gradD())
          + fvc::div
            (
                energies_.viscousPressure(rho(), waveSpeed_, gradD())
               *mesh().Sf()
            )
        );

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
        energies_.checkEnergies
        (
            rho(), U(), DD(), sigma(), gradD(), gradDD(), waveSpeed_
        );
    }
    while (mesh().update());

    blockLduMatrix::debug = 1;

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
