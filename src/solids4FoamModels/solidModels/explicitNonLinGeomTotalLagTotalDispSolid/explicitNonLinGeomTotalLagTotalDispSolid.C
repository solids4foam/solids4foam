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

#include "explicitNonLinGeomTotalLagTotalDispSolid.H"
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

defineTypeNameAndDebug(explicitNonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, explicitNonLinGeomTotalLagTotalDispSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, explicitNonLinGeomTotalLagTotalDispSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitNonLinGeomTotalLagTotalDispSolid::
explicitNonLinGeomTotalLagTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
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

void explicitNonLinGeomTotalLagTotalDispSolid::setDeltaT(Time& runTime)
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
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    Info<< "Setting deltaT to " << requiredDeltaT << nl << endl;

    runTime.setDeltaT(newDeltaT);
}


    bool explicitNonLinGeomTotalLagTotalDispSolid::evolve()
{
    blockLduMatrix::debug = 0;

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    // Momentum equation total displacement total Lagrangian form
    // To-do: check viscousPressure is correct for finite strain
    fvVectorMatrix DEqn
    (
        rho()*fvm::d2dt2(D())
      - eta_*rho()*fvm::ddt(D())
     == fvc::div(J_*Finv_ & sigma(), "div(sigma)")
      + rho()*g()
      + RhieChowScaleFactor_*mechanical().RhieChowCorrection(D(), gradD())
      + fvc::div
        (
            energies_.viscousPressure(rho(), waveSpeed_, gradD())*mesh().Sf()
        )
    );

    // Solve the linear system
    DEqn.solve();

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Update the momentum equation inverse diagonal field
    // This may be used by the mechanical law when calculating the
    // hydrostatic pressure
    const volScalarField DEqnA("DEqnA", DEqn.A());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    // Check energies
    energies_.checkEnergies
    (
        rho(), U(), DD(), sigma(), gradD(), gradDD(), waveSpeed_
    );

    blockLduMatrix::debug = 1;

    return true;
}


tmp<vectorField> explicitNonLinGeomTotalLagTotalDispSolid::
tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    vectorField nCurrent = Finv.T() & n;
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
