/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "explicitUnsLinGeomTotalDispSolid.H"
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

defineTypeNameAndDebug(explicitUnsLinGeomTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitUnsLinGeomTotalDispSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void explicitUnsLinGeomTotalDispSolid::updateStress()
{
    // Update increment of displacement
    DD() = D() - D().oldTime();

    // Interpolate D to pointD
    mechanical().interpolate(D(), pointD(), false);

    // Update gradient of displacement
    mechanical().grad(D(), pointD(), gradD(), gradDf_);

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigmaf_);
    mechanical().correct(sigma());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Increment of displacement
    DD() = D() - D().oldTime();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitUnsLinGeomTotalDispSolid::explicitUnsLinGeomTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    // RhieChowScaleFactor_
    // (
    //     solidModelDict().lookupOrDefault<scalar>
    //     (
    //         "RhieChowScale", 0.0
    //     )
    // ),
    JSTScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "JSTScaleFactor", 0.01
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
    energies_(mesh(), solidModelDict()),
    a_
    (
        IOobject
        (
            "a",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero", dimVelocity/dimTime, vector::zero
        ),
        "zeroGradient"
    )
{
    DisRequired();

    a_.oldTime();
    U().oldTime();

    // Update stress
    updateStress();

    // Update initial acceleration
#ifdef OPENFOAMESIORFOUNDATION
    a_.primitiveFieldRef() =
#else
    a_.internalField() =
#endif
        fvc::div(sigma(), "div(sigma)")().internalField()
       /(rho().internalField());
    a_.correctBoundaryConditions();

    // Set the printInfo
    physicsModel::printInfo() = bool
    (
        runTime.timeIndex() % infoFrequency() == 0
     || mag(runTime.value() - runTime.endTime().value()) < SMALL
    );

    Info<< "Frequency at which info is printed: every " << infoFrequency()
        << " time-steps" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitUnsLinGeomTotalDispSolid::setDeltaT(Time& runTime)
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
#ifdef OPENFOAMESIORFOUNDATION
            DimensionedField<scalar, Foam::surfaceMesh>
#else
            Field<scalar>
#endif
            (
                mesh().surfaceInterpolation::deltaCoeffs().internalField()
               *waveSpeed_.internalField()
            )
        );

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.7071);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    // Update print info
    physicsModel::printInfo() = bool
    (
        runTime.timeIndex() % infoFrequency() == 0
     || mag(runTime.value() - runTime.endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< nl << "Setting deltaT = " << newDeltaT
            << ", maxCo = " << maxCo << endl;
    }

    runTime.setDeltaT(newDeltaT);
}


bool explicitUnsLinGeomTotalDispSolid::evolve()
{
    // Mesh update loop
    do
    {
        if (printInfo())
        {
            Info<< "Solving the solid momentum equation for D" << endl;
        }

        // Central difference scheme

        const dimensionedScalar& deltaT = time().deltaT();
        const dimensionedScalar& deltaT0 = time().deltaT0();

        // Compute the velocity
        // Note: this is the velocity at the middle of the time-step
        U() = U().oldTime() + 0.5*(deltaT + deltaT0)*a_.oldTime();

        // Compute displacement
        D() = D().oldTime() + deltaT*U();

        // Enforce boundary conditions on the displacement field
        D().correctBoundaryConditions();

        // Update the stress field based on the latest D field
        updateStress();

        // Compute acceleration
        // Note the inclusion of a linear bulk viscosity pressure term to
        // dissipate high frequency energies, and a Rhie-Chow term to avoid
        // checker-boarding
#ifdef OPENFOAMESIORFOUNDATION
        a_.primitiveFieldRef() =
#else
        a_.internalField() =
#endif
            (
                fvc::div
                (
                    (mesh().Sf() & sigmaf_)
                  + mesh().Sf()*energies_.viscousPressure
                    (
                        rho(), waveSpeed_, gradD()
                    )
                )().internalField()
              - JSTScaleFactor_*fvc::laplacian
                (
                    mesh().magSf(),
                    fvc::laplacian
                    (
                        0.5*(deltaT + deltaT0)*impKf_, U(), "laplacian(DU,U)"
                    ),
                    "laplacian(DU,U)"
                )().internalField()
            )/rho().internalField()
          + g().value();
        a_.correctBoundaryConditions();

        // Check energies
        energies_.checkEnergies
        (
            rho(), U(), D(), DD(), sigma(), gradD(), gradDD(), waveSpeed_, g(),
            0.0, impKf_, printInfo()
        );
    }
    while (mesh().update());

    return true;
}


tmp<vectorField> explicitUnsLinGeomTotalDispSolid::tractionBoundarySnGrad
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
    const vectorField n(patch.nf());

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
