/*---------------------------------------------------------------------------*\
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

#include "explicitNonLinGeomTotalLagTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "patchCorrectionVectors.H"


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
    solidModel, explicitNonLinGeomTotalLagTotalDispSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void explicitNonLinGeomTotalLagTotalDispSolid::updateStress()
{
    // Increment of displacement
    DD() = D() - D().oldTime();

    // Interpolate D to pointD
    mechanical().interpolate(D(), pointD(), false);

    // Update gradient of displacement
    // mechanical().grad(D(), pointD(), gradD(), gradDf_);
    mechanical().grad(D(), gradD());

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Total deformation gradient
    F_ = I + gradD().T();
    // Ff_ = I + gradDf_.T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);
    // Finvf_ = inv(Ff_);

    // Jacobian of the deformation gradient
    J_ = det(F_);
    // Jf_ = det(Ff_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
    // mechanical().correct(sigmaf_);
}


void explicitNonLinGeomTotalLagTotalDispSolid::correctDOnTractionBoundaries()
{
    forAll(D().boundaryField(), patchI)
    {
        if (D().boundaryField()[patchI].type() == "solidTraction")
        {
            // Non-orthogonal correction vectors
            const vectorField k
            (
                patchCorrectionVectors(mesh().boundary()[patchI])
            );

#ifdef OPENFOAMESIORFOUNDATION
            D().boundaryFieldRef()[patchI] =
                D().boundaryField()[patchI].patchInternalField()
              + (k & gradD().boundaryField()[patchI].patchInternalField());
#else
            D().boundaryField()[patchI] =
                D().boundaryField()[patchI].patchInternalField()
              + (k & gradD().boundaryField()[patchI].patchInternalField());
#endif
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitNonLinGeomTotalLagTotalDispSolid::explicitNonLinGeomTotalLagTotalDispSolid
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
    Ff_
    (
        IOobject
        (
            "Ff",
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
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    Finvf_
    (
        IOobject
        (
            "Finvf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(Ff_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    Jf_
    (
        IOobject
        (
            "Jf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(Ff_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
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
            IOobject::READ_IF_PRESENT,
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

    // Force all required old-time fields to be created
    // fvc::d2dt2(D());

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        DD() = D() - D().oldTime();
        mechanical().grad(D(), gradD());
        mechanical().grad(D(), pointD(), gradD(), gradDf_);
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();
        Ff_ = I + gradDf_.T();
        Finv_ = inv(F_);
        Finvf_ = inv(Ff_);
        J_ = det(F_);
        Jf_ = det(Ff_);

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void explicitNonLinGeomTotalLagTotalDispSolid::setDeltaT(Time& runTime)
{
    // waveSpeed = cellWidth/deltaT
    waveSpeed_ = fvc::interpolate(Foam::sqrt(impK_/rho()));

    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // stableTimeStep/sqrt(2): we will default to this value
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
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    // Update print info
    physicsModel::printInfo() = false;
    // bool
    // (
    //     runTime.timeIndex() % infoFrequency() == 0
    //  || mag(runTime.value() - runTime.endTime().value()) < SMALL
    // );

    // if (physicsModel::printInfo())
    if (time().timeIndex() == 1)
    {
        Info<< nl << "Setting deltaT = " << newDeltaT
            << ", maxCo = " << maxCo << endl;
    }

    runTime.setDeltaT(newDeltaT);
}


bool explicitNonLinGeomTotalLagTotalDispSolid::evolve()
{
    // Mesh update loop
    do
    {
        if (time().timeIndex() == 1)
        {
            Info<< "Solving the solid momentum equation for D" << nl
                << "Simulation Time, Clock Time, Max Stress" << endl;
            
        }

        physicsModel::printInfo() = bool
        (
            time().timeIndex() % infoFrequency() == 0
         || mag(time().value() - time().endTime().value()) < SMALL
        );

        if (physicsModel::printInfo())
        {
            Info<< time().value() << " " << time().elapsedClockTime()
                << " " << max(mag(sigma())).value() << endl;

            physicsModel::printInfo() = false;
        }

        // Central difference scheme

        const dimensionedScalar& deltaT = time().deltaT();
        //const dimensionedScalar& deltaT0 = time().deltaT0();

        const bool firstOrder = true;

        if (firstOrder)
        {
            // Compute the velocity
            // Note: this is the velocity at the middle of the time-step
            // U() = U().oldTime() + 0.5*(deltaT + deltaT0)*a_.oldTime();
            U() = U().oldTime() + deltaT*a_.oldTime();

            // Compute displacement
            D() = D().oldTime() + deltaT*U();

            // Enforce boundary conditions on the displacement field
            D().correctBoundaryConditions();

            // Zero-gradient extrapolation on traction boundaries
            // Do I also need to fix the boundary gradient on gradD?
            //correctDOnTractionBoundaries();

            // Update the stress field based on the latest D field
            updateStress();

            // Compute acceleration
            // Note the inclusion of a linear bulk viscosity pressure term to
            // dissipate high frequency energies, and a Rhie-Chow or JST term to
            // suppress checker-boarding
    #ifdef OPENFOAMESIORFOUNDATION
            a_.primitiveFieldRef() =
    #else
            a_.internalField() =
    #endif
                (
                    fvc::div
                    (
                        // (Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_
                        // (mesh().Sf() & fvc::interpolate(J_*Finv_ & sigma()))
                        mesh().Sf()
                      & (
                            fvc::interpolate(J_)
                           *fvc::interpolate(Finv_)
                          & fvc::interpolate(sigma())
                        )
                      // + mesh().Sf()*energies_.viscousPressure
                      //   (
                      //       rho(), waveSpeed_, gradD()
                      //   )
                    )().internalField()
                  + stabilisation().stabilisation
                    (
                        D(), gradD(), impK_
                    )().internalField()
                )/rho().internalField()
    #ifdef OPENFOAMESIORFOUNDATION
              + g();
    #else
              + g().value();
    #endif

            a_.correctBoundaryConditions();
        }
        else
        {
            // Second-order two-stage Runge-Kutta as described in Lee et al 2013

            // Compute the intermediate velocity
            U() = U().oldTime() + deltaT*a_.oldTime();

            // Compute the intermediate displacement
            D() = D().oldTime() + deltaT*U().oldTime();

            // Enforce boundary conditions on the displacement field
            D().correctBoundaryConditions();

            // Zero-gradient extrapolation on traction boundaries
            // Do I also need to fix the boundary gradient on gradD?
            correctDOnTractionBoundaries();

            // Update the stress field based on the intermediate D field
            updateStress();

            // Compute the acceleration based on the intermediate D and stress
    #ifdef OPENFOAMESIORFOUNDATION
            a_.primitiveFieldRef() =
    #else
            a_.internalField() =
    #endif
                (
                    fvc::div
                    (
                        (mesh().Sf() & fvc::interpolate(J_*Finv_ & sigma()))
                    )().internalField()
                  + stabilisation().stabilisation
                    (
                        D(), gradD(), impK_
                    )().internalField()
                )/rho().internalField()
    #ifdef OPENFOAMESIORFOUNDATION
              + g();
    #else
              + g().value();
    #endif

            // Does a need non-orthogonal corrections at the boundary?
            a_.correctBoundaryConditions();

            // Compute the final displacement
            D() = 0.5*(D().oldTime() + D() + deltaT*U());

            // Zero-gradient extrapolation on traction boundaries
            correctDOnTractionBoundaries();

            // Compute the final velocity
            U() = 0.5*(U().oldTime() + U() + deltaT*a_);
        }

        // Check energies
        // energies_.checkEnergies
        // (
        //     rho(), U(), D(), DD(), sigma(), gradD(), gradDD(), waveSpeed_, g(),
        //     0.0, impKf_, physicsModel::printInfo()
        // );
    }
    while (mesh().update());

    if (time().outputTime())
    {
        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();
    }

    // Stop if stress explodes
    if (isnan(max(mag(a_)).value()))
    {
        FatalError
            << "max acceleration is NaN!" << abort(FatalError);
    }

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> explicitNonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
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
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
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

    // Test
    // const scalarField& J = J_.boundaryField()[patchID];
    // const scalarField SoS0(mag(J*(n & Finv)));
    // vectorField nCurrent(n & Finv);
    // nCurrent /= mag(nCurrent);
    // const vectorField tractionCauchy(traction - pressure*nCurrent);
    // const vectorField traction2PK((Finv & tractionCauchy)*SoS0);

    // // Return patch snGrad
    // return tmp<vectorField>
    // (
    //     new vectorField
    //     (
    //         (
    //             traction2PK - (n & (J*transform(Finv, pSigma)))
    //           + impK*(n & pGradD)
    //         )*rImpK
    //     )
    // );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
