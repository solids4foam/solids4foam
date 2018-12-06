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

#include "mechanicalEnergies.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanicalEnergies, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mechanicalEnergies::mechanicalEnergies
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    linearBulkViscosityCoeff_
    (
        dict.lookupOrDefault<scalar>
        (
            "linearBulkViscosityCoeff", 0.06
        )
    ),
    viscousPressurePtr_(),
    epsilonVolPtr_(),
    externalWorkField_
    (
        IOobject
        (
            "externalWorkField",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy, 0.0)
    ),
    internalEnergyField_
    (
        IOobject
        (
            "internalEnergyField",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy, 0.0)
    ),
    bulkViscosityEnergyField_
    (
        IOobject
        (
            "bulkViscosityEnergyField",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy, 0.0)
    )
{
    externalWorkField_.oldTime();
    internalEnergyField_.oldTime();
    bulkViscosityEnergyField_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceScalarField& mechanicalEnergies::viscousPressure
(
    const volScalarField& rho,
    const surfaceScalarField& waveSpeed,
    const volTensorField& gradD
)
{
    if (viscousPressurePtr_.empty())
    {
        viscousPressurePtr_.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "viscousPressure",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimPressure, 0.0)
            )
        );
    }

    viscousPressurePtr_() = linearBulkViscosityCoeff_*fvc::interpolate
    (
        rho*fvc::ddt(epsilonVol(gradD))
    )*waveSpeed/mesh_.deltaCoeffs();

    return viscousPressurePtr_();
}


const volScalarField& mechanicalEnergies::epsilonVol
(
    const volTensorField& gradD
)
{
    if (epsilonVolPtr_.empty())
    {
        epsilonVolPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "epsilonVol",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    epsilonVolPtr_() = tr(gradD)/3.0;

    return epsilonVolPtr_();
}

void mechanicalEnergies::checkEnergies
(
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& DD,
    const volSymmTensorField& sigma,
    const volTensorField& gradD,
    const volTensorField& gradDD,
    const surfaceScalarField& waveSpeed
)
{
    // Calculate kinetic energy
    const scalar kineticEnergy =
        gSum(0.5*rho.internalField()*mesh_.V()*(U & U));

    // Calculate internal energy
    internalEnergyField_.internalField() =
        internalEnergyField_.oldTime().internalField()
      + (sigma && symm(gradDD))*mesh_.V();

    const scalar internalEnergy = gSum(internalEnergyField_.internalField());

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
                        mesh_.Sf().boundaryField()[patchI]
                      & sigma.boundaryField()[patchI]
                    )
                  & DD.boundaryField()[patchI]
                );

            externalWork += sum(externalWorkField_.boundaryField()[patchI]);
        }
    }

    // Sync in parallel
    reduce(externalWork, sumOp<scalar>());


    // Calculate energy dissipated due to linear bulk viscosity term
    bulkViscosityEnergyField_.internalField() =
        bulkViscosityEnergyField_.oldTime().internalField()
      + (
            fvc::reconstruct(viscousPressurePtr_()*mesh_.Sf()) && gradDD
        )*mesh_.V();

    const scalar bulkViscosityEnergy =
        gSum(bulkViscosityEnergyField_.internalField());

    // Calculate energy dissipated due to velocity damping term
    // To-do
    //eta*rho()*fvm::ddt(D())

    // Calculate energy generated due to gravity forces
    // To-do

    // Calculate the energy imbalance
    // This should stay less than 1% of the max energy component
    const scalar energyImbalance =
        externalWork - internalEnergy - kineticEnergy - bulkViscosityEnergy;
    const scalar energyImbalancePercent =
        100.0*mag(energyImbalance)/max
        (
            SMALL, max(externalWork, max(internalEnergy, kineticEnergy))
        );

    Info<< "External work = " << externalWork << " J" << nl
        << "Internal energy = " << internalEnergy << " J" << nl
        << "Kinetic energy = " << kineticEnergy << " J" << nl
        << "Bulk viscosity energy = " << bulkViscosityEnergy << " J" << nl
        << "Energy imbalance = " << energyImbalance << " J" << nl
        << "Energy imbalance (% of max) = " << energyImbalancePercent << " %"
        << endl;

    if (energyImbalancePercent > 10.0)
    {
       WarningIn(type() + "::checkEnergies()")
           << "The energy imbalance is greater than 10%" << endl;
       // FatalErrorIn(type() + "::checkEnergies()")
       //     << "The energy imbalance is greater than 50%"
       //     << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
