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

#include "AitkenCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(AitkenCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, AitkenCouplingInterface, fluidSolidInteraction
);
addToRunTimeSelectionTable
(
    fluidSolidInterface, AitkenCouplingInterface, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AitkenCouplingInterface::AitkenCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    ),
    predictSolid_(fsiProperties().lookupOrDefault<bool>("predictSolid", true)),
    aitkenRelaxationFactor_(relaxationFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool AitkenCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    scalar residualNorm = 0;

    if (predictSolid_)
    {
        updateForce();

        solid().evolve();

        residualNorm =
            updateResidual();
    }

    do
    {
        outerCorr()++;

        // Transfer the displacement from the solid to the fluid
        updateDisplacement();

        // Move the fluid mesh
        moveFluidMesh();

        // Solve fluid
        fluid().evolve();

        // Transfer the force from the fluid to the solid
        updateForce();

        // Solve solid
        solid().evolve();

        // Calculate the FSI residual
        residualNorm = updateResidual();
    }
    while (residualNorm > outerCorrTolerance() && outerCorr() < nOuterCorr());

    solid().updateTotalFields();

    return 0;
}


void AitkenCouplingInterface::updateDisplacement()
{
    Info<< nl << "Time = " << fluid().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

    if (outerCorr() < 3)
    {
        Info<< "Current fsi under-relaxation factor: "
            << relaxationFactor_ << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        fluidZonePointsDispl() += relaxationFactor_*residual();
    }
    else
    {
        aitkenRelaxationFactor_ =
            -aitkenRelaxationFactor_
           *(
                sum
                (
                    residualPrev()
                  & (residual() - residualPrev())
                )
               /(
                    sum
                    (
                        (residual() - residualPrev())
                      & (residual() - residualPrev())
                    )
                )
            );

        if (Pstream::parRun())
        {
            if (!Pstream::master())
            {
                aitkenRelaxationFactor_ = 0.0;
            }

            // Pass to all procs
            reduce(aitkenRelaxationFactor_, sumOp<scalar>());
        }

        aitkenRelaxationFactor_ = mag(aitkenRelaxationFactor_);

        if (aitkenRelaxationFactor_ > 1)
        {
            // PC: in this case, would 1.0 be a better option?
            // Of course, the current option is more more stable
            aitkenRelaxationFactor_ = relaxationFactor_;
        }

        Info<< "Current fsi under-relaxation factor (Aitken): "
            << aitkenRelaxationFactor_ << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        fluidZonePointsDispl() += aitkenRelaxationFactor_*residual();
    }

    // Update movingWallPressure boundary conditions, if found
    fluidSolidInterface::updateMovingWallPressureAcceleration();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    fluidSolidInterface::syncFluidZonePointsDispl(fluidZonePointsDispl());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
