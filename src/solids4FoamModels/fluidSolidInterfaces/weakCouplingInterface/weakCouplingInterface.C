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

#include "weakCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weakCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, weakCouplingInterface, fluidSolidInteraction
);
addToRunTimeSelectionTable
(
    fluidSolidInterface, weakCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

weakCouplingInterface::weakCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool weakCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    solid().evolve();

    updateWeakDisplacement();

    moveFluidMesh();

    fluid().evolve();

    updateForce();

    solid().updateTotalFields();

    return 0;
}


void weakCouplingInterface::initializeFields()
{
    // forAll(solid().globalPatches(), interfaceI)
    // {
    //     predictedSolidZonesTractionPtrList_[interfaceI] =
    //         vectorField
    //         (
    //             solid().globalPatches()[interfaceI].globalPatch().size(),
    //             vector::zero
    //         );
    // }

    fluidSolidInterface::initializeFields();
}


void weakCouplingInterface::updateWeakDisplacement()
{
    // Update the residual
    updateResidual();

    forAll(fluid().globalPatches(), interfaceI)
    {
        fluidZonesPointsDisplsPrev()[interfaceI] =
            fluidZonesPointsDispls()[interfaceI];

        // No under-relaxation
        fluidZonesPointsDispls()[interfaceI] += residuals()[interfaceI];
    }

    // Update movingWallPressure boundary conditions, if found
    fluidSolidInterface::updateMovingWallPressureAcceleration();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    fluidSolidInterface::syncFluidZonePointsDispl(fluidZonesPointsDispls());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
