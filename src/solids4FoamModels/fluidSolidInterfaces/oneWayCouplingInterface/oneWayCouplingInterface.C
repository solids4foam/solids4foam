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

#include "oneWayCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "oneWayFsiFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, oneWayCouplingInterface, fluidSolidInteraction
);
addToRunTimeSelectionTable
(
    fluidSolidInterface, oneWayCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneWayCouplingInterface::oneWayCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    solidZonesTraction_(nGlobalPatches())
{
    // Initialize zone traction fields
    forAll(fluid().globalPatches(), interfaceI)
    {
        solidZonesTraction_[interfaceI] =
            vectorField
            (
                solid().globalPatches()[interfaceI].globalPatch().size(),
                vector::zero
            );
    }

    if (!isA<fluidModels::oneWayFsiFluid>(fluid()))
    {
        FatalErrorIn
        (
            "oneWayCouplingInterface::oneWayCouplingInterface\n"
            "(\n"
            "    Time& runTime,\n"
            "    const word& region\n"
            ")\n"
        )   << "The " << type() << " FSI coupling can only be used with the "
            << fluidModels::oneWayFsiFluid::typeName << " fluidModel"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool oneWayCouplingInterface::evolve()
{
    fluidSolidInterface::initializeFields();

    updateInterpolatorAndGlobalPatches();

    fluid().evolve();

    updateTraction();

    solid().evolve();

    solid().updateTotalFields();

    return 0;
}


void oneWayCouplingInterface::updateTraction()
{
    Info<< "Update traction on solid patch/patches" << endl;

    forAll(fluid().globalPatches(), interfaceI)
    {
        const vectorField fluidZoneTraction
        (
            fluid().faceZoneViscousForce(interfaceI)
        );

        const scalarField fluidZonePressure
        (
            fluid().faceZonePressureForce(interfaceI)
        );

        // Calculate fluid traction
        const vectorField& p =
            fluid().globalPatches()[interfaceI].globalPatch().localPoints();

        const faceList& f =
            fluid().globalPatches()[interfaceI].globalPatch().localFaces();

        vectorField n(f.size(), vector::zero);
        forAll(n, faceI)
        {
            n[faceI] = f[faceI].normal(p);
            n[faceI] /= mag(n[faceI]);
        }

        const vectorField fluidZonesTraction
        (
            fluidZoneTraction - fluidZonePressure*n
        );

        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluid().globalPatches()[interfaceI].globalPatch(), // fromZone
            solid().globalPatches()[interfaceI].globalPatch(), // toZone
            fluidZonesTraction,                                // fromField
            solidZonesTraction_[interfaceI]                    // toField
        );

        // Flip the sign
        solidZonesTraction_[interfaceI] = - solidZonesTraction_[interfaceI];

        if (coupled())
        {
            solid().setTraction
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZonesTraction_[interfaceI]
            );
        }

        // Print total force on solid and fluid interfaces
        Info<< "Total force on fluid interface " << interfaceI << ": "
            << fluidSolidInterface::totalForceOnInterface
            (
                fluid().globalPatches()[interfaceI].globalPatch(),
                fluidZonesTraction
            )
            << nl
            << "Total force on solid interface " << interfaceI << ": "
            << fluidSolidInterface::totalForceOnInterface
            (
                solid().globalPatches()[interfaceI].globalPatch(),
                solidZonesTraction_[interfaceI]
            ) << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
