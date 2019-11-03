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
    forAll(fluid().globalPatches(), i)
    {
        solidZonesTraction_[i] =
            vectorField
            (
                solid().globalPatches()[i].globalPatch().size(),
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

    const List<tmp<vectorField> > faceZonesViscousForce
    (
        fluid().faceZonesViscousForce()
    );

    const List<tmp<scalarField> > faceZonesPressureForce
    (
        fluid().faceZonesPressureForce()
    );

    List<vectorField> fluidZonesTractionAtSolid
    (
        nGlobalPatches(), vectorField()
    );

    List<vectorField> fluidZonesTraction
    (
        nGlobalPatches(), vectorField()
    );

    forAll(fluid().globalPatches(), i)
    {
        const vectorField& fluidZoneTraction = faceZonesViscousForce[i]();

        const scalarField& fluidZonePressure = faceZonesPressureForce[i]();

        // Calculate fluid traction
        const vectorField& p =
            fluid().globalPatches()[i].globalPatch().localPoints();

        const faceList& f =
            fluid().globalPatches()[i].globalPatch().localFaces();

        vectorField n(f.size(), vector::zero);
        forAll(n, faceI)
        {
            n[faceI] = f[faceI].normal(p);
            n[faceI] /= mag(n[faceI]);
        }

        fluidZonesTraction[i] =
            fluidZoneTraction - fluidZonePressure*n;

        fluidZonesTractionAtSolid[i] =
            vectorField
            (
                solid().globalPatches()[i].globalPatch().size(),
                vector::zero
            );

        fluidZonesTractionAtSolid[i] =
            ggiInterpolators()[i].masterToSlave
            (
                -fluidZonesTraction[i]
            );

        solidZonesTraction_[i] = fluidZonesTractionAtSolid[i];
    }

    if (coupled())
    {
        solid().setTraction
        (
            solidPatchIndices(),
            solidZonesTraction_
        );
    }

    // Total force at the fluid side of the interface
    {
        forAll(fluid().globalPatches(), i)
        {
            const vectorField& p =
                fluid().globalPatches()[i].globalPatch().localPoints();
            const faceList& f =
                fluid().globalPatches()[i].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
            }

            const vector totalTractionForce = sum(fluidZonesTraction[i]*mag(S));

            Info<< "Total force on interface patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << " (fluid) = " << totalTractionForce << endl;
        }
    }

    // Total force at the solid side of the interface
    {
        forAll(solid().globalPatches(), i)
        {
            const vectorField& p =
                solid().globalPatches()[i].globalPatch().localPoints();
            const faceList& f =
                solid().globalPatches()[i].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
            }

            const vector totalTractionForce = sum(fluidZonesTractionAtSolid[i]*mag(S));

            Info<< "Total force on interface patch "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[i].patch().index()
                   ].name()
                << " (solid) = " << totalTractionForce << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
