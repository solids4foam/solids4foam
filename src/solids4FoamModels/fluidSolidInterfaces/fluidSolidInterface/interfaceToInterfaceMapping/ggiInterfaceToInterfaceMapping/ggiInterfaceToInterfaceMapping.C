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

#include "ggiInterfaceToInterfaceMapping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ggiInterfaceToInterfaceMapping, 0);
addToRunTimeSelectionTable
(
    interfaceToInterfaceMapping, ggiInterfaceToInterfaceMapping, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ggiInterfaceToInterfaceMapping::makeInterpolator() const
{
    if (interpolatorPtr_.valid())
    {
        FatalErrorIn
        (
            "void ggiInterfaceToInterfaceMapping::calcGgiInterpolators() const"
        )   << "Pointer is already set!"
            << abort(FatalError);
    }

    Info<< "Create GGI zone-to-zone interpolator" << endl;

    interpolatorPtr_.set
    (
        new GGIInterpolation<standAlonePatch, standAlonePatch>
        (
            zoneA(),
            zoneB(),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            true,           // Patch data is complete on all processors
            SMALL,          // Non-overlapping face tolerances
            SMALL,
            true,           // Rescale weighting factors
            ggiInterpolation::BB_OCTREE
        )
    );

    checkZoneAToZoneBError();
    checkZoneBToZoneAError();

    Info<< "Number of uncovered master faces: "
        << interpolatorPtr_().uncoveredMasterFaces().size() << nl
        << "Number of uncovered slave faces: "
        << interpolatorPtr_().uncoveredSlaveFaces().size() << nl << endl;

    // Force point distance calculation
    //interpolatorPtr_().slavePointDistanceToIntersection();
    //interpolatorPtr_().masterPointDistanceToIntersection();
}

const GGIInterpolation<standAlonePatch, standAlonePatch>&
ggiInterfaceToInterfaceMapping::interpolator() const
{
    if (interpolatorPtr_.empty())
    {
        makeInterpolator();
    }

    return interpolatorPtr_();
}


void ggiInterfaceToInterfaceMapping::checkZoneAToZoneBError() const
{
    // Reference to patch face centres
    const vectorField& patchAFaceCentres = patchA().faceCentres();

    // Construct global zone field
    const vectorField zoneAFaceCentres =
        globalPatchA().patchFaceToGlobal(patchAFaceCentres);

    // Interpolate global zone field from A to B
    const vectorField zoneBFaceCentres =
        interpolator().masterToSlave(zoneAFaceCentres);

    // Extract local patch field
    const vectorField patchBFaceCentres =
        globalPatchB().globalFaceToPatch(zoneBFaceCentres);

    // Print maximum error
    Info<< "interface-to-interface face error: "
        << gMax(mag(patchBFaceCentres - patchB().faceCentres()))
        << endl;

}


void ggiInterfaceToInterfaceMapping::checkZoneBToZoneAError() const
{
    const vectorField zoneBPointsAtFluid =
        interpolator().slaveToMasterPointInterpolate(zoneB().localPoints());
            //(currentSolidZonesPoints()[interfaceI]); // double-check

    const vectorField& zoneAPoints = zoneA().localPoints();

    Info<< "interface-to-interface point error: "
        << gMax(mag(zoneAPoints - zoneBPointsAtFluid))
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ggiInterfaceToInterfaceMapping::ggiInterfaceToInterfaceMapping
(
    const word& type,
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalPatchA,
    const globalPolyPatch& globalPatchB
)
:
    interfaceToInterfaceMapping
    (
        type, dict, patchA, patchB, globalPatchA, globalPatchB
    ),
    interpolatorPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ggiInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFacesZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void ggiInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPointsZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void ggiInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFacesZoneToZone<vector>(fromZone, toZone, fromField, toField);
}


void ggiInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferPointsZoneToZone<vector>(fromZone, toZone, fromField, toField);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceToInterfaceMappings

} // End namespace Foam

// ************************************************************************* //
