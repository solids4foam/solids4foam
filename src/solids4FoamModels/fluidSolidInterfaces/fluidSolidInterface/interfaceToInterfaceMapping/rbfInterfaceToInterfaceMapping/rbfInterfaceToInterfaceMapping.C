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

#include "rbfInterfaceToInterfaceMapping.H"
#include "addToRunTimeSelectionTable.H"
#include "TPSFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rbfInterfaceToInterfaceMapping, 0);
addToRunTimeSelectionTable
(
    interfaceToInterfaceMapping,
    rbfInterfaceToInterfaceMapping,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void rbfInterfaceToInterfaceMapping::makeZoneAToZoneBInterpolator() const
{
    if (zoneAToZoneBInterpolatorPtr_ != NULL)
    {
        FatalErrorIn
        (
            "void rbfInterfaceToInterfaceMapping::"
            "makeZoneAToZoneBInterpolator() const"
        )   << "Pointer already set!"
            << abort(FatalError);
    }

    Info<< "Create RBF interpolator from " << globalPatchA().patchName()
        << " to " << globalPatchB().patchName() << endl;

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    zoneAToZoneBInterpolatorPtr_ =
        std::shared_ptr<RBFInterpolation>
        (
            new RBFInterpolation(rbfFunction)
        );

    const vectorField zoneBFaceCentres(zoneB().faceCentres());
    const vectorField zoneAFaceCentres(zoneA().faceCentres());

    matrix zoneAX(zoneAFaceCentres.size(), 3);
    matrix zoneBX(zoneBFaceCentres.size(), 3);

    forAll(zoneAFaceCentres, faceI)
    {
        zoneAX(faceI, 0) = zoneAFaceCentres[faceI].x();
        zoneAX(faceI, 1) = zoneAFaceCentres[faceI].y();
        zoneAX(faceI, 2) = zoneAFaceCentres[faceI].z();
    }

    forAll(zoneBFaceCentres, faceI)
    {
        zoneBX(faceI, 0) = zoneBFaceCentres[faceI].x();
        zoneBX(faceI, 1) = zoneBFaceCentres[faceI].y();
        zoneBX(faceI, 2) = zoneBFaceCentres[faceI].z();
    }

    zoneAToZoneBInterpolatorPtr_->compute(zoneAX, zoneBX);

    // Check interpolation error

    matrix zoneAXatZoneB(zoneBFaceCentres.size(), 3);

    zoneAToZoneBInterpolatorPtr_->interpolate(zoneAX, zoneAXatZoneB);

    vectorField zoneAFaceCentresAtZoneB(zoneBFaceCentres.size(), vector::zero);

    forAll(zoneAFaceCentresAtZoneB, faceI)
    {
        zoneAFaceCentresAtZoneB[faceI].x() = zoneAXatZoneB(faceI, 0);
        zoneAFaceCentresAtZoneB[faceI].y() = zoneAXatZoneB(faceI, 1);
        zoneAFaceCentresAtZoneB[faceI].z() = zoneAXatZoneB(faceI, 2);
    }

    const scalar maxDist = gMax
    (
        mag(zoneAFaceCentresAtZoneB - zoneBFaceCentres)
    );

    Info<< "    face interpolation error: " << maxDist << endl;
}


const std::shared_ptr<RBFInterpolation>&
rbfInterfaceToInterfaceMapping::zoneAToZoneBInterpolator() const
{
    if (zoneAToZoneBInterpolatorPtr_ == NULL)
    {
        makeZoneAToZoneBInterpolator();
    }

    return zoneAToZoneBInterpolatorPtr_;
}


void rbfInterfaceToInterfaceMapping::makeZoneBToZoneAInterpolator() const
{
    if (zoneBToZoneAInterpolatorPtr_ != NULL)
    {
        FatalErrorIn
        (
            "void rbfInterfaceToInterfaceMapping::"
            "makeZoneBToZoneAInterpolator() const"
        )   << "Pointer already set!"
            << abort(FatalError);
    }

    Info<< "Create RBF interpolator from " << globalPatchB().patchName()
        << " to " << globalPatchA().patchName() << endl;

    std::shared_ptr<RBFFunctionInterface> rbfFunction;
    rbfFunction = std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    zoneBToZoneAInterpolatorPtr_ =
        std::shared_ptr<RBFInterpolation>
        (
            new RBFInterpolation(rbfFunction)
        );

    const vectorField zoneBPoints = zoneB().localPoints();
    const vectorField zoneAPoints = zoneA().localPoints();

    matrix zoneAX(zoneAPoints.size(), 3);
    matrix zoneBX(zoneBPoints.size(), 3);

    forAll(zoneAPoints, faceI)
    {
        zoneAX(faceI, 0) = zoneAPoints[faceI].x();
        zoneAX(faceI, 1) = zoneAPoints[faceI].y();
        zoneAX(faceI, 2) = zoneAPoints[faceI].z();
    }

    forAll(zoneBPoints, faceI)
    {
        zoneBX(faceI, 0) = zoneBPoints[faceI].x();
        zoneBX(faceI, 1) = zoneBPoints[faceI].y();
        zoneBX(faceI, 2) = zoneBPoints[faceI].z();
    }

    zoneBToZoneAInterpolatorPtr_->compute(zoneBX, zoneAX);

    // Check interpolation error

    matrix zoneBXatZoneA(zoneAPoints.size(), 3);

    zoneBToZoneAInterpolatorPtr_->interpolate(zoneBX, zoneBXatZoneA);

    vectorField zoneAPointsInterp(zoneAPoints.size(), vector::zero);

    forAll(zoneAPoints, faceI)
    {
        zoneAPointsInterp[faceI].x() = zoneBXatZoneA(faceI, 0);
        zoneAPointsInterp[faceI].y() = zoneBXatZoneA(faceI, 1);
        zoneAPointsInterp[faceI].z() = zoneBXatZoneA(faceI, 2);
    }

    const scalar maxDist = gMax
    (
        mag(zoneAPointsInterp - zoneAPoints)
    );

    Info<< "    point interpolation error: " << maxDist << endl;
}


const std::shared_ptr<RBFInterpolation>&
rbfInterfaceToInterfaceMapping::zoneBToZoneAInterpolator() const
{
    if (zoneBToZoneAInterpolatorPtr_ == NULL)
    {
        makeZoneBToZoneAInterpolator();
    }

    return zoneBToZoneAInterpolatorPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rbfInterfaceToInterfaceMapping::rbfInterfaceToInterfaceMapping
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
    zoneAToZoneBInterpolatorPtr_(NULL),
    zoneBToZoneAInterpolatorPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rbfInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFacesZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void rbfInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPointsZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void rbfInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFacesZoneToZone<vector>(fromZone, toZone, fromField, toField);
}


void rbfInterfaceToInterfaceMapping::transferPointsZoneToZone
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
