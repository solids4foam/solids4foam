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

#include "directMapInterfaceToInterfaceMapping.H"
#include "addToRunTimeSelectionTable.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(directMapInterfaceToInterfaceMapping, 0);
addToRunTimeSelectionTable
(
    interfaceToInterfaceMapping,
    directMapInterfaceToInterfaceMapping,
    dictionary
);


const scalar directMapInterfaceToInterfaceMapping::relTol_ = 0.01;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void directMapInterfaceToInterfaceMapping::checkZoneSizes() const
{
    if
    (
        zoneA().size() != zoneB().size()
     || zoneA().nPoints() != zoneB().nPoints()
    )
    {
        FatalErrorIn(type() + "::checkZoneSizes() const")
            << "ZoneA and zoneB interfaces are not conformal (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")"
            << nl
            << "directMap method requires conformal interfaces!"
            << abort(FatalError);
    }
}


void directMapInterfaceToInterfaceMapping::calcZoneAToZoneBFaceMap() const
{
    if (zoneAToZoneBFaceMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcZoneAToZoneBFaceMap() const")
            << "List already set!" << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Map name
    const word mapName
    (
        "zoneAToZoneBFaceMap_" + globalPatchA().patchName()
      + "_to_" + globalPatchB().patchName()
    );

    // Check if map needs to be read from disk
    IOobject mapHeader
    (
        mapName,
        globalPatchA().mesh().time().timeName(),
        globalPatchA().mesh(),
        IOobject::MUST_READ
    );

#ifdef OPENFOAM_NOT_EXTEND
    if (mapHeader.typeHeaderOk<labelIOList>(true))
#else
    if (mapHeader.headerOk())
#endif
    {
        // Read map
        Info<< "Reading " << mapName << " from disk" << endl;

        zoneAToZoneBFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchA().mesh().time().timeName(),
                    globalPatchA().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );
    }
    else
    {
        // Initialise map
        zoneAToZoneBFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchA().mesh().time().timeName(),
                    globalPatchA().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(zoneB().size(), -1)
            )
        );
        labelList& zoneAToZoneBMap = zoneAToZoneBFaceMapPtr_();

        // Perform N^2 search for corresponding faces
        // We will take 0.1% of the minEdgeLength as the exact match
        // relative tolerance

        const vectorField zoneACf(zoneA().faceCentres());
        const vectorField zoneBCf(zoneB().faceCentres());
        const scalar tol = relTol_*gMin(minEdgeLengths());

        forAll(zoneBCf, zoneBFaceI)
        {
            const vector& curZoneBCf = zoneBCf[zoneBFaceI];
            forAll(zoneACf, zoneAFaceI)
            {
                if (mag(curZoneBCf - zoneACf[zoneAFaceI]) < tol)
                {
                    zoneAToZoneBMap[zoneBFaceI] = zoneAFaceI;
                    break;
                }
            }
        }
    }

    if (gMin(zoneAToZoneBFaceMapPtr_()) == -1)
    {
        labelList& zoneAToZoneBMap = zoneAToZoneBFaceMapPtr_();
        const vectorField zoneBCf(zoneB().faceCentres());
        Info<< "A direct map was not found for the following faces:" << endl;
        forAll(zoneBCf, zoneBFaceI)
        {
            if (zoneAToZoneBMap[zoneBFaceI] == -1)
            {
                Info<< "face " << zoneBFaceI << ": " << zoneBCf[zoneBFaceI]
                    << endl;
            }
        }

        FatalErrorIn(type() + "::calcZoneAToZoneBFaceMap() const")
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB interfaces are not conformal (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")"
            << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


void directMapInterfaceToInterfaceMapping::calcZoneBToZoneAFaceMap() const
{
    if (zoneBToZoneAFaceMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcZoneBToZoneAFaceMap() const")
            << "Map already set!"
            << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Map name
    const word mapName
    (
        "zoneBToZoneAFaceMap_" + globalPatchB().patchName()
      + "_to_" + globalPatchA().patchName()
    );

    // Check if map needs to be read from disk
    IOobject mapHeader
    (
        mapName,
        globalPatchB().mesh().time().timeName(),
        globalPatchB().mesh(),
        IOobject::MUST_READ
    );

#ifdef OPENFOAM_NOT_EXTEND
    if (mapHeader.typeHeaderOk<labelIOList>(true))
#else
    if (mapHeader.headerOk())
#endif
    {
        // Read map
        Info<< "Reading " << mapName << "from disk" << endl;

        zoneBToZoneAFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchB().mesh().time().timeName(),
                    globalPatchB().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );
    }
    else
    {
        // Initialise map
        zoneBToZoneAFaceMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchB().mesh().time().timeName(),
                    globalPatchB().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(zoneA().size(), -1)
            )
        );
        labelList& zoneBToZoneAMap = zoneBToZoneAFaceMapPtr_();

        // Perform N^2 search for corresponding faces
        // We will take 0.1% of the minEdgeLength as the exact match
        // tolerance

        const vectorField zoneACf(zoneA().faceCentres());
        const vectorField zoneBCf(zoneB().faceCentres());
        const scalar tol = relTol_*gMin(minEdgeLengths());

        forAll(zoneACf, zoneAFaceI)
        {
            const vector& curZoneACf = zoneACf[zoneAFaceI];

            forAll(zoneBCf, zoneBFaceI)
            {
                if (mag(curZoneACf - zoneBCf[zoneBFaceI]) < tol)
                {
                    zoneBToZoneAMap[zoneAFaceI] = zoneBFaceI;
                    break;
                }
            }
        }
    }

    if (gMin(zoneBToZoneAFaceMapPtr_()) == -1)
    {
        labelList& zoneBToZoneAMap = zoneBToZoneAFaceMapPtr_();
        const vectorField zoneACf(zoneA().faceCentres());
        Info<< "A direct map was not found for the following faces:" << endl;
        forAll(zoneACf, zoneAFaceI)
        {
            if (zoneBToZoneAMap[zoneAFaceI] == -1)
            {
                Info<< "face " << zoneAFaceI << ": " << zoneACf[zoneAFaceI]
                    << endl;
            }
        }

        FatalErrorIn(type() + "::calcZoneBToZoneAFaceMap() const")
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB interfaces are not conformal (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")"
            << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


void directMapInterfaceToInterfaceMapping::calcZoneAToZoneBPointMap() const
{
    if (zoneAToZoneBPointMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcZoneAToZoneBPointMap() const")
            << "Map already set!"
            << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Map name
    const word mapName
    (
        "zoneAToZoneBPointMap_" + globalPatchA().patchName()
      + "_to_" + globalPatchB().patchName()
    );

    // Check if map needs to be read from disk
    IOobject mapHeader
    (
        mapName,
        globalPatchA().mesh().time().timeName(),
        globalPatchA().mesh(),
        IOobject::MUST_READ
    );

#ifdef OPENFOAM_NOT_EXTEND
    if (mapHeader.typeHeaderOk<labelIOList>(true))
#else
    if (mapHeader.headerOk())
#endif
    {
        // Read map
        Info<< "Reading " << mapName << " from disk" << endl;

        zoneAToZoneBPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchA().mesh().time().timeName(),
                    globalPatchA().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );
    }
    else
    {
        // Initialise map
        zoneAToZoneBPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchA().mesh().time().timeName(),
                    globalPatchA().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(zoneB().nPoints(), -1)
            )
        );
        labelList& zoneAToZoneBMap = zoneAToZoneBPointMapPtr_();

        // Perform N^2 search for corresponding points
        // We will take 0.1% of the minEdgeLength as the exact match
        // tolerance

        const vectorField& zoneALP = zoneA().localPoints();
        const vectorField& zoneBLP = zoneB().localPoints();
        const scalar tol = relTol_*gMin(minEdgeLengths());

        forAll(zoneBLP, zoneBPointI)
        {
            const vector& curZoneBLP = zoneBLP[zoneBPointI];

            forAll(zoneALP, zoneAPointI)
            {
                if (mag(curZoneBLP - zoneALP[zoneAPointI]) < tol)
                {
                    zoneAToZoneBMap[zoneBPointI] = zoneAPointI;
                    break;
                }
            }
        }
    }

    if (gMin(zoneAToZoneBPointMapPtr_()) == -1)
    {
        FatalErrorIn(type() + "::calcZoneAToZoneBPointMap() const")
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB interfaces are not conformal (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")"
            << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


void directMapInterfaceToInterfaceMapping::calcZoneBToZoneAPointMap() const
{
    if (zoneBToZoneAPointMapPtr_.valid())
    {
        FatalErrorIn(type() + "::calcZoneBToZoneAPointMap() const")
            << "Map already set!"
            << abort(FatalError);
    }

    // Check zones are conformal
    checkZoneSizes();

    // Map name
    const word mapName
    (
        "zoneBToZoneAPointMap_" + globalPatchB().patchName()
      + "_to_" + globalPatchA().patchName()
    );

    // Check if map needs to be read from disk
    IOobject mapHeader
    (
        mapName,
        globalPatchB().mesh().time().timeName(),
        globalPatchB().mesh(),
        IOobject::MUST_READ
    );

#ifdef OPENFOAM_NOT_EXTEND
    if (mapHeader.typeHeaderOk<labelIOList>(true))
#else
    if (mapHeader.headerOk())
#endif
    {
        // Read map
        Info<< "Reading " << mapName << "from disk" << endl;

        zoneBToZoneAPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchB().mesh().time().timeName(),
                    globalPatchB().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );
    }
    else
    {
        // Initialise map
        zoneBToZoneAPointMapPtr_.set
        (
            new labelIOList
            (
                IOobject
                (
                    mapName,
                    globalPatchB().mesh().time().timeName(),
                    globalPatchB().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                labelList(zoneA().nPoints(), -1)
            )
        );
        labelList& zoneBToZoneAMap = zoneBToZoneAPointMapPtr_();

        // Perform N^2 search for corresponding points
        // We will take 0.1% of the minEdgeLength as the exact match
        // tolerance

        const vectorField& zoneALP = zoneA().localPoints();
        const vectorField& zoneBLP = zoneB().localPoints();
        const scalar tol = relTol_*gMin(minEdgeLengths());

        forAll(zoneALP, zoneAPointI)
        {
            const vector& curZoneALP = zoneALP[zoneAPointI];

            forAll(zoneBLP, zoneBPointI)
            {
                if (mag(curZoneALP - zoneBLP[zoneBPointI]) < tol)
                {
                    zoneBToZoneAMap[zoneAPointI] = zoneBPointI;
                    break;
                }
            }
        }
    }

    if (gMin(zoneBToZoneAPointMapPtr_()) == -1)
    {
        FatalErrorIn(type() + "::calcZoneBToZoneAPointMap() const")
            << "Cannot calculate the map between interfaces!" << nl
            << "ZoneA and zoneB interfaces are not conformal (zoneA patch ="
            << " " << globalPatchA().patchName() << ", zoneB "
            << "patch = " << globalPatchB().patchName() << ")"
            << nl
            << "Direct mapping can only be used with conformal interfaces!"
            << abort(FatalError);
    }
}


void directMapInterfaceToInterfaceMapping::calcMinEdgeLengths() const
{
    if (minEdgeLengthsPtr_.valid())
    {
        FatalErrorIn
        (
            "void zoneAZoneBInterface::calcMinEdgeLengths() const"
        )   << "Pointer already set!" << abort(FatalError);
    }

    minEdgeLengthsPtr_.set(new scalarField(zoneA().nPoints(), 0));
    scalarField& minEdgeLength = minEdgeLengthsPtr_();

    const edgeList& edges = zoneA().edges();
    const vectorField& points = zoneA().localPoints();
    const labelListList& pointEdges = zoneA().pointEdges();

    forAll(points, pointI)
    {
        const labelList& curPointEdges = pointEdges[pointI];

        scalar minLength = GREAT;

        forAll(curPointEdges, edgeI)
        {
            const edge& curEdge = edges[curPointEdges[edgeI]];

            scalar Le = curEdge.mag(points);

            if (Le < minLength)
            {
                minLength = Le;
            }
        }

        minEdgeLength[pointI] = minLength;
    }
}


const scalarField& directMapInterfaceToInterfaceMapping::minEdgeLengths() const
{
    if (minEdgeLengthsPtr_.empty())
    {
        calcMinEdgeLengths();
    }

    return minEdgeLengthsPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

directMapInterfaceToInterfaceMapping::directMapInterfaceToInterfaceMapping
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
    zoneAToZoneBFaceMapPtr_(),
    zoneBToZoneAFaceMapPtr_(),
    zoneAToZoneBPointMapPtr_(),
    zoneBToZoneAPointMapPtr_(),
    minEdgeLengthsPtr_()
{
    // Force calculation of maps upon constuction as they may need to be read
    // from disk to allow restarts
    zoneAToZoneBFaceMap();
    zoneBToZoneAFaceMap();
    zoneAToZoneBPointMap();
    zoneBToZoneAPointMap();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelIOList&
directMapInterfaceToInterfaceMapping::zoneBToZoneAFaceMap() const
{
    if (zoneBToZoneAFaceMapPtr_.empty())
    {
        calcZoneBToZoneAFaceMap();
    }

    return zoneBToZoneAFaceMapPtr_;
}


const labelIOList&
directMapInterfaceToInterfaceMapping::zoneAToZoneBFaceMap() const
{
    if (zoneAToZoneBFaceMapPtr_.empty())
    {
        calcZoneAToZoneBFaceMap();
    }

    return zoneAToZoneBFaceMapPtr_;
}


const labelIOList&
directMapInterfaceToInterfaceMapping::zoneAToZoneBPointMap() const
{
    if (zoneAToZoneBPointMapPtr_.empty())
    {
        calcZoneAToZoneBPointMap();
    }

    return zoneAToZoneBPointMapPtr_;
}


const labelIOList&
directMapInterfaceToInterfaceMapping::zoneBToZoneAPointMap() const
{
    if (zoneBToZoneAPointMapPtr_.empty())
    {
        calcZoneBToZoneAPointMap();
    }

    return zoneBToZoneAPointMapPtr_();

}


void directMapInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFacesZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void directMapInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPointsZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void directMapInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFacesZoneToZone<vector>(fromZone, toZone, fromField, toField);
}


void directMapInterfaceToInterfaceMapping::transferPointsZoneToZone
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
