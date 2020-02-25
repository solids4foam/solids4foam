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

#ifdef OPENFOAMESIORFOUNDATION

#include "amiInterfaceToInterfaceMapping.H"
#include "addToRunTimeSelectionTable.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(amiInterfaceToInterfaceMapping, 0);
addToRunTimeSelectionTable
(
    interfaceToInterfaceMapping, amiInterfaceToInterfaceMapping, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void amiInterfaceToInterfaceMapping::makeInterpolator() const
{
    if (interpolatorPtr_.valid())
    {
        FatalErrorIn
        (
            "void amiInterfaceToInterfaceMapping::makeInterpolators() const"
        )   << "Pointer is already set!"
            << abort(FatalError);
    }

    Info<< "Create AMI zone-to-zone interpolator" << endl;

    interpolatorPtr_.set
    (
#ifdef OPENFOAMESI
        new AMIInterpolation<standAlonePatch, standAlonePatch>
#else
        new newAMIInterpolation<standAlonePatch, standAlonePatch>
#endif
        (
            zoneA(),
            zoneB(),
            faceAreaIntersect::tmMesh, // triMode
            true,   // requireMatch
            amiZoneInterpolation::imFaceAreaWeight, // interpolationMethodNames
            -1,     // lowWeightCorrection
            false,  // reverseTarget
            true    // use globalPolyPatch
        )
    );

    checkZoneAToZoneBError();
    //checkZoneBToZoneAError();
}


#ifdef OPENFOAMESI
const AMIInterpolation<standAlonePatch, standAlonePatch>&
#else
const newAMIInterpolation<standAlonePatch, standAlonePatch>&
#endif
amiInterfaceToInterfaceMapping::interpolator() const
{
    if (interpolatorPtr_.empty())
    {
        makeInterpolator();
    }

    return interpolatorPtr_();
}


void amiInterfaceToInterfaceMapping::checkZoneAToZoneBError() const
{
    // Reference to patch face centres
    const vectorField& patchAFaceCentres = patchA().faceCentres();

    // Construct global zone field
    const vectorField zoneAFaceCentres =
        globalPatchA().patchFaceToGlobal(patchAFaceCentres);

    // Interpolate global zone field from A to B
    const vectorField zoneBFaceCentres =
        interpolator().interpolateToTarget(zoneAFaceCentres);

    // Extract local patch field
    const vectorField patchBFaceCentres =
        globalPatchB().globalFaceToPatch(zoneBFaceCentres);

    // Print maximum error
    Info<< "interface-to-interface face error: "
        << gMax(mag(patchBFaceCentres - patchB().faceCentres()))
        << endl;

}


// void amiInterfaceToInterfaceMapping::checkZoneBToZoneAError() const
// {
//     const vectorField zoneBPointsAtFluid =
//         interpolator().interpolateToSourcePointInterpolate(zoneB().localPoints());

//     const vectorField& zoneAPoints = zoneA().localPoints();

//     Info<< "interface-to-interface point error: "
//         << gMax(mag(zoneAPoints - zoneBPointsAtFluid))
//         << endl;
// }


void amiInterfaceToInterfaceMapping::calcZoneAPointAddressing() const
{
    if (zoneAPointAddressingPtr_)
    {
        FatalErrorIn(type())
            << "zoneA points addressing already exists"
            << abort(FatalError);
    }

    zoneAPointAddressingPtr_ =
        new List<labelPair>
        (
            zoneA().nPoints(), labelPair(-1,-1)
        );
    List<labelPair>& zoneAPointAddr = *zoneAPointAddressingPtr_;

    // zoneAPointDistancePtr_ =
    //     new scalarField
    //     (
    //         zoneA().nPoints(), GREAT
    //     );
    // scalarField& zoneAPointDist = *zoneAPointDistancePtr_;

    const labelListList& zoneAFaceAddr = interpolator().srcAddress();
    const labelListList& pointFaces = zoneA().pointFaces();
    const faceList& zoneBFaces = zoneB().localFaces();
    const pointField& zoneBPoints = zoneB().localPoints();
    const pointField& zoneAPoints = zoneA().localPoints();

    forAll(zoneAPointAddr, pointI)
    {
        const point& P = zoneAPoints[pointI];
        labelHashSet possibleZoneBFacesSet;
        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            const label curFace = curPointFaces[faceI];
            const labelList& curZoneBFaces = zoneAFaceAddr[curFace];

            forAll(curZoneBFaces, fI)
            {
                if (!possibleZoneBFacesSet.found(curZoneBFaces[fI]))
                {
                    possibleZoneBFacesSet.insert(curZoneBFaces[fI]);
                }
            }
        }

        const labelList possibleZoneBFaces = possibleZoneBFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        //scalar distance = GREAT;

        forAll(possibleZoneBFaces, faceI)
        {
            const label curZoneBFace = possibleZoneBFaces[faceI];
            const face& f = zoneBFaces[curZoneBFace];
            const point ctr = Foam::average(f.points(zoneBPoints));
            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = zoneBPoints[f.nextLabel(pI)];

                const triPointRef t
                (
                    zoneBPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                const scalar A = mag(n);
                n /= A;

                // Intersection point
                const point I = P + n*(n & (t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                const scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curZoneBFace;
                    faceTriangle.second() = pI;

                    //distance = ((P - I) & n);
                }
            }
        }

        zoneAPointAddr[pointI] = faceTriangle;
        //zoneAPointDist[pointI] = distance;
    }

    // Info<< "zoneA point distance, max: "
    //     << max(zoneAPointDist)
    //     << ", avg: " << average(zoneAPointDist)
    //     << ", min: " << min(zoneAPointDist) << endl;


    // Check orientation

    const pointField& zoneAPointNormals = zoneA().pointNormals();

    const vectorField& zoneBFaceNormals = zoneB().faceNormals();

    scalarField orientation(zoneAPointAddr.size(), 0);

    label nIncorrectPoints = 0;

    forAll(zoneAPointAddr, pointI)
    {
        orientation[pointI] =
            (
                zoneAPointNormals[pointI]
              & zoneBFaceNormals[zoneAPointAddr[pointI].first()]
            );

        if (orientation[pointI] > -SMALL)
        {
            nIncorrectPoints++;
        }
    }

    Info<< "zoneA point orientation (< 0), max: "
        << max(orientation)
        << ", min: " << min(orientation) << ", nIncorrectPoints: "
        << nIncorrectPoints << "/" << zoneAPointAddr.size() << endl;
}


void amiInterfaceToInterfaceMapping::calcZoneAPointWeights() const
{
    if (zoneAPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void amiInterfaceToInterfaceMapping::calcZoneAPointWeights() const"
        )   << "pointer already set" << abort(FatalError);
    }

    // zoneAPointWeightsPtr_ =
    //     new FieldField<Field, scalar>(zoneA().nPoints());
    // FieldField<Field, scalar>& zoneAPointWeights = *zoneAPointWeightsPtr_;
    zoneAPointWeightsPtr_ =
        new FieldField<Field, scalar>(zoneA().nPoints());
    FieldField<Field, scalar>& zoneAPointWeights = *zoneAPointWeightsPtr_;

    const faceList& zoneBFaces = zoneB().localFaces();
    const pointField& zoneBPoints = zoneB().localPoints();
    const pointField& zoneAPoints = zoneA().localPoints();

    const List<labelPair>& addr = this->zoneAPointAddr();

    forAll(zoneAPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = zoneAPoints[pointI];
            const face& hitFace = zoneBFaces[addr[pointI].first()];
            const point ctr = Foam::average(hitFace.points(zoneBPoints));
            const label pI = addr[pointI].second();

            const triPointRef t
            (
                zoneBPoints[hitFace[pI]],
                zoneBPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            const point I = P + n*(n&(t.a() - P));

#ifdef OPENFOAMESI
            zoneAPointWeights[pointI] = List<scalar>(3);
#else
            zoneAPointWeights.set(pointI, new scalarField(3));
#endif

            // zoneAPointWeights[pointI][0] = t.Ni(0, I);
            // zoneAPointWeights[pointI][1] = t.Ni(1, I);
            // zoneAPointWeights[pointI][2] = t.Ni(2, I);

            zoneAPointWeights[pointI][0] = t.pointToBarycentric(I).a();
            zoneAPointWeights[pointI][1] = t.pointToBarycentric(I).b();
            zoneAPointWeights[pointI][2] = t.pointToBarycentric(I).c();
        }
        else
        {
#ifdef OPENFOAMESI
            zoneAPointWeights[pointI] = List<scalar>(0);
#else
            zoneAPointWeights.set(pointI, new scalarField(0));
#endif
        }
    }
}


void amiInterfaceToInterfaceMapping::calcZoneBPointAddressing() const
{
    if (zoneBPointAddressingPtr_)
    {
        FatalErrorIn(type())
            << "zoneB points addressing already exists"
            << abort(FatalError);
    }

    zoneBPointAddressingPtr_ =
        new List<labelPair>
        (
            zoneB().nPoints(), labelPair(-1,-1)
        );
    List<labelPair>& zoneBPointAddr = *zoneBPointAddressingPtr_;

    // zoneBPointDistancePtr_ =
    //     new scalarField
    //     (
    //         zoneB().nPoints(), GREAT
    //     );
    // scalarField& zoneBPointDist = *zoneBPointDistancePtr_;

    const labelListList& zoneBFaceAddr = interpolator().tgtAddress();
    const labelListList& pointFaces = zoneB().pointFaces();
    const faceList& zoneAFaces = zoneA().localFaces();
    const pointField& zoneAPoints = zoneA().localPoints();
    const pointField& zoneBPoints = zoneB().localPoints();

    forAll(zoneBPointAddr, pointI)
    {
        const point& P = zoneBPoints[pointI];
        labelHashSet possibleZoneAFacesSet;
        const labelList& curPointFaces = pointFaces[pointI];
        forAll(curPointFaces, faceI)
        {
            const label curFace = curPointFaces[faceI];
            const labelList& curZoneAFaces = zoneBFaceAddr[curFace];

            forAll(curZoneAFaces, fI)
            {
                if (!possibleZoneAFacesSet.found(curZoneAFaces[fI]))
                {
                    possibleZoneAFacesSet.insert(curZoneAFaces[fI]);
                }
            }
        }

        const labelList possibleZoneAFaces = possibleZoneAFacesSet.toc();

        scalar MinEta = -GREAT;
        labelPair faceTriangle(-1, -1);
        //scalar distance = GREAT;

        forAll(possibleZoneAFaces, faceI)
        {
            const label curZoneAFace = possibleZoneAFaces[faceI];
            const face& f = zoneAFaces[curZoneAFace];
            const point ctr = Foam::average(f.points(zoneAPoints));
            point nextPoint = ctr;

            for (label pI = 0; pI < f.size(); pI++)
            {
                nextPoint = zoneAPoints[f.nextLabel(pI)];

                const triPointRef t
                (
                    zoneAPoints[f[pI]],
                    nextPoint,
                    ctr
                );

                vector n = t.normal();
                const scalar A = mag(n);
                n /= A;

                // Intersection point
                const point I = P + n*(n & (t.a() - P));

                // Areal coordinates
                scalarField eta(3, 0);

                eta[0] = (triPointRef(I, t.b(), t.c()).normal() & n)/A;
                eta[1] = (triPointRef(I, t.c(), t.a()).normal() & n)/A;
                eta[2] = (triPointRef(I, t.a(), t.b()).normal() & n)/A;

                const scalar minEta = min(eta);

                if (minEta > MinEta)
                {
                    MinEta = minEta;
                    faceTriangle.first() = curZoneAFace;
                    faceTriangle.second() = pI;

                    //distance = ((P - I) & n);
                }
            }
        }

        zoneBPointAddr[pointI] = faceTriangle;
        //zoneBPointDist[pointI] = distance;
    }

    // Info<< "zoneB point distance, max: "
    //     << max(zoneBPointDist)
    //     << ", avg: " << average(zoneBPointDist)
    //     << ", min: " << min(zoneBPointDist) << endl;


    // Check orientation

    const pointField& zoneBPointNormals = zoneB().pointNormals();

    const vectorField& zoneAFaceNormals = zoneA().faceNormals();

    scalarField orientation(zoneBPointAddr.size(), 0);

    label nIncorrectPoints = 0;

    forAll(zoneBPointAddr, pointI)
    {
        orientation[pointI] =
            (
                zoneBPointNormals[pointI]
              & zoneAFaceNormals[zoneBPointAddr[pointI].first()]
            );

        if (orientation[pointI] > -SMALL)
        {
            nIncorrectPoints++;
        }
    }

    Info<< "zoneB point orientation (< 0), max: "
        << max(orientation)
        << ", min: " << min(orientation) << ", nIncorrectPoints: "
        << nIncorrectPoints << "/" << zoneBPointAddr.size() << endl;
}


void amiInterfaceToInterfaceMapping::calcZoneBPointWeights() const
{
    if (zoneBPointWeightsPtr_)
    {
        FatalErrorIn
        (
            "void amiInterfaceToInterfaceMapping::calcZoneBPointWeights() const"
        )   << "pointer already set" << abort(FatalError);
    }

    // zoneBPointWeightsPtr_ =
    //     new FieldField<Field, scalar>(zoneB().nPoints());
    // FieldField<Field, scalar>& zoneBPointWeights = *zoneBPointWeightsPtr_;
    zoneBPointWeightsPtr_ =
        new FieldField<Field, scalar>(zoneB().nPoints());
    FieldField<Field, scalar>& zoneBPointWeights = *zoneBPointWeightsPtr_;

    const faceList& zoneAFaces = zoneA().localFaces();
    const pointField& zoneAPoints = zoneA().localPoints();
    const pointField& zoneBPoints = zoneB().localPoints();
    const List<labelPair>& addr = this->zoneBPointAddr();

    forAll(zoneBPointWeights, pointI)
    {
        if (addr[pointI].first() != -1)
        {
            const point& P = zoneBPoints[pointI];
            const face& hitFace = zoneAFaces[addr[pointI].first()];
            const point ctr = Foam::average(hitFace.points(zoneAPoints));
            const label pI = addr[pointI].second();

            const triPointRef t
            (
                zoneAPoints[hitFace[pI]],
                zoneAPoints[hitFace.nextLabel(pI)],
                ctr
            );

            vector n = t.normal();
            n /= mag(n);

            // Intersection point
            const point I = P + n*(n&(t.a() - P));

#ifdef OPENFOAMESI
            zoneBPointWeights[pointI] = List<scalar>(3);
#else
            zoneBPointWeights.set(pointI, new scalarField(3));
#endif

            // zoneBPointWeights[pointI][0] = t.Ni(0, I);
            // zoneBPointWeights[pointI][1] = t.Ni(1, I);
            // zoneBPointWeights[pointI][2] = t.Ni(2, I);

            zoneBPointWeights[pointI][0] = t.pointToBarycentric(I).a();
            zoneBPointWeights[pointI][1] = t.pointToBarycentric(I).b();
            zoneBPointWeights[pointI][2] = t.pointToBarycentric(I).c();
        }
        else
        {
#ifdef OPENFOAMESI
            zoneBPointWeights[pointI] = List<scalar>(0);
#else
            zoneBPointWeights.set(pointI, new scalarField(0));
#endif
        }
    }
}

const List<labelPair>&
amiInterfaceToInterfaceMapping::zoneAPointAddr() const
{
    if (!zoneAPointAddressingPtr_)
    {
        calcZoneAPointAddressing();
    }

    return *zoneAPointAddressingPtr_;
}


const FieldField<Field, scalar>&
amiInterfaceToInterfaceMapping::zoneAPointWeights() const
{
    if (!zoneAPointWeightsPtr_)
    {
        calcZoneAPointWeights();
    }

    return *zoneAPointWeightsPtr_;
}


const List<labelPair>&
amiInterfaceToInterfaceMapping::zoneBPointAddr() const
{
    if (!zoneBPointAddressingPtr_)
    {
        calcZoneBPointAddressing();
    }

    return *zoneBPointAddressingPtr_;
}


const FieldField<Field, scalar>&
amiInterfaceToInterfaceMapping::zoneBPointWeights() const
{
    if (!zoneBPointWeightsPtr_)
    {
        calcZoneBPointWeights();
    }

    return *zoneBPointWeightsPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

amiInterfaceToInterfaceMapping::amiInterfaceToInterfaceMapping
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
    interpolatorPtr_(nullptr),
    zoneAPointAddressingPtr_(nullptr),
    zoneAPointWeightsPtr_(nullptr),
    zoneBPointAddressingPtr_(nullptr),
    zoneBPointWeightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void amiInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFacesZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void amiInterfaceToInterfaceMapping::transferPointsZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPointsZoneToZone<scalar>(fromZone, toZone, fromField, toField);
}


void amiInterfaceToInterfaceMapping::transferFacesZoneToZone
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFacesZoneToZone<vector>(fromZone, toZone, fromField, toField);
}


void amiInterfaceToInterfaceMapping::transferPointsZoneToZone
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

#endif // end of #ifdef OPENFOAMESIORFOUNDATION

// ************************************************************************* //
