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

#include "fsiImmersedInterface.H"
#include "HashTable.H"
#include "EdgeMap.H"
#include "DynamicList.H"
#include "compatibilityFunctions.H"
#ifdef OPENFOAM_ORG
    #include "polygonTriangulate.H"
#endif

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Append the triangles of face f, with the orientation of f, to tris and
    // return their number
    static label appendTriangles
    (
        const face& f,
        const pointField& points,
        DynamicList<face>& tris
    )
    {
#ifdef OPENFOAM_ORG
        // OpenFOAM.org has no face::triangles
        polygonTriangulate triEngine;
        triEngine.triangulate(UIndirectList<point>(points, f));
        const List<triFace> triPoints(triEngine.triPoints(f));

        for (const triFace& t : triPoints)
        {
            face tri(3);
            tri[0] = t[0];
            tri[1] = t[1];
            tri[2] = t[2];
            tris.append(tri);
        }

        return triPoints.size();
#else
        faceList faceTris(f.nTriangles(points));
        label triI = 0;
        f.triangles(points, triI, faceTris);

        for (const face& tri : faceTris)
        {
            tris.append(tri);
        }

        return faceTris.size();
#endif
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fsiImmersedInterface::buildSurface
(
    const standAlonePatch& zone,
    const pointField& zonePoints,
    const List<const standAlonePatch*>& closureZones,
    const List<pointField>& closurePoints,
    const Vector<label>& fluidSolutionD,
    const boundBox& fluidBounds
)
{
    if (zonePoints.size() != zone.nPoints())
    {
        FatalErrorInFunction
            << "The interface zone has " << zone.nPoints() << " points but "
            << zonePoints.size() << " reference points are given"
            << abort(FatalError);
    }

    // Merge the points of the zones by their coordinates: the zones are
    // patches of the same mesh, so the points they share have the same
    // coordinates up to round-off (in parallel, a shared point is the
    // average of its copies on the processors, whose number may differ
    // between the zones), and are merged within a tolerance that is small
    // relative to the extent of the interface. The interface zone is first,
    // so that a shared point follows the interface
    const scalar mergeTol = 1e-8*mag(boundBox(zonePoints, false).span());

    // Bins of the merged points, of the width of the tolerance, keyed by a
    // hash of the bin indices (a collision only costs a distance check)
    const scalar binWidth = max(mergeTol, VSMALL);
    HashTable<DynamicList<label>, label, Hash<label>>
        pointBins(2*(zonePoints.size() + 1));

    auto binKey = [&](const point& p, const label di, const label dj,
        const label dk)
    {
        const long i = long(std::floor(p.x()/binWidth)) + di;
        const long j = long(std::floor(p.y()/binWidth)) + dj;
        const long k = long(std::floor(p.z()/binWidth)) + dk;
        return label((i*73856093L) ^ (j*19349663L) ^ (k*83492791L));
    };

    DynamicList<point> points(zonePoints.size());
    DynamicList<label> pointZones(zonePoints.size());
    DynamicList<label> pointIndices(zonePoints.size());

    // Merged point within the tolerance of p, or -1
    auto findPoint = [&](const point& p)
    {
        for (label di = -1; di <= 1; ++di)
        {
            for (label dj = -1; dj <= 1; ++dj)
            {
                for (label dk = -1; dk <= 1; ++dk)
                {
                    const auto iter = pointBins.find(binKey(p, di, dj, dk));

                    if (iter == pointBins.end())
                    {
                        continue;
                    }

                    for (const label pointi : iter())
                    {
                        if (mag(points[pointi] - p) <= mergeTol)
                        {
                            return pointi;
                        }
                    }
                }
            }
        }

        return label(-1);
    };

    // Map from the points of each zone to the surface points
    labelList zoneMap(zonePoints.size(), -1);
    List<labelList> closureMaps(closureZones.size());

    auto addPoints = [&]
    (
        const pointField& pts,
        const label zonei,
        labelList& map
    )
    {
        map.setSize(pts.size(), -1);
        forAll(pts, i)
        {
            const label found = findPoint(pts[i]);

            if (found >= 0)
            {
                map[i] = found;
            }
            else
            {
                map[i] = points.size();
                pointBins(binKey(pts[i], 0, 0, 0)).append(points.size());
                points.append(pts[i]);
                pointZones.append(zonei);
                pointIndices.append(i);
            }
        }
    };

    addPoints(zonePoints, -1, zoneMap);
    forAll(closureZones, zonei)
    {
        if (closurePoints[zonei].size() != closureZones[zonei]->nPoints())
        {
            FatalErrorInFunction
                << "Closure zone " << zonei << " has "
                << closureZones[zonei]->nPoints() << " points but "
                << closurePoints[zonei].size() << " reference points are "
                << "given" << abort(FatalError);
        }

        addPoints(closurePoints[zonei], zonei, closureMaps[zonei]);
    }

    // Triangulate the faces of the zones
    DynamicList<face> triangles(2*zone.size());
    DynamicList<label> triangleFaces(2*zone.size());
    const pointField mergedPoints(points);

    auto addFaces = [&]
    (
        const faceList& faces,
        const labelList& map,
        const bool interfaceZone
    )
    {
        forAll(faces, facei)
        {
            face f(faces[facei]);
            forAll(f, i)
            {
                f[i] = map[f[i]];
            }

            const label nTris =
                appendTriangles(f, mergedPoints, triangles);

            for (label i = 0; i < nTris; ++i)
            {
                triangleFaces.append(interfaceZone ? facei : -1);
            }
        }
    };

    addFaces(zone.localFaces(), zoneMap, true);
    forAll(closureZones, zonei)
    {
        addFaces(closureZones[zonei]->localFaces(), closureMaps[zonei], false);
    }

    pointOffsets_.setSize(points.size(), vector::zero);

    // Two-dimensional fluid mesh: extend the surface through the mesh in
    // the empty direction and cap its boundary loops
    label emptyDir = -1;
    label nSolutionD = 0;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (fluidSolutionD[d] == 1)
        {
            ++nSolutionD;
        }
        else
        {
            emptyDir = d;
        }
    }

    if (nSolutionD == 2)
    {
        const label e = emptyDir;
        const scalar zMid = 0.5*(fluidBounds.min()[e] + fluidBounds.max()[e]);
        const scalar span = fluidBounds.span()[e];

        // The points on either side of the mid-plane are moved a mesh
        // thickness beyond the bounds of the mesh
        label nBelow = 0;
        forAll(points, i)
        {
            const scalar z = points[i][e];
            const scalar target =
            (
                z < zMid
              ? fluidBounds.min()[e] - span
              : fluidBounds.max()[e] + span
            );
            pointOffsets_[i][e] = target - z;

            if (z < zMid)
            {
                ++nBelow;
            }
        }

        if (nBelow == 0 || nBelow == points.size())
        {
            FatalErrorInFunction
                << "The surface of immersed body " << bodyName_
                << " is on one side of the mid-plane of the fluid mesh in "
                << "its empty direction " << e << ": the solid mesh must "
                << "span the mid-plane" << abort(FatalError);
        }

        // Boundary edges: those of a single triangle
        EdgeMap<label> edgeCount(4*triangles.size());
        for (const face& tri : triangles)
        {
            for (label i = 0; i < 3; ++i)
            {
                const edge ed(tri[i], tri[(i + 1) % 3]);
                EdgeMap<label>::iterator iter = edgeCount.find(ed);

                if (iter != edgeCount.end())
                {
                    ++iter();
                }
                else
                {
                    edgeCount.insert(ed, 1);
                }
            }
        }

        // Boundary edges of each point
        List<DynamicList<label>> pointEdges(points.size());
        DynamicList<edge> boundaryEdges(edgeCount.size());
        forAllConstIter(EdgeMap<label>, edgeCount, iter)
        {
            if (iter() == 1)
            {
                pointEdges[iter.key()[0]].append(boundaryEdges.size());
                pointEdges[iter.key()[1]].append(boundaryEdges.size());
                boundaryEdges.append(iter.key());
            }
        }

        forAll(pointEdges, pointi)
        {
            const label nEdges = pointEdges[pointi].size();

            if (nEdges != 0 && nEdges != 2)
            {
                FatalErrorInFunction
                    << "The boundary of the surface of immersed body "
                    << bodyName_ << " is not a set of simple loops: point "
                    << points[pointi] << " has " << nEdges
                    << " boundary edges" << abort(FatalError);
            }
        }

        // Walk the boundary loops
        boolList edgeVisited(boundaryEdges.size(), false);
        label nCaps = 0;

        forAll(boundaryEdges, starti)
        {
            if (edgeVisited[starti])
            {
                continue;
            }

            DynamicList<label> loop(64);
            label edgei = starti;
            label pointi = boundaryEdges[starti][0];

            while (!edgeVisited[edgei])
            {
                edgeVisited[edgei] = true;
                loop.append(pointi);
                pointi = boundaryEdges[edgei].otherVertex(pointi);

                const DynamicList<label>& pe = pointEdges[pointi];
                edgei = (pe[0] == edgei ? pe[1] : pe[0]);
            }

            if (loop.size() < 3)
            {
                FatalErrorInFunction
                    << "A boundary loop of the surface of immersed body "
                    << bodyName_ << " has " << loop.size() << " points"
                    << abort(FatalError);
            }

            // The loop must lie on one side of the mid-plane
            const bool below = points[loop[0]][e] < zMid;
            for (const label pointi : loop)
            {
                if ((points[pointi][e] < zMid) != below)
                {
                    FatalErrorInFunction
                        << "The surface of immersed body " << bodyName_
                        << " is open: a boundary loop crosses the mid-plane "
                        << "of the fluid mesh in its empty direction. The "
                        << "union of the interface and closure patches must "
                        << "be closed apart from the patches of the empty "
                        << "direction" << abort(FatalError);
                }
            }

            // Cap the loop, with the normal pointing away from the mesh
            face polygon(loop);
            const pointField offsetPoints(pointField(points) + pointOffsets_);
#ifdef OPENFOAM_ORG
            const vector n(polygon.area(offsetPoints));
#else
            const vector n(polygon.normal(offsetPoints));
#endif

            if ((n[e] < 0) != below)
            {
                polygon = polygon.reverseFace();
            }

            const label nTris =
                appendTriangles(polygon, offsetPoints, triangles);

            for (label i = 0; i < nTris; ++i)
            {
                triangleFaces.append(-1);
            }

            ++nCaps;
        }

        Info<< "    Immersed interface " << bodyName_ << ": "
            << nCaps << " boundary loops capped in direction " << e << endl;
    }

    pointZones_.transfer(pointZones);
    pointIndices_.transfer(pointIndices);
    triangleFaces_.transfer(triangleFaces);

    Info<< "    Immersed interface " << bodyName_ << ": " << points.size()
        << " points and " << triangleFaces_.size() << " triangles from "
        << zone.size() << " interface faces and " << closureZones.size()
        << " closure patches" << endl;

    const pointField surfacePoints(pointField(points) + pointOffsets_);
    ib_.setReferenceSurface(bodyi_, surfacePoints, faceList(triangles));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fsiImmersedInterface::fsiImmersedInterface
(
    fsiImmersedBoundary& ib,
    const label bodyi,
    const word& bodyName,
    const standAlonePatch& zone,
    const pointField& zonePoints,
    const List<const standAlonePatch*>& closureZones,
    const List<pointField>& closurePoints,
    const Vector<label>& fluidSolutionD,
    const boundBox& fluidBounds
)
:
    ib_(ib),
    bodyi_(bodyi),
    bodyName_(bodyName),
    pointZones_(),
    pointIndices_(),
    pointOffsets_(),
    triangleFaces_(),
    nZoneFaces_(zone.size()),
    zonePoints0_(zonePoints),
    closurePoints0_(closurePoints),
    nFacesWithoutPoints_(0)
{
    buildSurface
    (
        zone,
        zonePoints,
        closureZones,
        closurePoints,
        fluidSolutionD,
        fluidBounds
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fsiImmersedInterface::setPoints0
(
    const pointField& zonePoints,
    const List<pointField>& closurePoints
)
{
    if
    (
        zonePoints.size() != zonePoints0_.size()
     || closurePoints.size() != closurePoints0_.size()
    )
    {
        FatalErrorInFunction
            << "The numbers of points have changed" << abort(FatalError);
    }

    zonePoints0_ = zonePoints;
    closurePoints0_ = closurePoints;
}


void Foam::fsiImmersedInterface::move
(
    const vectorField& zoneDisplacement,
    const vectorField& zoneVelocity
)
{
    if
    (
        zoneDisplacement.size() != zonePoints0_.size()
     || zoneVelocity.size() != zonePoints0_.size()
    )
    {
        FatalErrorInFunction
            << "The interface zone has " << zonePoints0_.size()
            << " points but the displacement has " << zoneDisplacement.size()
            << " and the velocity " << zoneVelocity.size()
            << abort(FatalError);
    }

    pointField points(pointZones_.size());
    vectorField velocities(pointZones_.size(), vector::zero);

    forAll(pointZones_, i)
    {
        const label zonei = pointZones_[i];
        const label pointi = pointIndices_[i];

        if (zonei < 0)
        {
            points[i] = zonePoints0_[pointi] + zoneDisplacement[pointi];
            velocities[i] = zoneVelocity[pointi];
        }
        else
        {
            points[i] = closurePoints0_[zonei][pointi];
        }

        points[i] += pointOffsets_[i];
    }

    ib_.moveSurface(bodyi_, points, velocities);
}


Foam::tmp<Foam::vectorField> Foam::fsiImmersedInterface::zoneTraction() const
{
    scalarField areas;
    const vectorField triangleTraction(ib_.surfaceTraction(bodyi_, areas));

    if (triangleTraction.size() != triangleFaces_.size())
    {
        FatalErrorInFunction
            << "The immersed boundary gives " << triangleTraction.size()
            << " tractions for " << triangleFaces_.size() << " triangles"
            << abort(FatalError);
    }

    tmp<vectorField> tzoneTraction(new vectorField(nZoneFaces_, vector::zero));
    vectorField& zoneTraction = tmpRef(tzoneTraction);
    scalarField zoneAreas(nZoneFaces_, 0.0);

    forAll(triangleFaces_, i)
    {
        const label facei = triangleFaces_[i];

        if (facei >= 0)
        {
            zoneTraction[facei] += areas[i]*triangleTraction[i];
            zoneAreas[facei] += areas[i];
        }
    }

    nFacesWithoutPoints_ = 0;
    forAll(zoneTraction, facei)
    {
        if (zoneAreas[facei] > VSMALL)
        {
            zoneTraction[facei] /= zoneAreas[facei];
        }
        else
        {
            zoneTraction[facei] = vector::zero;
            ++nFacesWithoutPoints_;
        }
    }

    return tzoneTraction;
}


// ************************************************************************* //
