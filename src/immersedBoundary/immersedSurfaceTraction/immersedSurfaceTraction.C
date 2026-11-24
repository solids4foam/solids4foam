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

#include "immersedSurfaceTraction.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "scalarMatrices.H"
#include "SVD.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

const Foam::Enum<Foam::immersedSurfaceTraction::quadratureType>
Foam::immersedSurfaceTraction::quadratureTypeNames
({
    { quadratureType::greedy, "greedy" },
    { quadratureType::none, "none" },
    { quadratureType::clustered, "clustered" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalarField Foam::immersedSurfaceTraction::cellWidths
(
    const pointField& pts
) const
{
    // polyMesh::findCell constructs the tet base points collectively
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();

    scalarField widths(pts.size(), Zero);

    forAll(pts, i)
    {
        const label celli = mesh_.findCell(pts[i]);

        if (celli >= 0)
        {
            widths[i] = immersedBody::cellWidth(mesh_, celli);
        }
    }

    // Native all-reduce of the contiguous values
    reduce
    (
        widths.data(),
        widths.size(),
        maxOp<scalar>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    return widths;
}


void Foam::immersedSurfaceTraction::createPoints()
{
    const triSurface& surf = body_.surface();
    const boundBox& meshBb = mesh_.bounds();

    // The candidates of the clustered quadrature are made in the reference
    // configuration of the surface, those of the others in the current one
    const bool clustered = (quadrature_ == quadratureType::clustered);
    const pointField& pts = clustered ? body_.points0() : surf.points();

    label nD = 0;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solD_[d] == 1)
        {
            ++nD;
        }
    }

    // Candidate points, kept if they are inside the mesh
    DynamicList<label> faces;
    DynamicList<FixedList<scalar, 3>> barycentric;
    DynamicList<scalar> areas;
    DynamicList<point> positions;

    if (nD == 2)
    {
        // Intersect the surface with the mid-plane of the empty direction
        direction e = 0;
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            if (solD_[d] != 1)
            {
                e = d;
            }
        }

        const scalar zMid = 0.5*(meshBb.min()[e] + meshBb.max()[e]);
        const scalar thickness = meshBb.span()[e];

        // Mean width of the cells, for the spacing of the points
        const scalar hMean =
            Foam::sqrt(gAverage(mesh_.V())/thickness);

        forAll(surf, facei)
        {
            const labelledTri& f = surf[facei];

            DynamicList<FixedList<scalar, 3>> ends(2);
            DynamicList<point> endPoints(2);

            for (label k = 0; k < 3; ++k)
            {
                const label k1 = (k + 1) % 3;
                const scalar dk = pts[f[k]][e] - zMid;
                const scalar dk1 = pts[f[k1]][e] - zMid;

                if ((dk < 0) != (dk1 < 0))
                {
                    const scalar t = dk/(dk - dk1);

                    FixedList<scalar, 3> b(Zero);
                    b[k] = 1 - t;
                    b[k1] = t;

                    ends.append(b);
                    endPoints.append((1 - t)*pts[f[k]] + t*pts[f[k1]]);
                }
            }

            if (ends.size() != 2)
            {
                continue;
            }

            const scalar L = mag(endPoints[1] - endPoints[0]);

            if (L < VSMALL)
            {
                continue;
            }

            // Segments much longer than the cells (e.g. a plane surface
            // extending beyond the mesh) are subdivided at the mean width
            const label n = max(label(1), label(std::ceil(L/hMean)));

            for (label k = 0; k < n; ++k)
            {
                const scalar t = (k + 0.5)/n;

                FixedList<scalar, 3> b;
                for (label j = 0; j < 3; ++j)
                {
                    b[j] = (1 - t)*ends[0][j] + t*ends[1][j];
                }

                faces.append(facei);
                barycentric.append(b);
                areas.append(L/n*thickness);
                positions.append
                (
                    b[0]*pts[f[0]] + b[1]*pts[f[1]] + b[2]*pts[f[2]]
                );
            }
        }
    }
    else
    {
        const scalar hMean = Foam::cbrt(gAverage(mesh_.V()));

        // Uniform subdivision of each face
        forAll(surf, facei)
        {
            const labelledTri& f = surf[facei];
            const scalar A = f.mag(pts);

            scalar maxEdge = 0;
            for (label k = 0; k < 3; ++k)
            {
                maxEdge =
                    max(maxEdge, mag(pts[f[(k + 1) % 3]] - pts[f[k]]));
            }

            const label n = max(label(1), label(std::ceil(maxEdge/hMean)));

            // Centroids of the n^2 sub-triangles
            auto add = [&](const scalar b0, const scalar b1)
            {
                FixedList<scalar, 3> b;
                b[0] = b0;
                b[1] = b1;
                b[2] = 1 - b0 - b1;
                faces.append(facei);
                barycentric.append(b);
                areas.append(A/scalar(n*n));
                positions.append
                (
                    b[0]*pts[f[0]] + b[1]*pts[f[1]] + b[2]*pts[f[2]]
                );
            };

            for (label i = 0; i < n; ++i)
            {
                for (label j = 0; i + j < n; ++j)
                {
                    add((i + 1.0/3.0)/n, (j + 1.0/3.0)/n);

                    if (i + j < n - 1)
                    {
                        add((i + 2.0/3.0)/n, (j + 2.0/3.0)/n);
                    }
                }
            }
        }
    }

    // Thin the points to a spacing of about 0.7 times the mean cell width,
    // for surfaces whose triangles are smaller than the cells: the points are
    // kept in order if no kept point is closer, and the area of each point
    // is given to the nearest kept point
    if (quadrature_ == quadratureType::greedy)
    {
        scalar thickness = 1;
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            if (solD_[d] != 1)
            {
                thickness = meshBb.span()[d];
            }
        }
        const scalar hMean =
        (
            nD == 2
          ? Foam::sqrt(gAverage(mesh_.V())/thickness)
          : Foam::cbrt(gAverage(mesh_.V()))
        );
        const scalar spacing = 0.7*hMean;

        auto key = [&](const point& x, const label di, const label dj,
            const label dk)
        {
            const int64_t i = int64_t(std::floor(x.x()/spacing)) + di;
            const int64_t j = int64_t(std::floor(x.y()/spacing)) + dj;
            const int64_t k = int64_t(std::floor(x.z()/spacing)) + dk;
            return (i*73856093) ^ (j*19349663) ^ (k*83492791);
        };

        // The points inside and outside the mesh are thinned separately, so
        // that the area of the surface in the mesh is kept
        const scalarField candidateWidths(cellWidths(pointField(positions)));
        auto inMesh = [&](const label i) { return candidateWidths[i] > 0; };

        // Kept points (indices into kept) by the cell of a grid of the
        // spacing
        HashTable<DynamicList<label>, int64_t> grid;
        DynamicList<label> kept(positions.size());

        // Nearest kept point within the spacing, on the same side of the
        // mesh boundary as the point i, -1 if none
        auto nearestKept = [&](const label pi, const scalar maxDist)
        {
            const point& x = positions[pi];
            label best = -1;
            scalar bestDistSqr = sqr(maxDist);
            for (label di = -1; di <= 1; ++di)
            {
                for (label dj = -1; dj <= 1; ++dj)
                {
                    for (label dk = -1; dk <= 1; ++dk)
                    {
                        const auto iter = grid.cfind(key(x, di, dj, dk));
                        if (!iter.good())
                        {
                            continue;
                        }
                        for (const label k : iter.val())
                        {
                            if (inMesh(kept[k]) != inMesh(pi))
                            {
                                continue;
                            }
                            const scalar dSqr = magSqr(positions[kept[k]] - x);
                            if (dSqr < bestDistSqr)
                            {
                                bestDistSqr = dSqr;
                                best = k;
                            }
                        }
                    }
                }
            }
            return best;
        };

        forAll(positions, i)
        {
            if (nearestKept(i, spacing) < 0)
            {
                grid(key(positions[i], 0, 0, 0)).append(kept.size());
                kept.append(i);
            }
        }

        scalarField keptAreas(kept.size(), Zero);
        forAll(positions, i)
        {
            const label k = nearestKept(i, 2*spacing);
            keptAreas[k >= 0 ? k : 0] += areas[i];
        }

        faces_.setSize(kept.size());
        barycentric_.setSize(kept.size());
        scalar meshArea = 0;
        forAll(kept, k)
        {
            faces_[k] = faces[kept[k]];
            barycentric_[k] = barycentric[kept[k]];
            if (inMesh(kept[k]))
            {
                meshArea += keptAreas[k];
            }
        }
        areas_.transfer(keptAreas);

        Info<< "    Immersed body " << body_.name() << ": area in the mesh "
            << meshArea << endl;
    }
    else if (clustered)
    {
        // Clusters of the candidates by the cell of a uniform grid of
        // spacing s, the dominant signed axis of the normal and whether they
        // are inside the mesh, numbered in the order of these keys, so that
        // they do not depend on the order of the triangles. The candidates
        // (positions), the cell widths from which s is set, the inside or
        // outside of the mesh, the grid and the choice of the points are
        // those of the reference configuration (pts), so that the clusters
        // move with the body and are a function of the reference surface
        // and the mesh only: they do not depend on the time at which they
        // are made (e.g. a restart), if the mesh is the same
        const scalarField candidateWidths(cellWidths(pointField(positions)));

        scalar hMin = GREAT;
        forAll(candidateWidths, i)
        {
            if (candidateWidths[i] > 0)
            {
                hMin = min(hMin, candidateWidths[i]);
            }
        }

        // Without a candidate in the mesh, the mean width of the cells
        scalar thickness = 1;
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            if (solD_[d] != 1)
            {
                thickness = meshBb.span()[d];
            }
        }
        const scalar hMean =
        (
            nD == 2
          ? Foam::sqrt(gAverage(mesh_.V())/thickness)
          : Foam::cbrt(gAverage(mesh_.V()))
        );
        const scalar s = spacingFactor_*(hMin < GREAT ? hMin : hMean);

        // The grid uses the solved directions only: in two dimensions every
        // candidate is on the mid-plane of the mesh, which may lie on a plane
        // of the grid, where rounding would split the candidates
        auto cell = [&](const point& x, const direction d)
        {
            return solD_[d] == 1 ? int64_t(std::floor(x[d]/s)) : int64_t(0);
        };

        typedef FixedList<int64_t, 5> clusterKey;
        auto keyOf = [&](const label i)
        {
            const point& x = positions[i];
            const vector nf(surf[faces[i]].areaNormal(pts));
            direction a = 0;
            for (direction d = 1; d < vector::nComponents; ++d)
            {
                if (mag(nf[d]) > (1 + 1e-12)*mag(nf[a]))
                {
                    a = d;
                }
            }

            clusterKey key;
            key[0] = cell(x, vector::X);
            key[1] = cell(x, vector::Y);
            key[2] = cell(x, vector::Z);
            key[3] = int64_t(2*a + (nf[a] < 0 ? 1 : 0));
            key[4] = int64_t(candidateWidths[i] > 0 ? 1 : 0);
            return key;
        };

        // Distances and positions are compared on a grid of 1e-9 s, so that
        // rounding differences cannot change the choice between equidistant
        // candidates or the order of the clusters
        const scalar q = 1e-9*s;
        auto quantise = [q](const scalar x)
        {
            return int64_t(std::llround(x/q));
        };

        // Groups of the candidates with the same key, numbered in the
        // (lexicographic) order of the keys
        labelList groupOf(positions.size());
        label nGroups = 0;
        {
            List<clusterKey> keys(positions.size());
            forAll(positions, i)
            {
                keys[i] = keyOf(i);
            }

            const labelList order(sortedOrder(keys));

            forAll(order, j)
            {
                if (j == 0 || keys[order[j]] != keys[order[j - 1]])
                {
                    ++nGroups;
                }
                groupOf[order[j]] = nGroups - 1;
            }
        }
        List<DynamicList<label>> groups(nGroups);
        forAll(groupOf, i)
        {
            groups[groupOf[i]].append(i);
        }

        // Each group is split into the connected components of the
        // triangles of its candidates (edge neighbours within the group), so
        // that a cluster never joins surface regions that are not connected
        // locally (the layers of a folded or the faces of nearby bodies),
        // ordered by their smallest candidate position
        const labelListList& faceFaces = surf.faceFaces();
        labelList root(surf.size(), -1);
        auto findRoot = [&root](label f)
        {
            while (root[f] != f)
            {
                root[f] = root[root[f]];
                f = root[f];
            }
            return f;
        };

        DynamicList<labelList> clusters(groups.size());
        for (const DynamicList<label>& g : groups)
        {
            for (const label i : g)
            {
                root[faces[i]] = faces[i];
            }
            for (const label i : g)
            {
                for (const label nbr : faceFaces[faces[i]])
                {
                    if (root[nbr] >= 0)
                    {
                        const label a = findRoot(faces[i]);
                        const label b = findRoot(nbr);
                        root[max(a, b)] = min(a, b);
                    }
                }
            }

            // The components by their root, each with its candidates in the
            // order of the group (the sort is stable)
            labelList roots(g.size());
            forAll(g, j)
            {
                roots[j] = findRoot(faces[g[j]]);
            }
            const labelList byRoot(sortedOrder(roots));

            DynamicList<labelList> parts;
            DynamicList<FixedList<int64_t, 4>> partKeys;
            for (label j0 = 0; j0 < byRoot.size(); )
            {
                label j1 = j0 + 1;
                while
                (
                    j1 < byRoot.size()
                 && roots[byRoot[j1]] == roots[byRoot[j0]]
                )
                {
                    ++j1;
                }

                // The smallest candidate position of the component, then its
                // root (which only separates components at the same
                // position, i.e. duplicated triangles)
                labelList part(j1 - j0);
                FixedList<int64_t, 4> partKey;
                for (label j = j0; j < j1; ++j)
                {
                    const label i = g[byRoot[j]];
                    part[j - j0] = i;

                    FixedList<int64_t, 4> x;
                    x[0] = quantise(positions[i].x());
                    x[1] = quantise(positions[i].y());
                    x[2] = quantise(positions[i].z());
                    x[3] = roots[byRoot[j0]];
                    if (j == j0 || x < partKey)
                    {
                        partKey = x;
                    }
                }
                parts.append(part);
                partKeys.append(partKey);

                j0 = j1;
            }

            for (const label k : sortedOrder(partKeys))
            {
                clusters.append(parts[k]);
            }

            for (const label i : g)
            {
                root[faces[i]] = -1;
            }
        }
        const label nClusters = clusters.size();
        members_.transfer(clusters);

        // The parent face of each point is that of its candidate nearest to
        // the area centroid of the cluster; the point itself is placed by
        // updatePoints

        faces_.setSize(nClusters);
        barycentric_.setSize(nClusters);
        scalarField clusterAreas(nClusters, Zero);
        scalar meshArea = 0;
        forAll(members_, k)
        {
            const labelList& m = members_[k];

            scalar A = 0;
            vector xc(Zero);
            for (const label i : m)
            {
                A += areas[i];
                xc += areas[i]*positions[i];
            }
            xc /= max(A, VSMALL);

            label best = -1;
            FixedList<int64_t, 4> bestKey(int64_t(0));
            for (const label i : m)
            {
                const point& x = positions[i];
                FixedList<int64_t, 4> key;
                key[0] = quantise(mag(x - xc));
                key[1] = quantise(x.x());
                key[2] = quantise(x.y());
                key[3] = quantise(x.z());
                if (best < 0 || key < bestKey)
                {
                    best = i;
                    bestKey = key;
                }
            }

            faces_[k] = faces[best];
            barycentric_[k] = barycentric[best];
            clusterAreas[k] = A;
            if (candidateWidths[best] > 0)
            {
                meshArea += A;
            }
        }

        candidateFaces_ = faces;
        candidateBarycentric_ = barycentric;
        candidateAreas0_ = areas;
        candidateAreas_ = areas;
        candidateFaceAreas0_.setSize(faces.size());
        forAll(faces, i)
        {
            candidateFaceAreas0_[i] = mag(surf[faces[i]].areaNormal(pts));
        }
        areas_.transfer(clusterAreas);

        Info<< "    Immersed body " << body_.name() << ": clustered "
            << positions.size() << " candidate points at a spacing of " << s
            << ", area in the mesh " << meshArea << endl;
    }
    else
    {
        faces_.transfer(faces);
        barycentric_.transfer(barycentric);
        areas_.transfer(areas);
    }

    Info<< "    Immersed body " << body_.name() << ": " << faces_.size()
        << " traction points, area " << sum(areas_) << endl;

    pressureScale_.setSize(faces_.size());
    pressureScale_ = 1;

    areas0_ = areas_;
    faceAreas0_.setSize(faces_.size());
    forAll(faces_, i)
    {
        faceAreas0_[i] = mag(surf[faces_[i]].areaNormal(surf.points()));
    }
}


void Foam::immersedSurfaceTraction::updatePoints()
{
    const triSurface& surf = body_.surface();
    const pointField& pts = surf.points();
    const vectorField& faceNormals = surf.faceNormals();

    points_.setSize(faces_.size());
    normals_.setSize(faces_.size());

    forAll(faces_, i)
    {
        const labelledTri& f = surf[faces_[i]];
        const FixedList<scalar, 3>& b = barycentric_[i];

        points_[i] = b[0]*pts[f[0]] + b[1]*pts[f[1]] + b[2]*pts[f[2]];
        normals_[i] = faceNormals[faces_[i]];
    }

    if (members_.size())
    {
        // Clustered: the area, the vector area and the area centroid of each
        // cluster from the current triangles of its candidates (whose areas
        // change with those of their triangles); the point is the centroid
        // projected on the current plane of its parent face, its normal the
        // direction of the vector area
        if (!body_.rigid())
        {
            forAll(candidateAreas_, i)
            {
                candidateAreas_[i] =
                    candidateAreas0_[i]
                   *mag(surf[candidateFaces_[i]].areaNormal(pts))
                   /max(candidateFaceAreas0_[i], VSMALL);
            }
        }

        forAll(faces_, k)
        {
            scalar A = 0;
            vector a(Zero);
            vector xc(Zero);
            for (const label i : members_[k])
            {
                const label facei = candidateFaces_[i];
                const labelledTri& f = surf[facei];
                const FixedList<scalar, 3>& b = candidateBarycentric_[i];
                const scalar Ai = candidateAreas_[i];
                A += Ai;
                a += Ai*faceNormals[facei];
                xc += Ai*(b[0]*pts[f[0]] + b[1]*pts[f[1]] + b[2]*pts[f[2]]);
            }

            if (A > VSMALL)
            {
                xc /= A;
                const vector& nf = faceNormals[faces_[k]];
                points_[k] = xc - ((xc - points_[k]) & nf)*nf;
            }

            const scalar magA = mag(a);
            areas_[k] = A;
            pressureScale_[k] = (A > VSMALL ? magA/A : 1);
            if (magA > VSMALL)
            {
                normals_[k] = a/magA;
            }
        }
    }
    else if (!body_.rigid())
    {
        // The area of each point of a deforming body changes with that of
        // its parent face
        forAll(faces_, i)
        {
            areas_[i] =
                areas0_[i]*mag(surf[faces_[i]].areaNormal(pts))
               /max(faceAreas0_[i], VSMALL);
        }
    }

    // The widths are found again only if the body moves
    if (body_.moving() || mesh_.moving() || h_.size() != points_.size())
    {
        h_ = cellWidths(points_);
    }
}


void Foam::immersedSurfaceTraction::wallValues
(
    const volVectorField& U,
    const volScalarField& p,
    const scalar nu,
    scalarField& pw,
    vectorField& tau
)
{
    updatePoints();

    // Surface and vertex velocities of a deforming body
    const triSurface& surf = body_.surface();
    const pointField& pts = surf.points();
    const vectorField& Ubp = body_.pointVelocities();
    const bool deforming = !body_.rigid() && Ubp.size() == pts.size();

    const label n = faces_.size();
    const vectorField& C = mesh_.C();

    // Solved directions
    labelList dirs;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solD_[d] == 1)
        {
            dirs.append(d);
        }
    }
    const label nD = dirs.size();

    // Quadratic basis relative to the wall point, without the constant for
    // the velocity (the velocity equals the body velocity at the wall point)
    const label nQ = nD*(nD + 1)/2;
    const label nU = nD + nQ;
    const label nP = 1 + nU;

    auto basis = [&](const vector& dx, const scalar h)
    {
        FixedList<scalar, 10> b(Zero);
        b[0] = 1;
        for (label k = 0; k < nD; ++k)
        {
            b[1 + k] = dx[dirs[k]]/h;
        }
        label q = 1 + nD;
        for (label k = 0; k < nD; ++k)
        {
            for (label l = k; l < nD; ++l)
            {
                b[q++] = b[1 + k]*b[1 + l];
            }
        }
        return b;
    };

    // Moments and right-hand sides of the pressure fit (full basis) and of
    // the velocity fit (without the constant) of each point, summed over the
    // processors
    const label nSum = nP*nP + nP + nU*nU + 3*nU;
    scalarField sums(n*nSum, Zero);

    const scalar radius = 3;

    // Local cells of the stencil of each point, found again only if the
    // body moves, and the rigid body velocity at their centres, from a single
    // evaluation
    if (body_.moving() || mesh_.moving() || stencils_.size() != n)
    {
        stencils_.setSize(n);
        forAll(faces_, i)
        {
            const scalar h = h_[i];

            if (h <= 0)
            {
                stencils_[i].clear();
                continue;
            }

            const point& xw = points_[i];
            const scalar r = radius*h;

            stencils_[i] =
                mesh_.cellTree().findBox
                (
                    treeBoundBox(xw - vector::one*r, xw + vector::one*r)
                );
        }
    }
    const List<labelList>& stencils = stencils_;

    labelList cellIndex(mesh_.nCells(), -1);
    DynamicList<label> stencilCells;
    forAll(stencils, i)
    {
        for (const label celli : stencils[i])
        {
            if (cellIndex[celli] < 0)
            {
                cellIndex[celli] = stencilCells.size();
                stencilCells.append(celli);
            }
        }
    }
    const vectorField stencilUb(body_.velocity(pointField(C, stencilCells)));

    forAll(faces_, i)
    {
        const scalar h = h_[i];

        if (h <= 0)
        {
            continue;
        }

        const point& xw = points_[i];
        const vector& nw = normals_[i];
        const scalar r = radius*h;

        const labelList& cells = stencils[i];

        scalar* Mp = &sums[i*nSum];
        scalar* Rp = Mp + nP*nP;
        scalar* Mu = Rp + nP;
        scalar* Ru = Mu + nU*nU;

        for (const label celli : cells)
        {
            const vector dx(C[celli] - xw);

            scalar rSqr = 0;
            for (const label d : dirs)
            {
                rSqr += sqr(dx[d]);
            }

            if (rSqr >= sqr(r))
            {
                continue;
            }

            // Fluid cells only: zero weight inside the body, rising to one
            // at h/2 from the tangent plane
            const scalar psi = dx & nw;
            const scalar s = min(max(psi/(0.5*h), scalar(0)), 1);
            const scalar w =
                sqr(s)*(3 - 2*s)*sqr(1 - rSqr/sqr(r))
               /(rSqr + sqr(0.5*h));

            if (w <= 0)
            {
                continue;
            }

            const FixedList<scalar, 10> b(basis(dx, h));

            // Velocity relative to the rigid body velocity at the cell
            // centre, so that a rigid rotation gives no traction
            const vector dU(U[celli] - stencilUb[cellIndex[celli]]);
            const scalar pc = p[celli];

            for (label j = 0; j < nP; ++j)
            {
                for (label k = 0; k < nP; ++k)
                {
                    Mp[j*nP + k] += w*b[j]*b[k];
                }
                Rp[j] += w*b[j]*pc;
            }

            for (label j = 0; j < nU; ++j)
            {
                for (label k = 0; k < nU; ++k)
                {
                    Mu[j*nU + k] += w*b[1 + j]*b[1 + k];
                }
                Ru[j] += w*b[1 + j]*dU.x();
                Ru[nU + j] += w*b[1 + j]*dU.y();
                Ru[2*nU + j] += w*b[1 + j]*dU.z();
            }
        }
    }

    // Native all-reduce of the contiguous sums
    reduce
    (
        sums.data(),
        sums.size(),
        sumOp<scalar>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    // Wall pressure (first) and viscous traction of each point, set on the
    // processor that solves its fits, and then summed over the processors
    scalarField values(4*n, Zero);

    // Solve a small weighted least squares system, whose matrix is
    // symmetric positive semi-definite, with a small regularisation of the
    // quadratic terms: by Cholesky decomposition, or with the pseudo-inverse
    // of the matrix if it is (nearly) singular, e.g. for collinear cell
    // centres
    auto solve = [](const scalar* M, const label m, const label nLin,
        List<List<scalar>>& rhs)
    {
        scalarRectangularMatrix A(m, m);
        scalar trace = 0;
        for (label j = 0; j < m; ++j)
        {
            for (label k = 0; k < m; ++k)
            {
                A(j, k) = M[j*m + k];
            }
            trace += M[j*m + j];
        }
        for (label j = nLin; j < m; ++j)
        {
            A(j, j) += 1e-6*trace/m;
        }

        // Cholesky decomposition A = L L^T, in the lower triangle of L
        scalarRectangularMatrix L(A);
        bool positive = true;
        for (label j = 0; j < m && positive; ++j)
        {
            scalar d = L(j, j);
            for (label k = 0; k < j; ++k)
            {
                d -= sqr(L(j, k));
            }
            if (d <= 1e-12*trace)
            {
                positive = false;
                break;
            }
            L(j, j) = Foam::sqrt(d);
            for (label i = j + 1; i < m; ++i)
            {
                scalar v = L(i, j);
                for (label k = 0; k < j; ++k)
                {
                    v -= L(i, k)*L(j, k);
                }
                L(i, j) = v/L(j, j);
            }
        }

        if (positive)
        {
            for (List<scalar>& b : rhs)
            {
                // L y = b, then L^T x = y
                for (label j = 0; j < m; ++j)
                {
                    for (label k = 0; k < j; ++k)
                    {
                        b[j] -= L(j, k)*b[k];
                    }
                    b[j] /= L(j, j);
                }
                for (label j = m - 1; j >= 0; --j)
                {
                    for (label k = j + 1; k < m; ++k)
                    {
                        b[j] -= L(k, j)*b[k];
                    }
                    b[j] /= L(j, j);
                }
            }
            return;
        }

        const scalarRectangularMatrix Ainv(SVD(A, 1e-10).VSinvUt());
        for (List<scalar>& b : rhs)
        {
            List<scalar> x(m, Zero);
            for (label j = 0; j < m; ++j)
            {
                for (label k = 0; k < m; ++k)
                {
                    x[j] += Ainv(j, k)*b[k];
                }
            }
            b = x;
        }
    };

    // The fits are solved for a share of the points on each processor, and
    // the tractions are then summed over the processors
    forAll(faces_, i)
    {
        const scalar h = h_[i];
        const scalar* Mp = &sums[i*nSum];
        const scalar* Rp = Mp + nP*nP;
        const scalar* Mu = Rp + nP;
        const scalar* Ru = Mu + nU*nU;

        // Too few fluid cells for a fit
        if
        (
            h <= 0
         || Mp[0] < SMALL
         || (i % Pstream::nProcs()) != Pstream::myProcNo()
        )
        {
            continue;
        }

        // Wall pressure
        List<List<scalar>> rp(1, List<scalar>(nP));
        for (label j = 0; j < nP; ++j)
        {
            rp[0][j] = Rp[j];
        }
        solve(Mp, nP, 1 + nD, rp);
        values[4*i] = rp[0][0];

        // Velocity gradient at the wall
        List<List<scalar>> ru(3, List<scalar>(nU));
        for (label c = 0; c < 3; ++c)
        {
            for (label j = 0; j < nU; ++j)
            {
                ru[c][j] = Ru[c*nU + j];
            }
        }
        solve(Mu, nU, nD, ru);

        tensor gradU(Zero);
        for (label c = 0; c < 3; ++c)
        {
            for (label k = 0; k < nD; ++k)
            {
                // d(U_c)/d(x_dirs[k])
                gradU(c, dirs[k]) = ru[c][k]/h;
            }
        }

        const vector& nw = normals_[i];

        if (deforming)
        {
            // On the surface of a deforming body, U equals the body
            // velocity Ub, whose surface gradient Gs is that of the linear
            // interpolation of the vertex velocities in the parent face, and
            // the fit gives the normal derivative a of U - Ub (Ub extended
            // along the normals). Then grad(U) = Gs + a n, continuity gives
            // a.n = -tr(Gs), and the viscous traction
            // nu*(grad(U) + grad(U)^T) & n is
            // nu*((I - nn) & a - 2 tr(Gs) n + Gs^T & n)
            const vector a(gradU & nw);
            vector tau((I - sqr(nw)) & a);

            const labelledTri& f = surf[faces_[i]];
            const vector e1(pts[f[1]] - pts[f[0]]);
            const vector e2(pts[f[2]] - pts[f[0]]);
            const vector A(e1 ^ e2);
            const scalar magSqrA = magSqr(A);

            if (magSqrA > VSMALL)
            {
                // Gradients of the barycentric coordinates in the face
                const vector g1((e2 ^ A)/magSqrA);
                const vector g2((A ^ e1)/magSqrA);
                const vector g0(-g1 - g2);

                // d(Ub_i)/dx_j, with the convention of gradU
                const tensor Gs
                (
                    (Ubp[f[0]]*g0) + (Ubp[f[1]]*g1) + (Ubp[f[2]]*g2)
                );

                tau += -2*tr(Gs)*nw + (Gs.T() & nw);
            }

            const vector taui(nu*tau);
            for (direction c = 0; c < vector::nComponents; ++c)
            {
                values[4*i + 1 + c] = taui[c];
            }
        }
        else
        {
            // Tangential viscous traction nu*(I - nn).(grad(U - Ub) & n): on
            // a no-slip wall the tangential derivatives of U - Ub vanish, and
            // the rigid body velocity Ub has no strain
            const vector taui(nu*((I - sqr(nw)) & (gradU & nw)));
            for (direction c = 0; c < vector::nComponents; ++c)
            {
                values[4*i + 1 + c] = taui[c];
            }
        }
    }

    reduce
    (
        values.data(),
        values.size(),
        sumOp<scalar>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    pw.setSize(n);
    tau.setSize(n);
    forAll(pw, i)
    {
        pw[i] = values[4*i];
        tau[i] = vector(values[4*i + 1], values[4*i + 2], values[4*i + 3]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedSurfaceTraction::immersedSurfaceTraction
(
    const immersedBody& body,
    const fvMesh& mesh,
    const quadratureType quadrature,
    const scalar spacingFactor
)
:
    mesh_(mesh),
    body_(body),
    solD_(mesh.solutionD()),
    quadrature_(quadrature),
    spacingFactor_(spacingFactor),
    faces_(),
    barycentric_(),
    areas_(),
    areas0_(),
    faceAreas0_(),
    h_(),
    points_(),
    normals_(),
    stencils_(),
    candidateFaces_(),
    candidateBarycentric_(),
    candidateAreas0_(),
    candidateAreas_(),
    candidateFaceAreas0_(),
    members_(),
    pressureScale_()
{
    if (spacingFactor_ <= 0)
    {
        FatalErrorInFunction
            << "The spacing factor of the clusters must be positive, not "
            << spacingFactor_ << exit(FatalError);
    }

    createPoints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::immersedSurfaceTraction::traction
(
    const volVectorField& U,
    const volScalarField& p,
    const scalar nu
)
{
    scalarField pw;
    vectorField tau;
    wallValues(U, p, nu, pw, tau);

    // The physical traction, with the unit normal of each point
    return -pw*normals_ + tau;
}


void Foam::immersedSurfaceTraction::faceValues
(
    const vectorField& pointValues,
    vectorField& faceValues,
    scalarField& faceAreas
) const
{
    if (pointValues.size() != faces_.size())
    {
        FatalErrorInFunction
            << "The field has " << pointValues.size() << " values, but "
            << "there are " << faces_.size() << " quadrature points"
            << abort(FatalError);
    }

    const label nFaces = body_.surface().size();

    faceValues.setSize(nFaces);
    faceValues = Zero;
    faceAreas.setSize(nFaces);
    faceAreas = 0;

    // The points outside the mesh have no value but keep their area, so
    // that the mean over a face partly outside the mesh gives the force of
    // its part in the mesh when it is applied to the whole face. The value
    // of a cluster applies to each of its candidates, with its area, on its
    // face
    if (members_.size())
    {
        forAll(members_, k)
        {
            for (const label i : members_[k])
            {
                const label facei = candidateFaces_[i];
                if (h_[k] > 0)
                {
                    faceValues[facei] += candidateAreas_[i]*pointValues[k];
                }
                faceAreas[facei] += candidateAreas_[i];
            }
        }
    }
    else
    {
        forAll(faces_, i)
        {
            if (h_[i] > 0)
            {
                faceValues[faces_[i]] += areas_[i]*pointValues[i];
            }
            faceAreas[faces_[i]] += areas_[i];
        }
    }

    forAll(faceValues, facei)
    {
        if (faceAreas[facei] > VSMALL)
        {
            faceValues[facei] /= faceAreas[facei];
        }
        else
        {
            faceValues[facei] = Zero;
            faceAreas[facei] = 0;
        }
    }
}


void Foam::immersedSurfaceTraction::force
(
    const volVectorField& U,
    const volScalarField& p,
    const scalar nu,
    const scalar rho,
    const point& CofR,
    vector& F,
    vector& T
)
{
    scalarField pw;
    vectorField tau;
    wallValues(U, p, nu, pw, tau);

    // The values are the same on every processor. The pressure is
    // integrated with the vector area of each point, |a| n = (|a|/A) A n,
    // and the viscous traction with its area A
    F = Zero;
    T = Zero;
    forAll(pw, i)
    {
        if (h_[i] > 0)
        {
            const vector t(-pw[i]*pressureScale_[i]*normals_[i] + tau[i]);
            const vector dF(rho*t*areas_[i]);
            F += dF;
            T += (points_[i] - CofR) ^ dF;
        }
    }
}


// ************************************************************************* //
