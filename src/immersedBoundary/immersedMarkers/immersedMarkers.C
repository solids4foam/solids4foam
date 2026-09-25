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

#include "immersedMarkers.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::immersedMarkers::phi(const scalar r)
{
    const scalar a = mag(r);

    if (a <= 0.5)
    {
        return (1 + Foam::sqrt(1 - 3*sqr(a)))/3;
    }
    else if (a <= 1.5)
    {
        return (5 - 3*a - Foam::sqrt(max(1 - 3*sqr(1 - a), scalar(0))))/6;
    }

    return 0;
}


Foam::scalar Foam::immersedMarkers::cellWidth(const label celli) const
{
    const vector span
    (
        boundBox(mesh_.points(), mesh_.cellPoints()[celli], false).span()
    );

    scalar w = GREAT;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solD_[d] == 1)
        {
            w = min(w, span[d]);
        }
    }

    return w;
}


Foam::scalarField Foam::immersedMarkers::cellWidths
(
    const pointField& pts
) const
{
    // polyMesh::findCell constructs the tet base points collectively
    (void)mesh_.tetBasePtIs();

    scalarField widths(pts.size(), Zero);

    forAll(pts, i)
    {
        const label celli = mesh_.findCell(pts[i]);

        if (celli >= 0)
        {
            widths[i] = cellWidth(celli);
        }
    }

    Pstream::listCombineReduce(widths, maxEqOp<scalar>());

    return widths;
}


void Foam::immersedMarkers::createMarkers()
{
    const triSurface& surf = body_.surface();
    const pointField& pts = surf.points();
    const boundBox& meshBb = mesh_.bounds();

    DynamicList<label> faces;
    DynamicList<FixedList<scalar, 3>> barycentric;
    DynamicList<scalar> areas;

    if (nD_ == 2)
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

        // Segment of each face, as the barycentric coordinates of its ends
        DynamicList<label> segFaces;
        DynamicList<FixedList<scalar, 3>> segA;
        DynamicList<FixedList<scalar, 3>> segB;
        DynamicList<point> segMid;
        DynamicList<scalar> segLength;

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
                    endPoints.append
                    (
                        (1 - t)*pts[f[k]] + t*pts[f[k1]]
                    );
                }
            }

            if (ends.size() == 2)
            {
                const scalar L = mag(endPoints[1] - endPoints[0]);

                if (L > VSMALL)
                {
                    segFaces.append(facei);
                    segA.append(ends[0]);
                    segB.append(ends[1]);
                    segMid.append(0.5*(endPoints[0] + endPoints[1]));
                    segLength.append(L);
                }
            }
        }

        const scalarField h(cellWidths(pointField(segMid)));

        forAll(segFaces, i)
        {
            if (h[i] <= 0)
            {
                continue;
            }

            const label n =
                max(label(1), label(std::ceil(segLength[i]/(spacing_*h[i]))));

            for (label k = 0; k < n; ++k)
            {
                const scalar t = (k + 0.5)/n;

                FixedList<scalar, 3> b;
                for (label j = 0; j < 3; ++j)
                {
                    b[j] = (1 - t)*segA[i][j] + t*segB[i][j];
                }

                faces.append(segFaces[i]);
                barycentric.append(b);
                areas.append(segLength[i]/n*thickness);
            }
        }
    }
    else
    {
        // Uniform subdivision of each face
        pointField centres(surf.size());
        forAll(surf, facei)
        {
            centres[facei] = surf[facei].centre(pts);
        }

        const scalarField h(cellWidths(centres));

        // Faces outside the mesh use the mean width of those inside
        scalar sumH = 0;
        label nH = 0;
        forAll(h, facei)
        {
            if (h[facei] > 0)
            {
                sumH += h[facei];
                ++nH;
            }
        }
        const scalar meanH = (nH > 0 ? sumH/nH : GREAT);

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

            const scalar hf = (h[facei] > 0 ? h[facei] : meanH);
            const label n =
                max(label(1), label(std::ceil(maxEdge/(spacing_*hf))));

            // Centroids of the n^2 sub-triangles
            for (label i = 0; i < n; ++i)
            {
                for (label j = 0; i + j < n; ++j)
                {
                    FixedList<scalar, 3> b;
                    b[0] = (i + 1.0/3.0)/n;
                    b[1] = (j + 1.0/3.0)/n;
                    b[2] = 1 - b[0] - b[1];
                    faces.append(facei);
                    barycentric.append(b);
                    areas.append(A/scalar(n*n));

                    if (i + j < n - 1)
                    {
                        b[0] = (i + 2.0/3.0)/n;
                        b[1] = (j + 2.0/3.0)/n;
                        b[2] = 1 - b[0] - b[1];
                        faces.append(facei);
                        barycentric.append(b);
                        areas.append(A/scalar(n*n));
                    }
                }
            }
        }
    }

    faces_.transfer(faces);
    barycentric_.transfer(barycentric);
    areas_.transfer(areas);

    Info<< "    Immersed body " << body_.name() << ": " << faces_.size()
        << " markers, total area " << sum(areas_) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedMarkers::immersedMarkers
(
    const immersedBody& body,
    const fvMesh& mesh,
    const scalar spacing,
    const scalar retraction
)
:
    mesh_(mesh),
    body_(body),
    spacing_(spacing),
    retraction_(retraction),
    solD_(mesh.solutionD()),
    nD_(mesh.nSolutionD()),
    faces_(),
    barycentric_(),
    areas_(),
    points_(),
    surfacePoints_(),
    normals_(),
    velocities_(),
    h_(),
    volumes_(),
    offsets_(),
    cells_(),
    weights_()
{
    createMarkers();
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedMarkers::update()
{
    const triSurface& surf = body_.surface();
    const pointField& pts = surf.points();
    const vectorField& faceNormals = surf.faceNormals();
    const label nMarkers = faces_.size();

    surfacePoints_.setSize(nMarkers);
    normals_.setSize(nMarkers);

    forAll(faces_, m)
    {
        const labelledTri& f = surf[faces_[m]];
        const FixedList<scalar, 3>& b = barycentric_[m];

        surfacePoints_[m] = b[0]*pts[f[0]] + b[1]*pts[f[1]] + b[2]*pts[f[2]];
        normals_[m] = faceNormals[faces_[m]];
    }

    h_ = cellWidths(surfacePoints_);
    points_ = surfacePoints_ - retraction_*h_*normals_;
    velocities_ = body_.velocity(surfacePoints_);
    volumes_ = areas_*h_;

    // Kernel stencils on this processor
    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();
    const boundBox& meshBb = mesh_.bounds();
    const label nP = 1 + nD_;

    DynamicList<label> offsets(nMarkers + 1);
    DynamicList<label> cells(27*nMarkers);
    DynamicList<scalar> weights(27*nMarkers);

    // Moment matrix of each marker, summed over the processors
    scalarField moments(nMarkers*nP*nP, Zero);

    // Active components of the position relative to a marker
    auto relative = [&](const vector& dx, const scalar h)
    {
        FixedList<scalar, 4> p(Zero);
        p[0] = 1;
        label k = 1;
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            if (solD_[d] == 1)
            {
                p[k++] = dx[d]/h;
            }
        }
        return p;
    };

    offsets.append(0);
    forAll(points_, m)
    {
        const scalar h = h_[m];

        if (h > 0)
        {
            point bbMin(points_[m] - 1.5*h*vector::one);
            point bbMax(points_[m] + 1.5*h*vector::one);
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                if (solD_[d] != 1)
                {
                    bbMin[d] = meshBb.min()[d];
                    bbMax[d] = meshBb.max()[d];
                }
            }

            const labelList candidates
            (
                mesh_.cellTree().findBox(treeBoundBox(bbMin, bbMax))
            );

            for (const label celli : candidates)
            {
                const vector dx(C[celli] - points_[m]);

                scalar w = V[celli];
                for (direction d = 0; d < vector::nComponents; ++d)
                {
                    if (solD_[d] == 1)
                    {
                        w *= phi(dx[d]/h);
                    }
                }

                if (w > 0)
                {
                    cells.append(celli);
                    weights.append(w);

                    const FixedList<scalar, 4> p(relative(dx, h));
                    for (label i = 0; i < nP; ++i)
                    {
                        for (label j = 0; j < nP; ++j)
                        {
                            moments[(m*nP + i)*nP + j] += w*p[i]*p[j];
                        }
                    }
                }
            }
        }

        offsets.append(cells.size());
    }

    Pstream::listCombineReduce(moments, plusEqOp<scalar>());

    // Reproducing kernel correction: w = w~ p.(M^-1 e0)
    forAll(points_, m)
    {
        if (h_[m] <= 0 || moments[m*nP*nP] <= VSMALL)
        {
            continue;
        }

        scalarSquareMatrix M(nP);
        for (label i = 0; i < nP; ++i)
        {
            for (label j = 0; j < nP; ++j)
            {
                M(i, j) = moments[(m*nP + i)*nP + j];
            }
        }

        const scalar M00 = moments[m*nP*nP];

        scalarField c(nP, Zero);
        c[0] = 1;
        LUsolve(M, c);

        // Constant correction if the moments are (nearly) singular, e.g. for
        // a stencil truncated by a boundary
        bool valid = true;
        for (const scalar ci : c)
        {
            if (!std::isfinite(ci))
            {
                valid = false;
            }
        }
        if (!valid || mag(c[0]*M00) > 10 || c[0]*M00 < 0.1)
        {
            c = Zero;
            c[0] = 1/M00;
        }

        for (label k = offsets[m]; k < offsets[m + 1]; ++k)
        {
            const FixedList<scalar, 4> p
            (
                relative(C[cells[k]] - points_[m], h_[m])
            );

            scalar cp = 0;
            for (label i = 0; i < nP; ++i)
            {
                cp += c[i]*p[i];
            }

            weights[k] *= cp;
        }
    }

    offsets_.transfer(offsets);
    cells_.transfer(cells);
    weights_.transfer(weights);
}


Foam::tmp<Foam::vectorField> Foam::immersedMarkers::interpolate
(
    const UList<vector>& vf
) const
{
    auto tvalues = tmp<vectorField>::New(faces_.size(), Zero);
    vectorField& values = tvalues.ref();

    forAll(values, m)
    {
        for (label k = offsets_[m]; k < offsets_[m + 1]; ++k)
        {
            values[m] += weights_[k]*vf[cells_[k]];
        }
    }

    Pstream::listCombineReduce(values, plusEqOp<vector>());

    return tvalues;
}


Foam::tmp<Foam::scalarField> Foam::immersedMarkers::interpolate
(
    const UList<scalar>& sf
) const
{
    auto tvalues = tmp<scalarField>::New(faces_.size(), Zero);
    scalarField& values = tvalues.ref();

    forAll(values, m)
    {
        for (label k = offsets_[m]; k < offsets_[m + 1]; ++k)
        {
            values[m] += weights_[k]*sf[cells_[k]];
        }
    }

    Pstream::listCombineReduce(values, plusEqOp<scalar>());

    return tvalues;
}


void Foam::immersedMarkers::spread
(
    const UList<vector>& F,
    UList<vector>& f
) const
{
    const scalarField& V = mesh_.V();

    forAll(F, m)
    {
        for (label k = offsets_[m]; k < offsets_[m + 1]; ++k)
        {
            const label celli = cells_[k];
            f[celli] += weights_[k]*F[m]*volumes_[m]/V[celli];
        }
    }
}


void Foam::immersedMarkers::spread
(
    const UList<scalar>& F,
    UList<scalar>& f
) const
{
    const scalarField& V = mesh_.V();

    forAll(F, m)
    {
        for (label k = offsets_[m]; k < offsets_[m + 1]; ++k)
        {
            const label celli = cells_[k];
            f[celli] += weights_[k]*F[m]*volumes_[m]/V[celli];
        }
    }
}


// ************************************************************************* //
