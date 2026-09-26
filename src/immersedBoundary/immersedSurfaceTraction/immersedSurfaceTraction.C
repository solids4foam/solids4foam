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
#include <unordered_map>
#include <vector>
#include <cstdint>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::immersedSurfaceTraction::cellWidth(const label celli) const
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
            widths[i] = cellWidth(celli);
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
    const pointField& pts = surf.points();
    const boundBox& meshBb = mesh_.bounds();

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
            const std::int64_t i = std::int64_t(std::floor(x.x()/spacing)) + di;
            const std::int64_t j = std::int64_t(std::floor(x.y()/spacing)) + dj;
            const std::int64_t k = std::int64_t(std::floor(x.z()/spacing)) + dk;
            return (i*73856093) ^ (j*19349663) ^ (k*83492791);
        };

        std::unordered_map<std::int64_t, std::vector<label>> grid;
        DynamicList<label> kept(positions.size());

        // Nearest kept point within the spacing, -1 if none
        auto nearestKept = [&](const point& x, const scalar maxDist)
        {
            label best = -1;
            scalar bestDistSqr = sqr(maxDist);
            for (label di = -1; di <= 1; ++di)
            {
                for (label dj = -1; dj <= 1; ++dj)
                {
                    for (label dk = -1; dk <= 1; ++dk)
                    {
                        const auto iter = grid.find(key(x, di, dj, dk));
                        if (iter == grid.end())
                        {
                            continue;
                        }
                        for (const label k : iter->second)
                        {
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
            if (nearestKept(positions[i], spacing) < 0)
            {
                grid[key(positions[i], 0, 0, 0)].push_back(kept.size());
                kept.append(i);
            }
        }

        scalarField keptAreas(kept.size(), Zero);
        forAll(positions, i)
        {
            const label k = nearestKept(positions[i], 2*spacing);
            keptAreas[k >= 0 ? k : 0] += areas[i];
        }

        faces_.setSize(kept.size());
        barycentric_.setSize(kept.size());
        forAll(kept, k)
        {
            faces_[k] = faces[kept[k]];
            barycentric_[k] = barycentric[kept[k]];
        }
        areas_.transfer(keptAreas);
    }

    Info<< "    Immersed body " << body_.name() << ": " << faces_.size()
        << " traction points, area " << sum(areas_) << endl;
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

    // The widths are found again only if the body moves
    if (body_.moving() || mesh_.moving() || h_.size() != points_.size())
    {
        h_ = cellWidths(points_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedSurfaceTraction::immersedSurfaceTraction
(
    const immersedBody& body,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    body_(body),
    solD_(mesh.solutionD()),
    faces_(),
    barycentric_(),
    areas_(),
    h_(),
    points_(),
    normals_(),
    stencils_()
{
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
    updatePoints();

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

    tmp<vectorField> ttraction(new vectorField(n, Zero));
    vectorField& t = ttraction.ref();

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
        const scalar pw = rp[0][0];

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

        // Tangential viscous traction nu*(I - nn).(grad(U - Ub) & n): on a
        // no-slip wall the tangential derivatives of U - Ub vanish, and the
        // rigid body velocity Ub has no strain
        const vector& nw = normals_[i];
        const vector tau(nu*((I - sqr(nw)) & (gradU & nw)));

        t[i] = -pw*nw + tau;
    }

    reduce
    (
        reinterpret_cast<scalar*>(t.data()),
        vector::nComponents*t.size(),
        sumOp<scalar>(),
        UPstream::msgType(),
        UPstream::worldComm
    );

    return ttraction;
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
    const vectorField t(traction(U, p, nu));

    // The traction is the same on every processor
    F = Zero;
    T = Zero;
    forAll(t, i)
    {
        if (h_[i] > 0)
        {
            const vector dF(rho*t[i]*areas_[i]);
            F += dF;
            T += (points_[i] - CofR) ^ dF;
        }
    }
}


// ************************************************************************* //
