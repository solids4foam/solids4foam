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

Application
    Test-immersedSurfaceTraction

Description
    Manufactured checks of the quadratures of immersedSurfaceTraction
    (tractionQuadrature none, greedy and clustered) on closed surfaces
    inside the mesh, without a flow solution: the pressure is set to
    constant and linear fields, whose least-squares wall fits are exact, and
    the velocity to that of the body, so that the force is the quadrature
    of -p n alone:

        F = -sum_k w_k p(x_k) n_k

    For a closed surface inside the mesh (in two dimensions, the contour of
    the mesh mid-plane times the thickness of the mesh), with V the volume
    it encloses, the exact integrals are

        integral(n dA) = 0,   integral(x_j n_i dA) = V delta_ij,

    so the force of p = 1 is zero, that of p = x_j - c_j is -V e_j, and a
    constant added to the pressure does not change the force. The quadrature
    none (every candidate point) is exact for them; the force of clustered
    is compared with it and with these invariants. The surfaces (a polygonal
    cylinder and an ellipse extruded across a two-dimensional mesh, an
    ellipsoid in a three-dimensional mesh) are written by the test in three
    facet orders (as generated, reversed, and permuted with cyclically
    rotated vertices), static, rotated, translated and bent (quadraticBend),
    together with two closely spaced surfaces (a folded thin strip and two
    stacked thin plates, 1.5 mm apart, against clusters of 5 mm),
    and the clustered quadrature is also made after the motion (as on a
    restart) and compared with the one made in the reference configuration.

    Checks (PASS/FAIL lines; tolerances relative to the scales given):
      - none and clustered: |F(p = 1)| <= 1e-12 A (vector area);
      - clustered: |F(p = x + 1e3) - F(p = x)| <= 1e-12 1e3 A (constant
        pressure invariance);
      - none: |M - V' I| <= 1e-10 V', M_ij = -F_i(p = x_j - c_j),
        V' = tr(M)/nD (isotropy), and V' = V to 1e-8 for the rigid motions
        (of the single-precision points of the STL file);
      - clustered: |M - M_none| <= 1e-2 V and |T - T_none| <= 1e-2 V L for
        p = x - c (second-order cluster error);
      - clustered: total and face areas of faceValues equal to those of
        none to 1e-12 A;
      - clustered, two dimensions: at most 0.6 times the points of none
        (the candidates are much closer than the clusters; the mesh
        mid-plane, z = 0, lies on a plane of the cluster grid);
      - none and clustered: forces, torques and point counts independent
        of the facet order to 1e-12 (relative to A, A l and A l^2, with l the
        span of the surface);
      - clustered made after the motion equal to clustered made before.
    The greedy quadrature is reported (INFO) only: its vector area and
    moments depend on the facet order.

Author
    Sairam Pamulaparthi Venkata, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "immersedBody.H"
#include "immersedSurfaceTraction.H"
#include "OFstream.H"
#include <map>
#include <numeric>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

label nPassed = 0;
label nFailed = 0;

void check(const string& name, const scalar error, const scalar tolerance)
{
    const bool ok = (error <= tolerance);
    Info<< (ok ? "PASS: " : "FAIL: ") << name.c_str() << ": " << error
        << " (tolerance " << tolerance << ")" << nl;
    if (ok)
    {
        ++nPassed;
    }
    else
    {
        ++nFailed;
    }
}

typedef FixedList<point, 3> tri3;

// A point rounded to single precision, as the STL reader keeps it, so that
// the volumes computed here are those of the surface read back
point stlPoint(const point& p)
{
    return point(float(p.x()), float(p.y()), float(p.z()));
}

// A polygon with its edges subdivided into pieces of at most ds
List<point> subdivided(const List<point>& corners, const scalar ds)
{
    DynamicList<point> poly;
    forAll(corners, i)
    {
        const point& a = corners[i];
        const point& b = corners[(i + 1) % corners.size()];
        const label n = max(label(1), label(std::ceil(mag(b - a)/ds)));
        for (label k = 0; k < n; ++k)
        {
            poly.append(a + (b - a)*scalar(k)/n);
        }
    }
    return List<point>(std::move(poly));
}

// Prism over a counter-clockwise polygon, from z0 to z1: two side triangles
// per edge, then the bottom and top cap triangles
List<tri3> prism(const List<point>& poly, const scalar z0, const scalar z1)
{
    point c(Zero);
    for (const point& p : poly)
    {
        c += p;
    }
    c /= poly.size();
    const point C0(stlPoint(point(c.x(), c.y(), z0)));
    const point C1(stlPoint(point(c.x(), c.y(), z1)));

    DynamicList<tri3> tris;
    forAll(poly, i)
    {
        const point a(stlPoint(poly[i]));
        const point b(stlPoint(poly[(i + 1) % poly.size()]));
        const point A0(a.x(), a.y(), z0), B0(b.x(), b.y(), z0);
        const point A1(a.x(), a.y(), z1), B1(b.x(), b.y(), z1);
        tris.append(tri3({A0, B0, B1}));
        tris.append(tri3({A0, B1, A1}));
        tris.append(tri3({C0, B0, A0}));
        tris.append(tri3({C1, A1, B1}));
    }
    return List<tri3>(std::move(tris));
}

// Area of a polygon (counter-clockwise > 0), of its single-precision
// vertices
scalar polygonArea(const List<point>& poly)
{
    scalar a = 0;
    forAll(poly, i)
    {
        const point p(stlPoint(poly[i]));
        const point q(stlPoint(poly[(i + 1) % poly.size()]));
        a += p.x()*q.y() - q.x()*p.y();
    }
    return 0.5*a;
}

// Ellipsoid from an icosahedron subdivided 'level' times, with the axes
// 'axes', rotated by R about 'centre'
List<tri3> ellipsoid
(
    const label level,
    const vector& axes,
    const tensor& R,
    const point& centre
)
{
    const scalar t = (1 + Foam::sqrt(5.0))/2;
    DynamicList<point> v
    ({
        point(-1, t, 0), point(1, t, 0), point(-1, -t, 0), point(1, -t, 0),
        point(0, -1, t), point(0, 1, t), point(0, -1, -t), point(0, 1, -t),
        point(t, 0, -1), point(t, 0, 1), point(-t, 0, -1), point(-t, 0, 1)
    });
    for (point& p : v)
    {
        p /= mag(p);
    }
    DynamicList<FixedList<label, 3>> f
    ({
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    });
    for (label l = 0; l < level; ++l)
    {
        std::map<std::pair<label, label>, label> mid;
        auto midpoint = [&](const label i, const label j)
        {
            const std::pair<label, label> k(min(i, j), max(i, j));
            auto iter = mid.find(k);
            if (iter != mid.end())
            {
                return iter->second;
            }
            point m(v[i] + v[j]);
            v.append(m/mag(m));
            mid[k] = v.size() - 1;
            return v.size() - 1;
        };
        DynamicList<FixedList<label, 3>> nf;
        for (const FixedList<label, 3>& g : f)
        {
            const label a = midpoint(g[0], g[1]);
            const label b = midpoint(g[1], g[2]);
            const label c = midpoint(g[2], g[0]);
            nf.append(FixedList<label, 3>({g[0], a, c}));
            nf.append(FixedList<label, 3>({g[1], b, a}));
            nf.append(FixedList<label, 3>({g[2], c, b}));
            nf.append(FixedList<label, 3>({a, b, c}));
        }
        f.transfer(nf);
    }
    List<tri3> tris(f.size());
    forAll(f, i)
    {
        for (label k = 0; k < 3; ++k)
        {
            const point& p = v[f[i][k]];
            tris[i][k] = stlPoint
            (
                centre
              + (R & vector(axes.x()*p.x(), axes.y()*p.y(), axes.z()*p.z()))
            );
        }
    }
    return tris;
}

// Volume enclosed by a closed triangulated surface
scalar enclosedVolume(const List<tri3>& tris)
{
    scalar V = 0;
    for (const tri3& t : tris)
    {
        V += (t[0] & (t[1] ^ t[2]))/6;
    }
    return V;
}

// Write the triangles in the order 'order', each with its vertices rotated
// cyclically by 'rotate[i]'
void writeStl
(
    const fileName& file,
    const List<tri3>& tris,
    const labelList& order,
    const labelList& rotate
)
{
    OFstream os(file);
    os.precision(17);
    os << "solid body" << nl;
    forAll(order, j)
    {
        const tri3& t = tris[order[j]];
        const label r = rotate[j];
        const point& a = t[r % 3];
        const point& b = t[(r + 1) % 3];
        const point& c = t[(r + 2) % 3];
        vector n((b - a) ^ (c - a));
        n /= max(mag(n), VSMALL);
        os  << "  facet normal " << n.x() << ' ' << n.y() << ' ' << n.z() << nl
            << "    outer loop" << nl;
        for (const point* p : {&a, &b, &c})
        {
            os  << "      vertex " << p->x() << ' ' << p->y() << ' '
                << p->z() << nl;
        }
        os  << "    endloop" << nl << "  endfacet" << nl;
    }
    os << "endsolid body" << nl;
}

// Forces of the pressure fields and the face areas of one quadrature
struct result
{
    label nPoints = 0;
    scalar area = 0;
    scalar size = 0;                // span of the quadrature points
    scalarField faceAreas;
    vector F1 = Zero;               // p = 1
    List<vector> Fx;                // p = x_j - c_j, j = 0, nD - 1
    vector FxOff = Zero;            // p = x_0 - c_0 + 1e3
    vector Tx = Zero;               // torque about c of p = x_0 - c_0
    tensor M = Zero;                // M_ij = -F_i(p = x_j - c_j)
};

result evaluate
(
    immersedSurfaceTraction& tr,
    const immersedBody& body,
    const fvMesh& mesh,
    const labelList& dirs,
    const point& c
)
{
    result r;
    r.nPoints = tr.size();
    r.size = mag(boundBox(tr.points(), true).span());

    volVectorField U
    (
        IOobject("U", mesh.time().timeName(), mesh),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    );
    U.primitiveFieldRef() = body.velocity(mesh.C());
    volScalarField p
    (
        IOobject("p", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(sqr(dimVelocity), Zero)
    );
    const vectorField& C = mesh.C();

    auto force = [&](const scalarField& pc, vector& F, vector& T)
    {
        p.primitiveFieldRef() = pc;
        tr.force(U, p, 0, 1, c, F, T);
    };

    vector T;
    force(scalarField(mesh.nCells(), 1), r.F1, T);

    r.Fx.setSize(dirs.size());
    forAll(dirs, j)
    {
        const scalarField pc(C.component(dirs[j]) - c[dirs[j]]);
        vector Tj;
        force(pc, r.Fx[j], Tj);
        if (j == 0)
        {
            r.Tx = Tj;
        }
        for (direction i = 0; i < vector::nComponents; ++i)
        {
            r.M(i, dirs[j]) = -r.Fx[j][i];
        }
    }
    force(C.component(dirs[0]) - c[dirs[0]] + 1e3, r.FxOff, T);

    // Face areas (and values) of a unit point value
    vectorField faceValues;
    tr.faceValues(vectorField(tr.size(), vector::one), faceValues, r.faceAreas);
    r.area = sum(r.faceAreas);

    return r;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label nD = mesh.nGeometricD();
    labelList dirs;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (mesh.solutionD()[d] == 1)
        {
            dirs.append(d);
        }
    }
    if (nD == 2 && mesh.solutionD()[vector::Z] != -1)
    {
        FatalErrorInFunction
            << "For a two-dimensional mesh the empty direction must be z"
            << abort(FatalError);
    }

    const boundBox& bb = mesh.bounds();
    const point c0(bb.centre());
    const scalar pi = constant::mathematical::pi;

    // Shapes: name, triangles, enclosed volume in the mesh, size
    DynamicList<word> names;
    DynamicList<List<tri3>> shapes;
    DynamicList<scalar> volumes;
    if (nD == 2)
    {
        const scalar z0 = bb.min().z() - 0.5*bb.span().z();
        const scalar z1 = bb.max().z() + 0.5*bb.span().z();
        const scalar thickness = bb.span().z();

        List<point> circle(128);
        forAll(circle, i)
        {
            const scalar a = 2*pi*i/circle.size();
            circle[i] = c0 + 0.05*vector(Foam::cos(a), Foam::sin(a), 0);
        }
        names.append("cylinder");
        shapes.append(prism(circle, z0, z1));
        volumes.append(polygonArea(circle)*thickness);

        List<point> ellipse(96);
        const scalar q = 20*pi/180;
        forAll(ellipse, i)
        {
            const scalar a = 2*pi*i/ellipse.size();
            const scalar x = 0.06*Foam::cos(a);
            const scalar y = 0.035*Foam::sin(a);
            ellipse[i] =
                c0
              + vector
                (
                    x*Foam::cos(q) - y*Foam::sin(q),
                    x*Foam::sin(q) + y*Foam::cos(q),
                    0
                );
        }
        names.append("ellipse");
        shapes.append(prism(ellipse, z0, z1));
        volumes.append(polygonArea(ellipse)*thickness);

        // Closely spaced surfaces, 1.5 mm apart (less than the cluster
        // size): a folded thin strip (hairpin), whose legs have faces of the
        // same orientation 3 mm apart, and two thin plates stacked in one
        // surface file
        const scalar t = 1.5e-3;
        const scalar L = 0.08;
        {
            const point o(c0 - vector(t + 0.5*t, 0.5*L, 0));
            const List<point> hairpin
            (
                subdivided
                (
                    {
                        o, o + vector(3*t, 0, 0), o + vector(3*t, L, 0),
                        o + vector(2*t, L, 0), o + vector(2*t, t, 0),
                        o + vector(t, t, 0), o + vector(t, L, 0),
                        o + vector(0, L, 0)
                    },
                    1e-3
                )
            );
            names.append("hairpin");
            shapes.append(prism(hairpin, z0, z1));
            volumes.append(polygonArea(hairpin)*thickness);
        }
        {
            List<tri3> plates;
            scalar V = 0;
            for (const scalar y : {-1.5*t, 0.5*t})
            {
                const point o(c0 + vector(-0.5*L, y, 0));
                const List<point> plate
                (
                    subdivided
                    (
                        {
                            o, o + vector(L, 0, 0), o + vector(L, t, 0),
                            o + vector(0, t, 0)
                        },
                        1e-3
                    )
                );
                plates.append(prism(plate, z0, z1));
                V += polygonArea(plate)*thickness;
            }
            names.append("stackedPlates");
            shapes.append(plates);
            volumes.append(V);
        }
    }
    else
    {
        const scalar a = 25*pi/180;
        const scalar b = 40*pi/180;
        const tensor Rz
        (
            Foam::cos(a), -Foam::sin(a), 0,
            Foam::sin(a), Foam::cos(a), 0,
            0, 0, 1
        );
        const tensor Rx
        (
            1, 0, 0,
            0, Foam::cos(b), -Foam::sin(b),
            0, Foam::sin(b), Foam::cos(b)
        );
        names.append("ellipsoid");
        shapes.append(ellipsoid(3, vector(0.06, 0.045, 0.035), Rz & Rx, c0));
        volumes.append(enclosedVolume(shapes.last()));
    }

    // Motions, evaluated at t = 1
    const scalar t = 1;
    DynamicList<word> motionNames;
    DynamicList<dictionary> motions;
    {
        dictionary d;
        d.add("type", word("static"));
        motionNames.append("static");
        motions.append(d);
    }
    {
        dictionary d;
        d.add("type", word("solidBodyRotation"));
        d.add("origin", c0 + vector(0.01, -0.02, 0));
        d.add("axis", nD == 2 ? vector(0, 0, 1) : vector(1, 1, 1));
        d.add("omega", 0.7);
        motionNames.append("rotation");
        motions.append(d);
    }
    if (nD == 2)
    {
        dictionary d;
        d.add("type", word("sinusoidalTranslation"));
        d.add("amplitude", 0.05);
        d.add("period", 4.0);
        d.add("direction", vector(1, 0.3, 0));
        motionNames.append("translation");
        motions.append(d);
    }
    {
        dictionary d;
        d.add("type", word("quadraticBend"));
        d.add("amplitude", 0.02);
        d.add("period", 4.0);
        d.add("yMin", c0.y() - 0.07);
        d.add("yMax", c0.y() + 0.07);
        d.add("direction", vector(1, 0, 0));
        motionNames.append("bend");
        motions.append(d);
    }

    const scalar exact = 1e-12;
    const scalar clusterTol = 1e-2;
    mkDir(runTime.constant()/"triSurface");

    forAll(shapes, si)
    {
        const List<tri3>& tris = shapes[si];
        const label nTris = tris.size();

        // Facet orders: as generated, reversed, and permuted (stride
        // coprime with the number of triangles) with rotated vertices
        const wordList orders({"original", "reversed", "permuted"});
        List<labelList> order(3, labelList(nTris));
        List<labelList> rotate(3, labelList(nTris, Zero));
        label stride = 37;
        while (std::gcd(stride, nTris) != 1)
        {
            ++stride;
        }
        for (label i = 0; i < nTris; ++i)
        {
            order[0][i] = i;
            order[1][i] = nTris - 1 - i;
            order[2][i] = (i*stride) % nTris;
            rotate[2][i] = i % 3;
        }
        forAll(orders, o)
        {
            writeStl
            (
                runTime.constant()/"triSurface"
               /(names[si] + "-" + orders[o] + ".stl"),
                tris,
                order[o],
                rotate[o]
            );
        }

        forAll(motions, mi)
        {
            const string tag(names[si] + ", " + motionNames[mi]);
            Info<< nl << "=== " << tag.c_str() << nl;

            const bool rigid = (motionNames[mi] != "bend");
            List<result> none(orders.size());
            List<result> clustered(orders.size());

            forAll(orders, o)
            {
                dictionary bodyDict;
                bodyDict.add
                (
                    "surface",
                    fileName(names[si] + "-" + orders[o] + ".stl")
                );
                bodyDict.add("motion", motions[mi]);
                immersedBody body
                (
                    word(names[si] + "_" + motionNames[mi] + "_" + orders[o]),
                    bodyDict,
                    mesh
                );

                // Clustered, made in the initial configuration and moved
                body.move(0);
                immersedSurfaceTraction clustered0
                (
                    body, mesh, word("clustered")
                );
                body.move(t);

                // The others made in the current configuration (clustered
                // too, as on a restart)
                immersedSurfaceTraction noneTr(body, mesh, word("none"));
                immersedSurfaceTraction greedyTr(body, mesh, word("greedy"));
                immersedSurfaceTraction clusteredT
                (
                    body, mesh, word("clustered")
                );

                const point c(gAverage(body.surface().points()));
                const result rn(evaluate(noneTr, body, mesh, dirs, c));
                const result rc(evaluate(clustered0, body, mesh, dirs, c));
                const result rcT(evaluate(clusteredT, body, mesh, dirs, c));
                const result rg(evaluate(greedyTr, body, mesh, dirs, c));
                none[o] = rn;
                clustered[o] = rc;

                const scalar A = rn.area;
                const scalar Vn = tr(rn.M)/nD;
                const scalar L =
                (
                    nD == 2
                  ? A/(2*pi*bb.span().z())
                  : Foam::sqrt(A/(4*pi))
                );
                const string ot(tag + ", " + orders[o]);

                if (o == 0)
                {
                    Info<< "points: none " << rn.nPoints << ", greedy "
                        << rg.nPoints << ", clustered " << rc.nPoints
                        << "; area " << A << ", V " << Vn << nl;
                }

                // none: exact invariants
                check(ot + ": none |F(p=1)|/A", mag(rn.F1)/A, exact);
                {
                    tensor dM(rn.M);
                    for (const label d : dirs)
                    {
                        dM(d, d) -= Vn;
                    }
                    check(ot + ": none |M - V I|/V", mag(dM)/Vn, 1e-10);
                }
                if (rigid)
                {
                    check
                    (
                        ot + ": none |V - V_exact|/V",
                        mag(Vn - volumes[si])/volumes[si],
                        1e-8
                    );
                }

                // clustered: exact invariants
                check(ot + ": clustered |F(p=1)|/A", mag(rc.F1)/A, exact);
                check
                (
                    ot + ": clustered |F(p=x+1e3) - F(p=x)|/(1e3 A)",
                    mag(rc.FxOff - rc.Fx[0])/(1e3*A),
                    exact
                );
                check
                (
                    ot + ": clustered |area - area_none|/A",
                    mag(rc.area - A)/A,
                    exact
                );
                check
                (
                    ot + ": clustered max face |area - area_none|/A",
                    max(mag(rc.faceAreas - rn.faceAreas))/A,
                    exact
                );

                // clustered: the candidate pieces of the two-dimensional
                // surfaces are much closer than the cluster size, so the
                // clusters must merge them
                if (nD == 2)
                {
                    check
                    (
                        ot + ": clustered points/none points",
                        scalar(rc.nPoints)/rn.nPoints,
                        0.6
                    );
                }

                // clustered: second-order quadrature error
                check
                (
                    ot + ": clustered |M - M_none|/V",
                    mag(rc.M - rn.M)/Vn,
                    clusterTol
                );
                check
                (
                    ot + ": clustered |T - T_none|/(V L), p = x",
                    mag(rc.Tx - rn.Tx)/(Vn*L),
                    clusterTol
                );

                // clustered made after the motion (restart)
                scalar dRestart = mag(rcT.F1 - rc.F1)/A;
                forAll(rc.Fx, j)
                {
                    dRestart = max(dRestart, mag(rcT.Fx[j] - rc.Fx[j])/Vn);
                }
                dRestart = max(dRestart, mag(rcT.Tx - rc.Tx)/(Vn*L));
                check
                (
                    ot + ": clustered made after the motion, max difference",
                    dRestart + mag(rcT.nPoints - rc.nPoints),
                    exact
                );

                Info<< "INFO: " << ot.c_str() << ": greedy |F(p=1)|/A "
                    << mag(rg.F1)/A << ", |M - M_none|/V "
                    << mag(rg.M - rn.M)/Vn << ", |T - T_none|/(V L) "
                    << mag(rg.Tx - rn.Tx)/(Vn*L) << "; clustered |M - "
                    << "M_none|/V " << mag(rc.M - rn.M)/Vn << nl;
            }

            // Independence of the facet order
            for (const bool isClustered : {false, true})
            {
                const List<result>& r = isClustered ? clustered : none;
                const word q(isClustered ? "clustered" : "none");
                // Relative to the force and torque scales of a unit pressure
                // and pressure gradient over the surface
                const scalar A = r[0].area;
                const scalar l = r[0].size;
                scalar d = 0;
                for (label o = 1; o < orders.size(); ++o)
                {
                    d = max(d, mag(r[o].F1 - r[0].F1)/A);
                    forAll(r[0].Fx, j)
                    {
                        d = max(d, mag(r[o].Fx[j] - r[0].Fx[j])/(A*l));
                    }
                    d = max(d, mag(r[o].Tx - r[0].Tx)/(A*sqr(l)));
                    d = max(d, scalar(mag(r[o].nPoints - r[0].nPoints)));
                }
                check(tag + ": " + q + " facet-order spread", d, exact);
            }
        }
    }

    Info<< nl << "Checks: " << nPassed << " passed, " << nFailed << " failed"
        << nl << endl;

    return (nFailed == 0 ? 0 : 1);
}


// ************************************************************************* //
