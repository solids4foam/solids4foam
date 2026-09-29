#!/usr/bin/env python3
"""Split the Hessenthaler et al. (2017) fluid-domain surface into patches.

The published fluid-domain surface (geometry/fluidDomain.vtk.gz, millimetres,
CC0) is one closed triangulation. This script writes a multi-solid ASCII STL
in metres with the patches

    upperInlet  inlet disk at z = -29.5 mm with y > 0
    lowerInlet  inlet disk at z = -29.5 mm with y < 0
    outlet      outlet disk, moved to the end of the outlet extension
    wall        the phantom wall and the outlet extension
    interface   the five faces of the flap cavity (11 x 2 x 65 mm)

and extends the outlet pipe (diameter 76.2 mm) by --extension metres, as in
Hessenthaler, Roehrle and Nordsletten (2017), who found that 50 mm changed
the flow and used 250 mm.

Only the Python 3 standard library is used.
"""

import argparse
import gzip
import math
from pathlib import Path


def read_vtk(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        tokens = handle.read().split()
    i = tokens.index("POINTS")
    n_points = int(tokens[i + 1])
    values = [float(v) for v in tokens[i + 3:i + 3 + 3 * n_points]]
    points = [tuple(values[3 * k:3 * k + 3]) for k in range(n_points)]
    i = tokens.index("POLYGONS")
    n_polys = int(tokens[i + 1])
    j = i + 3
    triangles = []
    for _ in range(n_polys):
        if tokens[j] != "3":
            raise SystemExit("Only triangles are supported")
        triangles.append(tuple(int(v) for v in tokens[j + 1:j + 4]))
        j += 4
    return points, triangles


def classify(points, triangles, tol=1.0e-3):
    """Return the patch name of every triangle."""
    z_in = min(p[2] for p in points)
    z_out = max(p[2] for p in points)
    names = []
    for tri in triangles:
        c = [sum(points[v][k] for v in tri) / 3.0 for k in range(3)]
        on_x = abs(abs(c[0]) - 5.5) < tol and abs(c[1]) <= 1.0 + tol
        on_y = abs(abs(c[1]) - 1.0) < tol and abs(c[0]) <= 5.5 + tol
        on_tip = abs(c[2] - 65.0) < tol and abs(c[0]) <= 5.5 + tol \
            and abs(c[1]) <= 1.0 + tol
        in_z = -tol <= c[2] <= 65.0 + tol
        if abs(c[2] - z_in) < tol:
            names.append("upperInlet" if c[1] > 0 else "lowerInlet")
        elif abs(c[2] - z_out) < tol:
            names.append("outlet")
        elif in_z and (on_x or on_y or on_tip):
            names.append("interface")
        else:
            names.append("wall")
    return names, z_out


def rim_loop(triangles, names, patch):
    """Ordered boundary loop of the vertices of one planar patch."""
    count = {}
    for tri, name in zip(triangles, names):
        if name != patch:
            continue
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (min(a, b), max(a, b))
            count[key] = count.get(key, 0) + 1
    edges = [e for e, n in count.items() if n == 1]
    neighbours = {}
    for a, b in edges:
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)
    if any(len(v) != 2 for v in neighbours.values()):
        raise SystemExit(f"The {patch} rim is not a single simple loop")
    start = edges[0][0]
    loop = [start]
    previous, current = None, start
    while True:
        nxt = [v for v in neighbours[current] if v != previous][0]
        if nxt == start:
            break
        loop.append(nxt)
        previous, current = current, nxt
    if len(loop) != len(neighbours):
        raise SystemExit(f"The {patch} rim has more than one loop")
    return loop


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", default=here / "geometry/fluidDomain.vtk.gz")
    parser.add_argument("--output",
                        default=here / "constant/triSurface/fluidDomain.stl")
    parser.add_argument("--extension", type=float, default=0.25,
                        help="outlet extension length in metres (default 0.25)")
    parser.add_argument("--segments", type=int, default=25,
                        help="axial triangle rows along the extension")
    args = parser.parse_args()

    points, triangles = read_vtk(args.input)
    names, z_out = classify(points, triangles)

    faces = {name: [] for name in
             ("upperInlet", "lowerInlet", "outlet", "wall", "interface")}
    for tri, name in zip(triangles, names):
        if name != "outlet":
            faces[name].append(tuple(points[v] for v in tri))
    if not faces["interface"] or len(faces["interface"]) < 10:
        raise SystemExit("Flap interface faces not found")

    # Outlet extension: a tube from the outlet rim, then a new outlet disk
    length = 1000.0 * args.extension
    loop = rim_loop(triangles, names, "outlet")
    rim = [points[v] for v in loop]
    # Keep the loop orientation consistent with the outward normal (+z)
    area_z = sum(a[0] * b[1] - b[0] * a[1]
                 for a, b in zip(rim, rim[1:] + rim[:1]))
    if area_z < 0.0:
        rim.reverse()
    rows = max(1, args.segments)
    rings = [[(p[0], p[1], z_out + length * k / rows) for p in rim]
             for k in range(rows + 1)]
    n = len(rim)
    if length > 0.0:
        for k in range(rows):
            lower, upper = rings[k], rings[k + 1]
            for i in range(n):
                j = (i + 1) % n
                # Outward normals point away from the pipe axis
                faces["wall"].append((lower[i], lower[j], upper[j]))
                faces["wall"].append((lower[i], upper[j], upper[i]))
    top = rings[-1]
    centre = (sum(p[0] for p in top) / n, sum(p[1] for p in top) / n,
              top[0][2])
    for i in range(n):
        faces["outlet"].append((centre, top[i], top[(i + 1) % n]))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        for name, tris in faces.items():
            handle.write(f"solid {name}\n")
            for tri in tris:
                a, b, c = ([1.0e-3 * x for x in p] for p in tri)
                u = [b[k] - a[k] for k in range(3)]
                v = [c[k] - a[k] for k in range(3)]
                nrm = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                       u[0] * v[1] - u[1] * v[0]]
                mag = math.sqrt(sum(x * x for x in nrm)) or 1.0
                handle.write("  facet normal {:.6e} {:.6e} {:.6e}\n"
                             .format(*(x / mag for x in nrm)))
                handle.write("    outer loop\n")
                for p in (a, b, c):
                    handle.write("      vertex {:.9e} {:.9e} {:.9e}\n"
                                 .format(*p))
                handle.write("    endloop\n  endfacet\n")
            handle.write(f"endsolid {name}\n")
    summary = ", ".join(f"{k} {len(v)}" for k, v in faces.items())
    print(f"Wrote {out} ({summary} triangles; outlet at "
          f"z = {1.0e-3 * top[0][2]:.4f} m)")


if __name__ == "__main__":
    main()
