#!/usr/bin/env python3
"""Generator of reference/results/controls/st-audit.csv: audit of the surface-
traction (ST) rotor torque of the cutLink method.

The recipe of the solver (immersedSurfaceTraction: the segments of the rotor
surface cut by the mid-plane of the mesh, thinned to a spacing of 0.7 h, the
weighted quadratic fits of the pressure and of the velocity relative to the
rigid rotation within 3 h of each point, the viscous traction
nu (I - n n) grad(U - U_rigid) n and the pressure -p n) is re-implemented here
in plain Python and applied to

 - the exact fields (u_theta = A r + B / r, the matching pressure, gauge
   p(r = 0.075) = 0 like the pressure reference of the tutorial), once on the
   inscribed polygon of the rotor surface (the wall of the solver) and once on
   the circle (the wall points moved to r = R1, normals radial);
 - the numerical fields at t = 10 s of the cutLink surface-traction run of the
   level (the case cutLink-ST-L<k> of ./Allverify);
 - the numerical fields with the cells whose centre is less than one cell
   width from the wall replaced by the exact fields (the first fluid layer).

Columns of st-audit.csv (all as a difference from the exact torque in percent,
100 (T / T_exact - 1)): recipe_exact_field_circle_wall_pct,
recipe_exact_field_polygon_wall_pct, recipe_numerical_field_pct,
solver_surface_traction_pct (the mean of the written torque of the run over
8 <= t <= 10 s) and first_fluid_layer_replaced_by_exact_pct. thinned_points
is the number of points kept by the thinning.

    ./Allverify run --runs L1 L2 L3 L4 --methods cutLink    # ME, ST and FO
    python3 scripts/taylor_couette_st_audit.py --out st-audit.csv

Levels whose case is missing in work/cases are skipped (Level 4 runs on four
ranks, 1.3 rank-h). Only the standard library and the driver are used; the
numerics are those of the research script that wrote the tracked file, with
the stencil search restricted to the cells near the point instead of a scan
of all cells (the same cells are selected). The exact-field columns need no
run at all.
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import immersed_taylor_couette_verification as v  # noqa: E402

P = v.parameters()
R1, R2, OMEGA, NU, RHO, LZ = (P["R1"], P["R2"], P["omega"], P["nu"], P["rho"],
                              P["Lz"])
A_, B_ = v.exact_AB(P)
T_EX = v.exact_rotor_torque(P)
GAUGE_RADIUS = 0.075              # pressure reference of the tutorial
COLUMNS = ["level", "thinned_points", "recipe_exact_field_circle_wall_pct",
           "recipe_exact_field_polygon_wall_pct", "recipe_numerical_field_pct",
           "solver_surface_traction_pct",
           "first_fluid_layer_replaced_by_exact_pct"]


def pc(torque: float) -> float:
    return 100.0 * (torque / T_EX - 1.0)


# ------------------------------------------------------------------ geometry
def read_triangles(path: Path) -> list:
    """[(v0, v1, v2)] of an ascii STL, each vertex (x, y, z)"""
    pts = [tuple(float(x) for x in m) for m in re.findall(
        r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", path.read_text())]
    return [tuple(pts[i:i + 3]) for i in range(0, len(pts), 3)]


def triangle_normal(tri: tuple) -> tuple:
    a, b, c = tri
    u = [b[k] - a[k] for k in range(3)]
    w = [c[k] - a[k] for k in range(3)]
    n = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
         u[0] * w[1] - u[1] * w[0])
    norm = math.sqrt(sum(x * x for x in n))
    return tuple(x / norm for x in n)


def midplane_points(tris: list, h: float, zmid: float = 0.05,
                    thick: float = 0.1) -> tuple:
    """The segments of the surface cut by z = zmid, split to a length of at
    most h: positions, facet indices and areas (length x thickness)"""
    pos, fac, area = [], [], []
    for fi, tri in enumerate(tris):
        ends = []
        for k in range(3):
            k1 = (k + 1) % 3
            dk, dk1 = tri[k][2] - zmid, tri[k1][2] - zmid
            if (dk < 0) != (dk1 < 0):
                t = dk / (dk - dk1)
                ends.append(tuple((1 - t) * tri[k][i] + t * tri[k1][i]
                                  for i in range(3)))
        if len(ends) != 2:
            continue
        length = math.dist(ends[0], ends[1])
        if length < 1e-300:
            continue
        n = max(1, math.ceil(length / h))
        for k in range(n):
            t = (k + 0.5) / n
            pos.append(tuple((1 - t) * ends[0][i] + t * ends[1][i]
                             for i in range(3)))
            fac.append(fi)
            area.append(length / n * thick)
    return pos, fac, area


def thin_points(pos: list, area: list, spacing: float) -> tuple:
    """Points scanned in order and kept if no kept point is closer than the
    spacing; the area of every candidate goes to the nearest kept point"""
    def d2(a, b):
        return sum((a[i] - b[i]) ** 2 for i in range(3))
    kept = []
    for i, x in enumerate(pos):
        if not kept or min(d2(pos[k], x) for k in kept) >= spacing ** 2:
            kept.append(i)
    kareas = [0.0] * len(kept)
    for i, x in enumerate(pos):
        j = min(range(len(kept)), key=lambda q: d2(pos[kept[q]], x))
        kareas[j] += area[i]
    return kept, kareas


def solver_points(tris: list, level: int) -> tuple:
    """Points, facet indices and areas of the solver (thinning 0.7 h)"""
    n, h, _ = v.cell_centres(level)
    pos, fac, area = midplane_points(tris, h)
    kept, kareas = thin_points(pos, area, 0.7 * h)
    return [pos[i] for i in kept], [fac[i] for i in kept], kareas, h


# -------------------------------------------------------------- exact fields
def u_exact(x: float, y: float) -> tuple:
    f = A_ + B_ / (x * x + y * y)             # u_theta / r
    return -f * y, f * x


def p_exact(r: float) -> float:
    return A_ ** 2 * r ** 2 / 2 + 2 * A_ * B_ * math.log(r) - B_ ** 2 / (2 * r ** 2)


def p_gauge(r: float) -> float:
    return p_exact(r) - p_exact(GAUGE_RADIUS)


def exact_fields(level: int) -> tuple:
    """Velocity (ux, uy) and the pressure with the gauge of the tutorial, at
    the cell centres"""
    _, _, centres = v.cell_centres(level)
    return ([u_exact(x, y) for x, y in centres],
            [p_gauge(math.hypot(x, y)) for x, y in centres])


# ------------------------------------------------------------ linear algebra
def solve(a: list, b: list) -> list:
    """Gaussian elimination with partial pivoting; b is a list of rows (one
    value per right-hand side)"""
    n = len(a)
    m = [list(a[i]) + list(b[i]) for i in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[piv] = m[piv], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for k in range(c, len(m[r])):
                m[r][k] -= f * m[c][k]
    x = [None] * n
    for r in range(n - 1, -1, -1):
        x[r] = [(m[r][n + q] - sum(m[r][k] * x[k][q] for k in range(r + 1, n)))
                / m[r][r] for q in range(len(b[0]))]
    return x


# -------------------------------------------------------------- traction fit
def stencil(level: int, point: tuple, radius_cells: float = 3.0) -> list:
    """Indices of the cells with a centre closer than radius_cells * h to the
    point (in the plane), in the order of the mesh (x runs fastest)"""
    n, h, _ = v.cell_centres(level)
    side = v.SIDE
    reach = radius_cells * h
    i0 = max(0, int(math.floor((point[0] - reach + side / 2) / h - 0.5)) - 1)
    i1 = min(n - 1, int(math.floor((point[0] + reach + side / 2) / h - 0.5)) + 2)
    j0 = max(0, int(math.floor((point[1] - reach + side / 2) / h - 0.5)) - 1)
    j1 = min(n - 1, int(math.floor((point[1] + reach + side / 2) / h - 0.5)) + 2)
    out = []
    for j in range(j0, j1 + 1):
        for i in range(i0, i1 + 1):
            x = -side / 2 + (i + 0.5) * h
            y = -side / 2 + (j + 0.5) * h
            if (x - point[0]) ** 2 + (y - point[1]) ** 2 < reach ** 2:
                out.append(j * n + i)
    return out


def traction_fit(level: int, vel: list, pres: list, pts: list, normals: list,
                 radius: float = 3.0) -> list:
    """The traction of the solver at every point (non-deforming branch):
    (-p_wall n, tau, p_wall, cells in the stencil)"""
    n, h, centres = v.cell_centres(level)
    out = []
    for xw, nw in zip(pts, normals):
        cells = stencil(level, xw, radius)
        rows = []
        for ci in cells:
            cx, cy = centres[ci]
            dx, dy = cx - xw[0], cy - xw[1]
            r2 = dx * dx + dy * dy
            psi = dx * nw[0] + dy * nw[1]
            s = min(max(psi / (0.5 * h), 0.0), 1.0)
            w = s * s * (3 - 2 * s) * (1 - r2 / (radius * h) ** 2) ** 2 \
                / (r2 + (0.5 * h) ** 2)
            if w > 0:
                rows.append((w, dx / h, dy / h, ci, cx, cy))
        basis = [(1.0, d[1], d[2], d[1] ** 2, d[1] * d[2], d[2] ** 2)
                 for d in rows]
        # pressure: full quadratic basis, quadratic terms regularised
        m = [[sum(d[0] * b[i] * b[j] for d, b in zip(rows, basis))
              for j in range(6)] for i in range(6)]
        tr = sum(m[i][i] for i in range(6))
        for j in range(3, 6):
            m[j][j] += 1e-6 * tr / 6
        rhs = [[sum(d[0] * b[i] * pres[d[3]] for d, b in zip(rows, basis))]
               for i in range(6)]
        pw = solve(m, rhs)[0][0]
        # velocity relative to the rigid rotation, quadratic basis without
        # the constant term (zero at the wall point)
        bu = [b[1:] for b in basis]
        mu = [[sum(d[0] * b[i] * b[j] for d, b in zip(rows, bu))
               for j in range(5)] for i in range(5)]
        tru = sum(mu[i][i] for i in range(5))
        for j in range(2, 5):
            mu[j][j] += 1e-6 * tru / 5
        du = [(vel[d[3]][0] - OMEGA * (-d[5]), vel[d[3]][1] - OMEGA * d[4])
              for d in rows]
        rhs = [[sum(d[0] * b[i] * q[c] for d, b, q in zip(rows, bu, du))
                for c in range(2)] for i in range(5)]
        ru = solve(mu, rhs)                       # 5 x 2
        # grad[c][d] = dU_c / dx_d = coefficient of x_d in component c over h
        grad = [[ru[d][c] / h for d in range(2)] for c in range(2)]
        a = [grad[c][0] * nw[0] + grad[c][1] * nw[1] for c in range(2)]
        an = a[0] * nw[0] + a[1] * nw[1]
        tau = (NU * (a[0] - an * nw[0]), NU * (a[1] - an * nw[1]))
        out.append(((-pw * nw[0], -pw * nw[1]), tau, pw, len(rows)))
    return out


def torque_from(res: list, pts: list, areas: list) -> float:
    """Total torque about the axis, rho sum (x t_y - y t_x) A"""
    total = 0.0
    for (tp, tv, _, _), x, area in zip(res, pts, areas):
        tx, ty = tp[0] + tv[0], tp[1] + tv[1]
        total += RHO * (x[0] * ty - x[1] * tx) * area
    return total


def wall_torque(level: int, tris: list, vel: list, pres: list,
                wall: str = "polygon", points: tuple = None) -> float:
    """The torque of the recipe on the polygon wall of the solver or on the
    circle (points moved to R1, radial normals)"""
    pts, fac, areas, _ = points or solver_points(tris, level)
    normals = [triangle_normal(tris[f]) for f in fac]
    if wall == "circle":
        radii = [math.hypot(q[0], q[1]) for q in pts]
        pts = [(q[0] * R1 / r, q[1] * R1 / r, q[2]) for q, r in zip(pts, radii)]
        normals = [(q[0] / R1, q[1] / R1, 0.0) for q in pts]
    return torque_from(traction_fit(level, vel, pres, pts, normals), pts, areas)


# ---------------------------------------------------------------- run fields
def read_scalar_field(path: Path) -> list:
    text = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)"
                  r"\s*\n?\(", text)
    if not m:
        v.fail(f"{path}: not a nonuniform scalar field")
    body = text[m.end():]
    return [float(x) for x in body[:body.index("\n)\n")].split()]


def solver_surface_traction(case: Path) -> float:
    """Mean of the written rotor torque over 8 <= t <= 10 s (plain mean of
    the samples)"""
    rotor = v.read_force_file(case / "postProcessing/immersedBoundary/0/rotor.dat")
    sel = [q[2] for t, q in zip(rotor["t"], rotor["torque"]) if 8 <= t <= 10]
    return sum(sel) / len(sel)


def audit_level(level: int, tris: list, case: Path = None) -> dict:
    """The row of a level; the columns that need the run are empty when the
    case is None"""
    points = solver_points(tris, level)
    ue, pe = exact_fields(level)
    row = {"level": level, "thinned_points": len(points[0]),
           "recipe_exact_field_circle_wall_pct": pc(wall_torque(
               level, tris, ue, pe, "circle", points)),
           "recipe_exact_field_polygon_wall_pct": pc(wall_torque(
               level, tris, ue, pe, "polygon", points))}
    if case is not None:
        un = [q[:2] for q in v.read_vector_field(case / "10/U")]
        pn = read_scalar_field(case / "10/p")
        row["recipe_numerical_field_pct"] = pc(wall_torque(
            level, tris, un, pn, "polygon", points))
        row["solver_surface_traction_pct"] = pc(solver_surface_traction(case))
        n, h, centres = v.cell_centres(level)
        near = [(math.hypot(x, y) - R1) / h < 1.0 for x, y in centres]
        mixed_u = [ue[i] if near[i] else un[i] for i in range(len(un))]
        mixed_p = [pe[i] if near[i] else pn[i] for i in range(len(pn))]
        row["first_fluid_layer_replaced_by_exact_pct"] = pc(wall_torque(
            level, tris, mixed_u, mixed_p, "polygon", points))
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="st-audit.csv (default: stdout)")
    ap.add_argument("--cases-dir", type=Path, default=v.CASES_DIR)
    ap.add_argument("--levels", nargs="+", type=int, default=[1, 2, 3, 4])
    args = ap.parse_args(argv)
    tris = read_triangles(v.TUTORIAL_DIR / "constant/triSurface/rotor.stl")
    rows = []
    for level in args.levels:
        case = args.cases_dir / f"cutLink-ST-L{level}"
        have = (case / "10/U").exists() and (case / "10/p").exists()
        if not have:
            print(f"Level {level}: no completed {case.name}, the columns that "
                  "need the run are empty", file=sys.stderr)
        rows.append(audit_level(level, tris, case if have else None))
    out = sys.stdout if args.out is None else open(args.out, "w", newline="")
    try:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            # three decimals, two for the last column (as the tracked file)
            w.writerow([(f"{r[c]:.2f}" if c == COLUMNS[-1] else f"{r[c]:.3f}")
                        if isinstance(r.get(c), float) else r.get(c, "")
                        for c in COLUMNS])
    finally:
        if args.out is not None:
            out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
