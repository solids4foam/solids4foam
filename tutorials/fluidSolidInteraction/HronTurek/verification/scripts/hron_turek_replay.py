#!/usr/bin/env python3
"""Replay a recorded coupled FSI3 flag trajectory on the CFD3 fluid meshes.

1. `hron_turek_replay_source.py` records the plate-patch point positions at
   every step of a coupled FSI3 run in its periodic state.
2. This script fits the recorded motion over an integer number of periods with
   a mean and harmonics 1..H in x and y, per point of the (open) polyline of
   the flag's reference surface, interpolates the coefficients along that
   polyline to the plate points of the target fluid mesh (any level whose
   points lie on the flag surface), writes them to constant/plateMotion.tab and
   builds a fluid-only case from the developed rigid CFD3 state of that level
   in which the plate point velocity is the finite difference of the replayed
   displacement (a smooth 1 s ramp from 25 s). Same mesh motion and fluid
   settings as the FSI3 and ALE tests, no solid, no coupling.

    python3 scripts/hron_turek_replay.py --level 4 --delta-t 0.00025 --cores 24
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402
import hron_turek_ale as ale  # noqa: E402

SOURCE = "replaysrc_iqnils_mesh_2x"
SOURCE_RUN = "iqnils_mesh_2x"
PERIODS = 5
HARMONICS = 4
START = 25.0
RAMP = 1.0
KEY = 1e5  # coefficient rows are keyed by the reference xy, rounded to 10 um


# ---------------------------------------------------------------------------
# OpenFOAM ASCII mesh readers
# ---------------------------------------------------------------------------

def _body(path: Path) -> str:
    text = path.read_text()
    start = text.index("\n(") + 2
    return text[start:text.rindex("\n)")]


def read_points(path: Path) -> np.ndarray:
    return np.array([[float(v) for v in m.groups()] for m in
                     re.finditer(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", _body(path))])


def read_faces(path: Path) -> list[list[int]]:
    return [[int(v) for v in m.group(2).split()] for m in
            re.finditer(r"(\d+)\(([\d ]+)\)", _body(path))]


def patch_range(boundary: Path, name: str) -> tuple[int, int]:
    text = boundary.read_text()
    m = re.search(rf"\b{name}\s*\{{[^}}]*?nFaces\s+(\d+);[^}}]*?startFace\s+(\d+);", text, re.S)
    if not m:
        driver.fail(f"Patch {name} not found in {boundary}")
    return int(m.group(2)), int(m.group(1))


def plate_polyline(mesh: Path) -> tuple[np.ndarray, dict]:
    """Ordered xy polyline of the flag surface, and xy-key -> vertex index."""
    pts = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    start, n = patch_range(mesh / "boundary", "plate")
    edges = set()
    for face in faces[start:start + n]:
        xy = {(round(pts[i][0] * KEY), round(pts[i][1] * KEY)) for i in face}
        if len(xy) == 2:
            edges.add(tuple(sorted(xy)))
    adjacency: dict = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    ends = [k for k, v in adjacency.items() if len(v) == 1]
    if len(ends) != 2:
        driver.fail(f"Plate polyline is not an open chain ({len(ends)} ends)")
    chain = [ends[0]]
    previous = None
    while True:
        nxt = [k for k in adjacency[chain[-1]] if k != previous]
        if not nxt:
            break
        previous = chain[-1]
        chain.append(nxt[0])
    return np.array(chain, float) / KEY, {k: i for i, k in enumerate(chain)}


# ---------------------------------------------------------------------------
# Source trajectory and fit
# ---------------------------------------------------------------------------

def read_source(case: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Times, reference polyline and displacement[t, vertex, 2] of the source."""
    mesh = case / "constant/fluid/polyMesh"
    ref = read_points(mesh / "points")
    path, index = plate_polyline(mesh)
    series: dict[int, np.ndarray] = {}
    times = None
    for processor in sorted(case.glob("processor*")):
        lines = (processor / "postProcessing/plateMotion.dat").read_text().splitlines()
        local = [int(v) for v in lines[0].split()[1:]]
        addressing = [int(v) for v in _body(processor / "constant/fluid/polyMesh/pointProcAddressing").split()]
        data = np.array([[float(v) for v in line.split()] for line in lines[1:]])
        if times is None:
            times = data[:, 0]
        n = min(len(times), len(data))
        times = times[:n]
        for k, lid in enumerate(local):
            gid = addressing[lid]
            series[gid] = data[:n, 1 + 2 * k:3 + 2 * k]
    disp = np.zeros((len(times), len(path), 2))
    count = np.zeros(len(path))
    for gid, xy in series.items():
        key = (round(ref[gid][0] * KEY), round(ref[gid][1] * KEY))
        if key not in index:
            continue
        v = index[key]
        disp[:, v, :] += xy[:len(times)] - ref[gid][:2]
        count[v] += 1
    if (count == 0).any():
        driver.fail("Some polyline vertices have no recorded points")
    disp /= count[None, :, None]
    return times, path, disp


def fit_trajectory(times: np.ndarray, path: np.ndarray, disp: np.ndarray) -> dict:
    # Tip: the polyline vertex nearest (0.6, 0.19), the lower tip corner
    tip = int(np.argmin(np.hypot(path[:, 0] - 0.6, path[:, 1] - 0.19)))
    signal = disp[:, tip, 1]
    keep = times > times[0] + 0.2  # skip the restart transient
    t, y = times[keep], signal[keep]
    stats = driver.periodic_statistics(list(t), list(y), t[-1] - t[0] - 1e-9)
    crossings = stats["crossings"]
    if len(crossings) < PERIODS + 1:
        driver.fail(f"Only {len(crossings) - 1} periods recorded")
    ta, tb = crossings[-PERIODS - 1], crossings[-1]
    omega = 2 * math.pi * PERIODS / (tb - ta)
    sel = (times >= ta) & (times <= tb)
    tau = times[sel] - ta
    basis = [np.ones_like(tau)]
    for h in range(1, HARMONICS + 1):
        basis += [np.cos(h * omega * tau), np.sin(h * omega * tau)]
    basis = np.array(basis).T
    coef = np.zeros((len(path), 2, 1 + 2 * HARMONICS))
    worst = 0.0
    worst_abs = 0.0
    for v in range(len(path)):
        for c in range(2):
            solution, *_ = np.linalg.lstsq(basis, disp[sel, v, c], rcond=None)
            coef[v, c] = solution
            residual = disp[sel, v, c] - basis @ solution
            worst = max(worst, residual.std() / max(disp[sel, v, c].std(), 1e-9))
            worst_abs = max(worst_abs, float(np.abs(residual).max()))
    amplitude = 0.5 * (y.max() - y.min())
    return {"coef": coef, "omega": omega, "frequency": omega / (2 * math.pi),
            "window": (ta, tb), "worst_relative_residual": worst, "max_abs_residual_m": worst_abs, "tip_vertex": tip,
            "tip_y_amplitude": 0.5 * (signal[sel].max() - signal[sel].min()),
            "tip_x_mean": float(disp[sel, tip, 0].mean())}


def interpolate(path: np.ndarray, coef: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Project target xy points on the polyline and interpolate coefficients."""
    out = np.zeros((len(target),) + coef.shape[1:])
    segment = path[1:] - path[:-1]
    length2 = (segment ** 2).sum(1)
    worst = 0.0
    for i, p in enumerate(target):
        w = np.clip(((p - path[:-1]) * segment).sum(1) / length2, 0, 1)
        foot = path[:-1] + w[:, None] * segment
        d = np.hypot(*(foot - p).T)
        k = int(np.argmin(d))
        worst = max(worst, d[k])
        out[i] = (1 - w[k]) * coef[k] + w[k] * coef[k + 1]
    if worst > 1e-5:
        driver.fail(f"Target plate point {worst:.2e} m off the source polyline")
    return out


def write_table(path: Path, reference: np.ndarray, coef: np.ndarray, omega: float) -> None:
    with path.open("w") as handle:
        handle.write(f"{omega:.12g} {len(reference)} {coef.shape[2]}\n")
        for xy, c in zip(reference, coef):
            handle.write(f"{round(xy[0] * KEY)} {round(xy[1] * KEY)} "
                         + " ".join(f"{v:.9g}" for v in c.reshape(-1)) + "\n")


# ---------------------------------------------------------------------------
# Case construction
# ---------------------------------------------------------------------------

BC_CODE = r"""
        static bool loaded = false;
        static std::map<std::pair<long, long>, std::vector<double> > table;
        static std::vector<point> ref;
        static scalar omega = 0;
        static label ncoef = 0;
        const scalar t = this->db().time().value();
        const scalar dt = this->db().time().deltaTValue();
        const pointField& pts = this->patch().localPoints();
        if (!loaded)
        {
            std::ifstream in
            (
                (this->db().time().globalPath()/"constant/plateMotion.tab").c_str()
            );
            long n = 0;
            in >> omega >> n >> ncoef;
            for (long i = 0; i < n; i++)
            {
                long kx, ky;
                in >> kx >> ky;
                std::vector<double> c(2*ncoef);
                for (label j = 0; j < 2*ncoef; j++) { in >> c[j]; }
                table[std::make_pair(kx, ky)] = c;
            }
            ref.resize(pts.size());
            forAll(pts, i) { ref[i] = pts[i]; }
            loaded = true;
        }
        const scalar t0 = %(t0)g, ramp = %(ramp)g;
        auto s = [&](scalar tt) -> scalar
        {
            scalar r = min(max((tt - t0)/ramp, 0.0), 1.0);
            return r*r*(3.0 - 2.0*r);
        };
        vectorField v(pts.size(), vector::zero);
        forAll(pts, i)
        {
            auto it = table.find
            (
                std::make_pair
                (
                    long(std::llround(ref[i].x()*%(key)g)),
                    long(std::llround(ref[i].y()*%(key)g))
                )
            );
            if (it == table.end())
            {
                FatalErrorInFunction << "No replay row for plate point "
                    << ref[i] << exit(FatalError);
            }
            const std::vector<double>& c = it->second;
            vector d1 = vector::zero, d0 = vector::zero;
            for (int comp = 0; comp < 2; comp++)
            {
                scalar a = c[comp*ncoef], b = c[comp*ncoef];
                const scalar tn = max(t - t0, 0.0), to = max(t - dt - t0, 0.0);
                for (label h = 1; h < ncoef/2 + 1; h++)
                {
                    a += c[comp*ncoef + 2*h - 1]*cos(h*omega*tn)
                       + c[comp*ncoef + 2*h]*sin(h*omega*tn);
                    b += c[comp*ncoef + 2*h - 1]*cos(h*omega*to)
                       + c[comp*ncoef + 2*h]*sin(h*omega*to);
                }
                d1[comp] = s(t)*a;
                d0[comp] = s(t - dt)*b;
            }
            v[i] = (d1 - d0)/dt;
        }
        operator==(v);
"""

POWER_FO = """
    platePower
    {
        type            coded;
        libs            ( utilityFunctionObjects );
        name            platePower;
        codeExecute
        #{
            const fvMesh& m = mesh();
            const label pid = m.boundaryMesh().findPatchID("plate");
            const volScalarField& p = m.lookupObject<volScalarField>("p");
            const volVectorField& U = m.lookupObject<volVectorField>("U");
            const scalarField& pp = p.boundaryField()[pid];
            const vectorField& Up = U.boundaryField()[pid];
            const vectorField& Sf = m.Sf().boundaryField()[pid];
            scalar w = 0;
            forAll(pp, i) { w += 1000.0*pp[i]*(Sf[i] & Up[i])/0.015; }
            reduce(w, sumOp<scalar>());
            Info<< "POWR " << m.time().value() << " " << w << endl;
        #};
    }
"""


def case_name(level: int, delta_t: float) -> str:
    return f"replay_{level}x_dt{delta_t:g}"


def build_case(level: int, delta_t: float, duration: float, cores: int,
               plate_pressure: str | None = None, tag: str = "_replay") -> Path:
    source = driver.WORK_ROOT / SOURCE
    times, path, disp = read_source(source)
    fit = fit_trajectory(times, path, disp)
    cfd3 = driver.WORK_ROOT / f"cfd3_{level}x_dt{delta_t:g}"
    target_path, _ = plate_polyline(cfd3 / "constant/polyMesh")
    target_coef = interpolate(path, fit["coef"], target_path)
    bc = ("plate\n    {\n        type            codedFixedValue;\n"
          "        value           uniform (0 0 0);\n        name            plateReplay;\n"
          "        codeInclude\n        #{\n            #include <fstream>\n"
          "            #include <map>\n            #include <vector>\n            #include <cmath>\n"
          "        #};\n"
          "        code\n        #{" + BC_CODE % dict(t0=START, ramp=RAMP, key=KEY)
          + "        #};\n    }")
    case = ale.build_case(level, delta_t, duration, cores, tag=tag,
                          plate_bc=bc, extra_functions=POWER_FO,
                          plate_pressure=plate_pressure)
    write_table(case / "constant/plateMotion.tab", target_path, target_coef, fit["omega"])
    summary = {k: v for k, v in fit.items() if k != "coef"}
    (case / "replay_fit.json").write_text(__import__("json").dumps(summary, indent=2, default=str))
    print(f"{case.name}: fit window {fit['window']}, f = {fit['frequency']:.4f} Hz, "
          f"tip dy amplitude {fit['tip_y_amplitude'] * 1e3:.2f} mm, "
          f"worst relative fit residual {fit['worst_relative_residual']:.3e}, max abs {fit['max_abs_residual_m'] * 1e3:.3f} mm")
    return case


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--level", type=int, required=True, choices=(1, 2, 4))
    parser.add_argument("--delta-t", type=float, required=True)
    parser.add_argument("--duration", type=float, default=8.0)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--setup-only", action="store_true")
    parser.add_argument("--pressure", default=None,
                        help="plate pressure condition instead of zeroGradient, "
                             "e.g. movingWallPressure (case tag _replay_<name>)")
    args = parser.parse_args()
    tag = "_replay" + (f"_{args.pressure}" if args.pressure else "")
    case = build_case(args.level, args.delta_t, args.duration, args.cores,
                      args.pressure, tag)
    if not args.setup_only:
        ale.run(case, args.cores)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
