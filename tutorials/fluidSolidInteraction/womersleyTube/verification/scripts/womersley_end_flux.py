#!/usr/bin/env python3
"""Mass-flux inconsistency on the tube ends of womersleyTube runs.

For each run under verification/work/temporal given on the command line,
reads the final fluid fields and reports, on the inlet and outlet patches,

    mean |phi_b - U_b . S_f|  /  max |U_b . S_f|

the difference between the face flux used by the continuity equation and
the flux of the boundary velocity used by the momentum equation. With the
tutorial's codedMixed end condition (valueFraction 0), OpenFOAM's
constrainHbyA sets HbyA_b = U_b (mixed patches are not assignable), so
phi_b = U_b . S_f - rAU_b snGrad(p) |S_f| and the difference is
rAU snGrad(p) |S_f|, with rAU ~ dt/1.5: an O(dt) inconsistency.

Writes verification/temporal/end_flux_inconsistency.csv.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

VERIFICATION = Path(__file__).resolve().parents[1]
RUN_ROOT = VERIFICATION / "work" / "temporal"


def boundary_values(path: Path, patch: str):
    text = path.read_text()
    text = text[text.index("boundaryField"):]
    start = text.index(f"\n    {patch}\n")
    block = text[start:text.index("\n    }\n", start)]
    match = re.search(r"value\s+nonuniform List<(\w+)>\s*\n?(\d+)\s*\n?\(",
                      block)
    if not match:
        raise ValueError(f"No nonuniform value for {patch} in {path}")
    n = int(match.group(2))
    body = block[match.end():]
    if match.group(1) == "scalar":
        return [float(x) for x in
                re.findall(r"[-+\d.eE]+", body.split(")")[0])[:n]]
    return [tuple(map(float, v.split()))
            for v in re.findall(r"\(([^)]*)\)", body)[:n]]


def face_areas(mesh: Path, final: Path, patch: str):
    points_file = final / "polyMesh" / "points"
    if not points_file.exists():
        points_file = mesh / "points"
    points = [tuple(map(float, v.split())) for v in re.findall(
        r"\(([-\d.eE+]+ [-\d.eE+]+ [-\d.eE+]+)\)", points_file.read_text())]
    faces = [list(map(int, f.split())) for f in
             re.findall(r"\d+\(([\d ]+)\)", (mesh / "faces").read_text())]
    boundary = (mesh / "boundary").read_text()
    block = boundary[boundary.index(f"\n    {patch}\n"):]
    first = int(re.search(r"startFace\s+(\d+)", block).group(1))
    count = int(re.search(r"nFaces\s+(\d+)", block).group(1))
    areas = []
    for face in faces[first:first + count]:
        p = [points[i] for i in face]
        s = [0.0, 0.0, 0.0]
        for i, a in enumerate(p):
            b = p[(i + 1) % len(p)]
            s[0] += a[1] * b[2] - a[2] * b[1]
            s[1] += a[2] * b[0] - a[0] * b[2]
            s[2] += a[0] * b[1] - a[1] * b[0]
        areas.append([x / 2 for x in s])
    return areas


def measure(case: Path) -> list[dict]:
    region = "fluid" if (case / "constant" / "fluid").is_dir() else ""
    times = sorted((float(d.name), d) for d in case.iterdir()
                   if re.fullmatch(r"[0-9.eE+-]+", d.name) and d.name != "0")
    final = times[-1][1] / region if region else times[-1][1]
    mesh = case / "constant" / region / "polyMesh" if region \
        else case / "constant" / "polyMesh"
    rows = []
    for patch in ("inlet", "outlet"):
        phi = boundary_values(final / "phi", patch)
        u = boundary_values(final / "U", patch)
        areas = face_areas(mesh, final, patch)
        flux_u = [sum(a * b for a, b in zip(uf, sf))
                  for uf, sf in zip(u, areas)]
        diff = [p - f for p, f in zip(phi, flux_u)]
        rows.append({"run": case.name, "patch": patch, "time": times[-1][0],
                     "mean_abs_mismatch_over_max_flux":
                     sum(abs(d) for d in diff) / len(diff)
                     / max(abs(f) for f in flux_u)})
    return rows


def main() -> int:
    rows = []
    for name in sys.argv[1:]:
        try:
            rows += measure(RUN_ROOT / name)
        except (ValueError, OSError) as error:
            print(f"{name}: {error}")
    out = VERIFICATION / "temporal" / "end_flux_inconsistency.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['run']:55s} {row['patch']:7s} "
              f"{row['mean_abs_mismatch_over_max_flux']:.3e}")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
