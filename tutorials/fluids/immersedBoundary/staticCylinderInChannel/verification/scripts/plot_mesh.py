#!/usr/bin/env python3
"""Close-up of the cylinder on the background mesh of a generated case.

Draws the actual cells of the case (the front faces of its polyMesh, the
same mesh for all the methods at a given level), the section of the
cylinder surface (constant/triSurface/cylinder.stl) and, with --fields, the
cells that the cutLink method penalises, from its recorded penalisation
rate immersedBoundary:kappa: the cells whose centre is inside the body
(kappa = kappa_s, the largest rate) and their fluid neighbours across the
cut links (0 < kappa < kappa_s). It reads existing files only and runs no
OpenFOAM application.

Usage (from verification/):
    python3 scripts/plot_mesh.py                   # work/final/cutLink-L2
    python3 scripts/plot_mesh.py --case work/final/cutLink-L3 --name mesh_L3

Writes postProcessing/figures/mesh/<name>.pdf and .png (300 dpi) and,
with --preview, a 150 dpi PNG (the tracked images/<name>.png).
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

VERIFICATION_DIR = Path(__file__).resolve().parent.parent

# Cylinder of the tutorial (Schafer and Turek 2D-1)
CENTRE = (0.2, 0.2)
DIAMETER = 0.1


def foam_list(path: Path, kind: str):
    """Entries of the (ASCII) list of an OpenFOAM file"""
    text = path.read_text()
    body = text[text.index("\n(", text.index("FoamFile")) + 2:]
    body = body[:body.rindex("\n)")]
    if kind == "vector":
        return [tuple(map(float, m))
                for m in re.findall(r"\(([^ ()]+) ([^ ()]+) ([^ ()]+)\)",
                                    body)]
    if kind == "face":
        return [tuple(map(int, m.split()))
                for m in re.findall(r"\d+\(([^()]*)\)", body)]
    return [int(v) for v in body.split()]


def patch_range(boundary: Path, patch: str):
    text = boundary.read_text()
    block = re.search(rf"\b{patch}\s*\{{([^}}]*)\}}", text)
    if not block:
        sys.exit(f"ERROR: no {patch} patch in {boundary}")
    n = int(re.search(r"nFaces\s+(\d+);", block.group(1)).group(1))
    start = int(re.search(r"startFace\s+(\d+);", block.group(1)).group(1))
    return start, n


def cell_field(path: Path):
    text = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(",
                  text)
    if not m:
        sys.exit(f"ERROR: {path} has no nonuniform internalField")
    n = int(m.group(1))
    values = text[m.end():].split(")", 1)[0].split()
    if len(values) != n:
        sys.exit(f"ERROR: {path}: expected {n} values")
    return [float(v) for v in values]


def stl_section(path: Path, z0: float):
    """Segments of the intersection of an ASCII STL with the plane z = z0"""
    vertices = [tuple(map(float, m)) for m in re.findall(
        r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", path.read_text())]
    segments = []
    for i in range(0, len(vertices), 3):
        tri = vertices[i:i + 3]
        cut = []
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            if (a[2] - z0)*(b[2] - z0) < 0:
                s = (z0 - a[2])/(b[2] - a[2])
                cut.append((a[0] + s*(b[0] - a[0]), a[1] + s*(b[1] - a[1])))
        if len(cut) == 2:
            segments.append(cut)
    return len(vertices)//3, segments


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--case", default="work/final/cutLink-L2",
                        help="case directory, relative to verification/")
    parser.add_argument("--time", default="",
                        help="time directory of the fields (default latest)")
    parser.add_argument("--name", default="mesh_L2_closeup")
    parser.add_argument("--half-width", type=float, default=0.7,
                        help="half-width of the view, in diameters")
    parser.add_argument("--no-fields", dest="fields", action="store_false",
                        help="draw the mesh and surface only")
    parser.add_argument("--preview", action="store_true",
                        help="also write images/<name>.png (150 dpi)")
    args = parser.parse_args()

    case = (VERIFICATION_DIR / args.case).resolve()
    mesh = case / "constant/polyMesh"
    if not (mesh / "points").exists():
        sys.exit(f"ERROR: no mesh in {case} (run blockMesh in the case)")

    points = foam_list(mesh / "points", "vector")
    faces = foam_list(mesh / "faces", "face")
    owner = foam_list(mesh / "owner", "label")
    start, n = patch_range(mesh / "boundary", "frontAndBack")
    z0 = min(p[2] for p in points)
    z1 = max(p[2] for p in points)

    # Cells as their faces on the back plane, in diameters about the centre
    hw = args.half_width
    polys, cells = [], []
    for facei in range(start, start + n):
        face = faces[facei]
        if any(abs(points[p][2] - z0) > 1e-9 for p in face):
            continue
        xy = [((points[p][0] - CENTRE[0])/DIAMETER,
               (points[p][1] - CENTRE[1])/DIAMETER) for p in face]
        if all(abs(x) > hw + 0.1 or abs(y) > hw + 0.1 for x, y in xy):
            continue
        polys.append(xy)
        cells.append(owner[facei])
    widths = sorted({round(max(x for x, _ in p) - min(x for x, _ in p), 9)
                     for p in polys})
    h = widths[0]*DIAMETER

    # Penalised cells from the recorded cutLink penalisation rate
    colours = ["white"]*len(polys)
    counts = None
    if args.fields:
        times = sorted((d for d in case.iterdir() if d.is_dir()
                        and re.fullmatch(r"[0-9.eE+-]+", d.name)
                        and (d / "immersedBoundary:kappa").exists()),
                       key=lambda d: float(d.name))
        if args.time:
            times = [case / args.time]
        if not times:
            sys.exit(f"ERROR: no immersedBoundary:kappa in {case}; run "
                     "with --no-fields or reconstruct the case")
        kappa = cell_field(times[-1] / "immersedBoundary:kappa")
        kmax = max(kappa)
        n_pen = sum(1 for k in kappa if k == kmax)
        n_link = sum(1 for k in kappa if 0 < k < kmax)
        counts = (n_pen, n_link, times[-1].name)
        for i, c in enumerate(cells):
            if kappa[c] == kmax:
                colours[i] = "#c9ccd1"
            elif kappa[c] > 0:
                colours[i] = "#9ecae1"

    n_tri, segments = stl_section(case / "constant/triSurface/cylinder.stl",
                                  0.5*(z0 + z1))
    segments = [[((x - CENTRE[0])/DIAMETER, (y - CENTRE[1])/DIAMETER)
                 for x, y in s] for s in segments]

    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 10, "axes.linewidth": 0.8,
    })
    fig, ax = plt.subplots(figsize=(4.6, 5.3))
    ax.add_collection(PolyCollection(polys, facecolors=colours,
                                     edgecolors="#4d4d4d", linewidths=0.35))
    ax.add_collection(LineCollection(segments, colors="#c0392b",
                                     linewidths=1.6, capstyle="round"))
    ax.set_xlim(-hw, hw)
    ax.set_ylim(-hw, hw)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$(x - x_c)/D$")
    ax.set_ylabel(r"$(y - y_c)/D$")
    ax.tick_params(direction="in", top=True, right=True, length=3)

    level = re.search(r"-L(\d)", case.name)
    ax.set_title(
        (f"L{level.group(1)}: " if level else "")
        + rf"$D/h = {DIAMETER/h:.0f}$, $h = {h:g}$ m", fontsize=10.5)

    handles = [Line2D([], [], color="#c0392b", lw=1.6,
                      label=f"cylinder surface (STL section, {len(segments)} "
                            "segments)")]
    if counts:
        handles += [
            Patch(facecolor="#c9ccd1", edgecolor="#4d4d4d", lw=0.35,
                  label=rf"centre inside: $\kappa = \kappa_s$ "
                        f"({counts[0]} cells)"),
            Patch(facecolor="#9ecae1", edgecolor="#4d4d4d", lw=0.35,
                  label=rf"cut-link fluid cells: $0 < \kappa < \kappa_s$ "
                        f"({counts[1]})"),
        ]
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.13), frameon=False, fontsize=9,
              handlelength=1.6)
    fig.tight_layout()

    out = VERIFICATION_DIR / "postProcessing/figures/mesh"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{args.name}.pdf")
    fig.savefig(out / f"{args.name}.png", dpi=300)
    written = [out / f"{args.name}.pdf", out / f"{args.name}.png"]
    if args.preview:
        images = VERIFICATION_DIR / "images"
        images.mkdir(exist_ok=True)
        fig.savefig(images / f"{args.name}.png", dpi=150)
        written.append(images / f"{args.name}.png")
    for path in written:
        print(path.relative_to(VERIFICATION_DIR))
    if counts:
        print(f"fields at t = {counts[2]}: {counts[0]} penalised cells, "
              f"{counts[1]} cut-link fluid cells")


if __name__ == "__main__":
    main()
