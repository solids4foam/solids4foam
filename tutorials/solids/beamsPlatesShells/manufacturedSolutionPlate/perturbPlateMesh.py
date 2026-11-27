#!/usr/bin/env python3
"""Perturb the interior points of the one-layer plate mesh in its plane.

Every interior point is moved by a uniformly distributed random fraction of
the local mesh spacing, by the same vector on the front and back surfaces so
that the cells remain straight prisms, while the points on the plate edges
stay in place. The fraction is given on the command line (default 0.15, as a
fraction of the shortest in-plane edge at the point) together with the seed
of the random generator, so the perturbed mesh is reproducible and the
distortion scales with the mesh spacing when the mesh is refined. Run from
the case directory after blockMesh; the points in constant/polyMesh are
overwritten.

perturbMeshPoints is not used because it perturbs the two surfaces of the
one-layer mesh independently and the twisted side faces then fail its mesh
checks.
"""

import random
import re
import sys
from pathlib import Path

POINTS = Path("constant") / "polyMesh" / "points"
BOUNDARY = Path("constant") / "polyMesh" / "boundary"
FACES = Path("constant") / "polyMesh" / "faces"


def list_body(path: Path) -> tuple[str, str, str]:
    """Split the file into text before the list, the list body and the rest."""
    text = path.read_text()
    header_end = text.index("}") + 1
    body = re.sub(r"/\*.*?\*/", "", text[header_end:], flags=re.DOTALL)
    body = re.sub(r"//.*", "", body)
    match = re.search(r"(\d+)\s*\(", body)
    start = header_end + match.start()
    open_paren = text.index("(", start)
    close_paren = text.rindex(")")
    return text[:open_paren + 1], text[open_paren + 1:close_paren], text[close_paren:]


def main() -> None:
    fraction = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    if not POINTS.is_file():
        raise SystemExit("run blockMesh before perturbPlateMesh.py")

    head, body, tail = list_body(POINTS)
    points = [
        [float(v) for v in m.split()] for m in re.findall(r"\(([^()]*)\)", body)
    ]
    _, faces_body, _ = list_body(FACES)
    faces = [
        [int(v) for v in m.split()]
        for m in re.findall(r"\d+\s*\(([^()]*)\)", faces_body)
    ]

    # Points on the side patches of the plate (every patch except the two
    # plate surfaces, recognised by their large number of faces) stay fixed
    boundary = re.sub(r"//.*", "", BOUNDARY.read_text())
    patches = re.findall(
        r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", boundary, re.DOTALL
    )
    sizes = sorted(int(n) for _, n, _ in patches)
    fixed = set()
    for _, n_faces, start in patches:
        if int(n_faces) >= sizes[-2]:
            continue
        for f in range(int(start), int(start) + int(n_faces)):
            fixed.update(faces[f])

    # Shortest in-plane edge at each point, from the faces
    spacing = [float("inf")] * len(points)
    for face in faces:
        for a, b in zip(face, face[1:] + face[:1]):
            pa, pb = points[a], points[b]
            length = ((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2) ** 0.5
            if length > 1e-12:
                spacing[a] = min(spacing[a], length)
                spacing[b] = min(spacing[b], length)

    # One displacement per in-plane position, shared by the two surfaces
    rng = random.Random(seed)
    displacement: dict[tuple[float, float], tuple[float, float]] = {}
    moved = 0
    for i, (x, y, z) in enumerate(points):
        if i in fixed:
            continue
        key = (round(x, 9), round(y, 9))
        if key not in displacement:
            displacement[key] = (
                fraction * spacing[i] * rng.uniform(-1.0, 1.0),
                fraction * spacing[i] * rng.uniform(-1.0, 1.0),
            )
        dx, dy = displacement[key]
        points[i] = [x + dx, y + dy, z]
        moved += 1

    POINTS.write_text(
        head
        + "\n"
        + "\n".join("(%.12g %.12g %.12g)" % tuple(p) for p in points)
        + "\n"
        + tail
    )
    print(f"Moved {moved} points by up to {fraction:g} of the local spacing (seed {seed})")


if __name__ == "__main__":
    main()
