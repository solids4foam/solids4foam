#!/usr/bin/env python3
"""Write the load p = p0 cos(phi) into 0/finite-area/p of the annular plate.

Run from the case directory after makeFaMesh. The finite area mesh has no
field-expression utility, so the face centres are read from the ascii mesh
files and the load is written as a nonuniform list in the order of the finite
area faces.
"""

import math
import re
import sys
from pathlib import Path

P0 = 1000.0

CASE = Path(".")
TEMPLATE = CASE / "0" / "finite-area" / "p.uniform"
OUTPUT = CASE / "0" / "finite-area" / "p"


def list_body(path: Path) -> str:
    text = path.read_text()
    text = text[text.index("}") + 1:]
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//.*", "", text)
    match = re.search(r"(\d+)\s*\((.*)\)", text, re.DOTALL)
    if not match:
        raise SystemExit(f"no list found in {path}")
    return match.group(2)


def main() -> None:
    poly = CASE / "constant" / "polyMesh"
    if not (poly / "points").is_file():
        raise SystemExit("run blockMesh and makeFaMesh before writeCosineLoad.py")
    points = [
        tuple(float(v) for v in m.split())
        for m in re.findall(r"\(([^()]*)\)", list_body(poly / "points"))
    ]
    faces = [
        [int(v) for v in m.split()]
        for m in re.findall(r"\d+\s*\(([^()]*)\)", list_body(poly / "faces"))
    ]
    labels = [
        int(v)
        for v in list_body(CASE / "constant" / "finite-area" / "faMesh" / "faceLabels").split()
    ]

    values = []
    for label in labels:
        face = faces[label]
        cx = sum(points[i][0] for i in face) / len(face)
        cy = sum(points[i][1] for i in face) / len(face)
        values.append(P0 * math.cos(math.atan2(cy, cx)))

    text = TEMPLATE.read_text()
    entry = "internalField   nonuniform List<scalar>\n%d\n(\n%s\n)\n;" % (
        len(values),
        "\n".join("%.10g" % v for v in values),
    )
    text, count = re.subn(r"internalField\s+uniform\s+[-+0-9.eE]+\s*;", entry, text)
    if count != 1:
        raise SystemExit(f"no uniform internalField found in {TEMPLATE}")
    if OUTPUT.is_symlink() or OUTPUT.exists():
        OUTPUT.unlink()
    OUTPUT.write_text(text)
    print(f"Wrote p = p0 cos(phi) on {len(values)} faces to {OUTPUT}")


if __name__ == "__main__":
    main()
