#!/usr/bin/env python3
"""Write the manufactured load and the exact fields of the plate case.

The manufactured deflection is

    w = w0 [a (x/L)^2 + b (x/L)^3 + (x/L)^4] cos(pi y / W),

with a and b chosen so that the bending moment and the Kirchhoff shear force
vanish on the free edge x = L; w and its slope vanish on the clamped edge
x = 0, and the cosine makes y = 0 and y = W symmetry lines. The load
p = D laplacian(laplacian(w)) that produces this deflection is written to
0/finite-area/p from the template 0/finite-area/p.orig, and the exact
deflection and moment sum are written to 0/finite-area/wExact and MExact for
comparison. Run from the case directory after makeFaMesh.
"""

import math
import re
from pathlib import Path

# Plate dimensions and material, as in constant/ and system/blockMeshDict
L = 10.0
W = 10.0
H = 0.1
E = 2e11
NU = 0.3
D = E * H**3 / (12.0 * (1.0 - NU**2))

# Deflection amplitude of the manufactured solution
W0 = 0.05

CASE = Path(".")
FA = CASE / "0" / "finite-area"
TEMPLATE = FA / "p.orig"


class Manufactured:
    """w = w0 f(x) cos(k y) with f a quartic satisfying the edge conditions."""

    def __init__(self, length: float, width: float, nu: float, w0: float):
        self.L, self.k, self.nu = length, math.pi / width, nu
        k, lam = self.k, length
        # f = a x^2 + b x^3 + c x^4 with c = 1 / L^4, and
        #   f'' - nu k^2 f = 0 and f''' - (2 - nu) k^2 f' = 0 at x = L
        c = 1.0 / lam**4
        a11 = 2.0 - nu * k * k * lam * lam
        a12 = 6.0 * lam - nu * k * k * lam**3
        r1 = -(12.0 * lam * lam - nu * k * k * lam**4) * c
        a21 = -(2.0 - nu) * k * k * 2.0 * lam
        a22 = 6.0 - (2.0 - nu) * k * k * 3.0 * lam * lam
        r2 = -(24.0 * lam - (2.0 - nu) * k * k * 4.0 * lam**3) * c
        det = a11 * a22 - a12 * a21
        a = (r1 * a22 - a12 * r2) / det
        b = (a11 * r2 - a21 * r1) / det
        fmax = max(abs(a * x * x + b * x**3 + c * x**4)
                   for x in (i * lam / 2000.0 for i in range(2001)))
        scale = w0 / fmax
        self.a, self.b, self.c = scale * a, scale * b, scale * c

    def f(self, x: float) -> float:
        return self.a * x * x + self.b * x**3 + self.c * x**4

    def df(self, x: float) -> float:
        return 2.0 * self.a * x + 3.0 * self.b * x * x + 4.0 * self.c * x**3

    def d2f(self, x: float) -> float:
        return 2.0 * self.a + 6.0 * self.b * x + 12.0 * self.c * x * x

    def d4f(self) -> float:
        return 24.0 * self.c

    def w(self, x: float, y: float) -> float:
        return self.f(x) * math.cos(self.k * y)

    def grad_w(self, x: float, y: float) -> tuple[float, float]:
        return (
            self.df(x) * math.cos(self.k * y),
            -self.k * self.f(x) * math.sin(self.k * y),
        )

    def laplacian_w(self, x: float, y: float) -> float:
        return (self.d2f(x) - self.k * self.k * self.f(x)) * math.cos(self.k * y)

    def M(self, x: float, y: float, D: float) -> float:
        return -D * self.laplacian_w(x, y)

    def p(self, x: float, y: float, D: float) -> float:
        k2 = self.k * self.k
        return D * (self.d4f() - 2.0 * k2 * self.d2f(x) + k2 * k2 * self.f(x)) * math.cos(self.k * y)


def list_body(path: Path) -> str:
    text = path.read_text()
    text = text[text.index("}") + 1:]
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//.*", "", text)
    match = re.search(r"(\d+)\s*\((.*)\)", text, re.DOTALL)
    if not match:
        raise SystemExit(f"no list found in {path}")
    return match.group(2)


def face_centre(points, face) -> tuple[float, float]:
    """Area-weighted centre of a polygon, as OpenFOAM computes it."""
    n = len(face)
    px = sum(points[i][0] for i in face) / n
    py = sum(points[i][1] for i in face) / n
    if n == 3:
        return px, py
    sum_a = sum_ax = sum_ay = 0.0
    for i in range(n):
        x1, y1 = points[face[i]][0], points[face[i]][1]
        x2, y2 = points[face[(i + 1) % n]][0], points[face[(i + 1) % n]][1]
        area = abs((x1 - px) * (y2 - py) - (x2 - px) * (y1 - py))
        sum_a += area
        sum_ax += area * (x1 + x2 + px) / 3.0
        sum_ay += area * (y1 + y2 + py) / 3.0
    return sum_ax / sum_a, sum_ay / sum_a


def write_field(template: str, name: str, cls_object: str, dimensions: str, values) -> None:
    text = template.replace("object      p;", f"object      {cls_object};")
    text = re.sub(r"dimensions\s+\[[^\]]*\];", f"dimensions      {dimensions};", text)
    text = re.sub(r"//\s*Template for the manufactured load.*?\n\n", "", text, flags=re.DOTALL)
    entry = "internalField   nonuniform List<scalar>\n%d\n(\n%s\n)\n;" % (
        len(values),
        "\n".join("%.10g" % v for v in values),
    )
    text, count = re.subn(r"internalField\s+uniform\s+[-+0-9.eE]+\s*;", entry, text)
    if count != 1:
        raise SystemExit(f"no uniform internalField found in {TEMPLATE}")
    (FA / name).write_text(text)
    print(f"Wrote {name} on {len(values)} faces")


def find_mesh_file(*candidates: Path) -> Path:
    """A mesh written to the start time by perturbMeshPoints, and the finite
    area mesh built from it, take precedence over the ones in constant."""
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit(
        "create the mesh and run makeFaMesh before writeManufacturedSolution.py"
    )


def main() -> None:
    poly = find_mesh_file(
        CASE / "0" / "polyMesh" / "points", CASE / "constant" / "polyMesh" / "points"
    ).parent
    labels_file = find_mesh_file(
        CASE / "0" / "finite-area" / "faMesh" / "faceLabels",
        CASE / "constant" / "finite-area" / "faMesh" / "faceLabels",
    )
    points = [
        tuple(float(v) for v in m.split())
        for m in re.findall(r"\(([^()]*)\)", list_body(poly / "points"))
    ]
    faces = [
        [int(v) for v in m.split()]
        for m in re.findall(r"\d+\s*\(([^()]*)\)", list_body(poly / "faces"))
    ]
    labels = [int(v) for v in list_body(labels_file).split()]
    centres = [face_centre(points, faces[label]) for label in labels]

    exact = Manufactured(L, W, NU, W0)
    template = TEMPLATE.read_text()
    write_field(template, "p", "p", "[1 -1 -2 0 0 0 0]", [exact.p(x, y, D) for x, y in centres])
    write_field(template, "wExact", "wExact", "[0 1 0 0 0 0 0]", [exact.w(x, y) for x, y in centres])
    write_field(template, "MExact", "MExact", "[1 1 -2 0 0 0 0]", [exact.M(x, y, D) for x, y in centres])


if __name__ == "__main__":
    main()
