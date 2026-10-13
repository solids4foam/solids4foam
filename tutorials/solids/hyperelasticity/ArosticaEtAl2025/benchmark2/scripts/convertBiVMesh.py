#!/usr/bin/env python3
"""
Convert the biventricular mesh and fibres distributed with the Arostica et al.
(2025) benchmark (cardiac-benchmark-toolkit, data/biv_mesh_and_fibers_*_p1)
to an OpenFOAM polyMesh and cell-centred f0, s0, n0 fields.

Inputs (from the toolkit):
    bi_ventricular.h5      tetrahedra, nodes and facet markers
                           (10 base, 20 endo RV, 30 endo LV, 40 epi)
    fiber_directions.h5    nodal fibre, sheet and sheet-normal vectors (P1)
Outputs:
    constant/polyMesh/{points,faces,owner,neighbour,boundary}
    0.orig/f0, 0.orig/s0, 0.orig/n0   (cell and face values)

The nodal vectors are averaged to cells and faces and re-orthonormalised
(f0 first, then s0 by Gram-Schmidt, n0 = f0 x s0), which is what a P1
interpolation evaluated at the centroid gives, up to normalisation. The sheet
s0 is the toolkit's "sheet" vector, as used by the reference implementation.

Usage:
    python3 convertBiVMesh.py <toolkit>/data/biv_mesh_and_fibers_coarse_p1 [caseDir]
Needs numpy and h5py.
"""
import os
import sys

import h5py
import numpy as np

MARKERS = {10: "base", 20: "endo_rv", 30: "endo_lv", 40: "epi"}
ORDER = ["base", "endo_lv", "endo_rv", "epi"]

HEADER = r"""/*--------------------------------*- C++ -*----------------------------------*\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       %s;
    location    "%s";
    object      %s;
}
// Converted by scripts/convertBiVMesh.py from the cardiac-benchmark-toolkit
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def orthonormal(f, s):
    f = unit(f)
    s = unit(s - np.sum(s * f, axis=-1, keepdims=True) * f)
    return f, s, unit(np.cross(f, s))


def main():
    src = sys.argv[1]
    case = sys.argv[2] if len(sys.argv) > 2 else "."
    with h5py.File(os.path.join(src, "bi_ventricular.h5"), "r") as h:
        X = h["Mesh/mesh/geometry"][:]
        T = h["Mesh/mesh/topology"][:].astype(np.int64)
        tri = h["MeshFunction/0/mesh/topology"][:].astype(np.int64)
        val = h["MeshFunction/0/values"][:].ravel().astype(int)
    with h5py.File(os.path.join(src, "fiber_directions.h5"), "r") as h:
        if np.abs(h["Mesh/0/mesh/geometry"][:] - X).max() > 1e-12:
            sys.exit("Fibre and mesh node orderings differ")
        fN = h["VisualisationVector/0"][:]
        sN = h["VisualisationVector/1"][:]

    # Right-handed tetrahedra
    vol6 = np.einsum("ij,ij->i", np.cross(X[T[:, 1]] - X[T[:, 0]],
                                          X[T[:, 2]] - X[T[:, 0]]),
                     X[T[:, 3]] - X[T[:, 0]])
    T[vol6 < 0] = T[vol6 < 0][:, [0, 2, 1, 3]]
    if (np.abs(vol6) < 1e-30).any():
        sys.exit("Degenerate tetrahedra")

    # Outward faces of a right-handed tet
    TF = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
    faceMap = {}
    for ci, t in enumerate(T):
        for f in TF:
            verts = [int(t[k]) for k in f]
            key = tuple(sorted(verts))
            e = faceMap.get(key)
            if e is None:
                faceMap[key] = [ci, -1, verts]
            else:
                e[1] = ci

    marker = {tuple(sorted(map(int, tr))): MARKERS[v]
              for tr, v in zip(tri, val) if v in MARKERS}

    internal, patches = [], {k: [] for k in ORDER}
    for key, (o, n, verts) in faceMap.items():
        if n >= 0:
            internal.append((o, n, verts))
        else:
            name = marker.get(key)
            if name is None:
                sys.exit(f"Boundary face {key} has no marker")
            patches[name].append((o, verts))
    internal.sort(key=lambda f: (f[0], f[1]))
    for k in patches:
        patches[k].sort(key=lambda f: f[0])

    # Check outward orientation
    cc = X[T].mean(1)
    for k in ORDER:
        for o, v in patches[k]:
            P = X[v]
            if np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P.mean(0) - cc[o]) <= 0:
                sys.exit("Inward boundary face")

    poly = os.path.join(case, "constant/polyMesh")
    os.makedirs(poly, exist_ok=True)
    faces = [v for _, _, v in internal] + [v for k in ORDER for _, v in patches[k]]
    owner = [o for o, _, _ in internal] + [o for k in ORDER for o, _ in patches[k]]
    neigh = [n for _, n, _ in internal]

    def wr(name, cls, body, loc="constant/polyMesh"):
        with open(os.path.join(case, loc, name), "w") as fh:
            fh.write(HEADER % (cls, loc, name))
            fh.write(body)

    wr("points", "vectorField", f"{len(X)}\n(\n" + "".join(
        f"({p[0]:.12g} {p[1]:.12g} {p[2]:.12g})\n" for p in X) + ")\n")
    wr("faces", "faceList", f"{len(faces)}\n(\n" + "".join(
        f"3({f[0]} {f[1]} {f[2]})\n" for f in faces) + ")\n")
    wr("owner", "labelList", f"{len(owner)}\n(\n" + "".join(
        f"{x}\n" for x in owner) + ")\n")
    wr("neighbour", "labelList", f"{len(neigh)}\n(\n" + "".join(
        f"{x}\n" for x in neigh) + ")\n")
    b = f"{len(ORDER)}\n(\n"
    start = len(internal)
    for k in ORDER:
        b += (f"    {k}\n    {{\n        type            patch;\n"
              f"        nFaces          {len(patches[k])};\n"
              f"        startFace       {start};\n    }}\n")
        start += len(patches[k])
    wr("boundary", "polyBoundaryMesh", b + ")\n")

    # Fibres at cells and boundary faces
    fC, sC, nC = orthonormal(fN[T].mean(1), sN[T].mean(1))
    os.makedirs(os.path.join(case, "0.orig"), exist_ok=True)
    for name, vc, vN in (("f0", fC, fN), ("s0", sC, sN), ("n0", nC, None)):
        body = "dimensions      [0 0 0 0 0 0 0];\n\ninternalField   nonuniform List<vector>\n"
        body += f"{len(vc)}\n(\n" + "".join(
            f"({v[0]:.8g} {v[1]:.8g} {v[2]:.8g})\n" for v in vc) + ")\n;\n\nboundaryField\n{\n"
        for k in ORDER:
            idx = [v for _, v in patches[k]]
            fF, sF, nF = orthonormal(fN[idx].mean(1), sN[idx].mean(1))
            vf = {"f0": fF, "s0": sF, "n0": nF}[name]
            body += f"    {k}\n    {{\n        type            calculated;\n        value           nonuniform List<vector>\n"
            body += f"{len(vf)}\n(\n" + "".join(
                f"({v[0]:.8g} {v[1]:.8g} {v[2]:.8g})\n" for v in vf) + ")\n;\n    }\n"
        body += "}\n"
        wr(name, "volVectorField", body, loc="0.orig")

    print(f"{len(T)} tetrahedra, {len(X)} points, {len(internal)} internal faces")
    print("patches: " + ", ".join(f"{k} {len(patches[k])}" for k in ORDER))
    print(f"volume {np.abs(vol6).sum()/6*1e6:.2f} ml")


if __name__ == "__main__":
    main()
