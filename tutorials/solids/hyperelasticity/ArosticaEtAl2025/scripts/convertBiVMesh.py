#!/usr/bin/env python3
"""
Convert the biventricular mesh and fibres distributed with the Arostica et al.
(2025) benchmark (cardiac-benchmark-toolkit, data/biv_mesh_and_fibers_*_p1)
to an OpenFOAM polyMesh and cell-centred f0, s0, n0 fields.

Inputs (from the toolkit, downloaded at a fixed commit and checked against
their SHA-256 sums when a mesh name is given):
    bi_ventricular.h5      tetrahedra, nodes and facet markers
                           (10 base, 20 endo RV, 30 endo LV, 40 epi)
    fiber_directions.h5    nodal fibre, sheet and sheet-normal vectors (P1)
Outputs:
    constant/polyMesh/{points,faces,owner,neighbour,boundary}
    0/f0, 0/s0, 0/n0   (cell and face values)

The nodal vectors are averaged to cells and faces and re-orthonormalised
(f0 first, then s0 by Gram-Schmidt, n0 = f0 x s0), which is what a P1
interpolation evaluated at the centroid gives, up to normalisation. The sheet
s0 is the toolkit's "sheet" vector, as used by the reference implementation.

Usage, from the case directory:
    python3 ../scripts/convertBiVMesh.py coarse|fine
        downloads the toolkit's files into data/<mesh> (once) and converts
    python3 ../scripts/convertBiVMesh.py <dir with the two .h5 files>
        converts a local copy, e.g. of a clone of the toolkit
Needs numpy and h5py.
"""
import hashlib
import os
import sys
import urllib.request

import h5py
import numpy as np

# The toolkit's data, at a fixed commit, with the SHA-256 sums of the files
TOOLKIT = ("https://raw.githubusercontent.com/Reidmen/cardiac_benchmark_toolkit/"
           "e8d47553cfc83eb274eba3e177de33148e7f441c/data")
DATA = {
    "coarse": {
        "bi_ventricular.h5":
            "0711b37423ab4f87b9b8fe68e02dc8261e8df023621a61202756df008c8ddd53",
        "fiber_directions.h5":
            "8859c8bc10ec32aca97069f2bf513b31c2cacbecdfae7dfbdd3bdd99684d1d69",
    },
    "fine": {
        "bi_ventricular.h5":
            "230e0f7c67b03066035e34fd69f8d3cd572addd9ad99f1538c3c9e11fc3ff0ef",
        "fiber_directions.h5":
            "8f9f1ba8491986580a8ced7964a5eaa1b821ef2aa1201a976ed1e1764569765d",
    },
}

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


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(mesh, case):
    """Download the toolkit's files of a mesh into data/<mesh>, once."""
    dest = os.path.join(case, "data", mesh)
    os.makedirs(dest, exist_ok=True)
    for name, digest in DATA[mesh].items():
        path = os.path.join(dest, name)
        if os.path.isfile(path) and sha256(path) == digest:
            continue
        url = f"{TOOLKIT}/biv_mesh_and_fibers_{mesh}_p1/{name}"
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, path)
        if sha256(path) != digest:
            os.remove(path)
            sys.exit(f"{name}: SHA-256 mismatch, download rejected")
    return dest


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    case = sys.argv[2] if len(sys.argv) > 2 else "."
    src = sys.argv[1]
    if src in DATA:
        src = download(src, case)
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
    os.makedirs(os.path.join(case, "0"), exist_ok=True)
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
        wr(name, "volVectorField", body, loc="0")

    print(f"{len(T)} tetrahedra, {len(X)} points, {len(internal)} internal faces")
    print("patches: " + ", ".join(f"{k} {len(patches[k])}" for k in ORDER))
    print(f"volume {np.abs(vol6).sum()/6*1e6:.2f} ml")


if __name__ == "__main__":
    main()
