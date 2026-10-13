#!/usr/bin/env python3
"""
Report J = det(F) over the written times: minimum, maximum and where the
minimum is. With kappa = 1e6 Pa and stresses of ~1e5 Pa, |J - 1| of a few
per cent is expected; a cell far below that marks where a run is losing
accuracy and, if it keeps falling, where it will invert.

Usage (case directory, after a run):  python3 scripts/checkJacobian.py [times...]
Needs numpy; reads the ascii constant/polyMesh and <time>/F.
"""
import glob
import os
import re
import sys

import numpy as np


def read_list(path, ncomp):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", s) \
        or re.search(r"\n(\d+)\s*\n\(", s)
    n = int(m.group(1))
    body = s[m.end():]
    nums = np.array(re.findall(r"[-+]?[0-9.]+(?:[eE][-+]?\d+)?",
                               body[:body.index("\n)")]), float)
    return nums[:n * ncomp].reshape(n, ncomp) if ncomp > 1 else nums[:n]


def cell_centres():
    pts = read_list("constant/polyMesh/points", 3)
    s = open("constant/polyMesh/faces").read()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    nf = int(m.group(1))
    faces = [list(map(int, f.split()))
             for f in re.findall(r"\d+\(([\d ]+)\)", s[m.end():])][:nf]
    own = read_list("constant/polyMesh/owner", 1).astype(int)
    nei = read_list("constant/polyMesh/neighbour", 1).astype(int)
    nc = own.max() + 1
    cc = np.zeros((nc, 3))
    cnt = np.zeros(nc)
    for i, f in enumerate(faces):
        c = pts[f].mean(0)
        cc[own[i]] += c
        cnt[own[i]] += 1
        if i < len(nei):
            cc[nei[i]] += c
            cnt[nei[i]] += 1
    return cc / cnt[:, None]


def main():
    times = sys.argv[1:] or sorted(
        (d for d in glob.glob("[0-9]*") if os.path.isfile(f"{d}/F")),
        key=float)
    cc = cell_centres()
    print(f"{'time':>8s} {'Jmin':>8s} {'Jmax':>8s}   Jmin at [mm]")
    for t in times:
        F = read_list(f"{t}/F", 9).reshape(-1, 3, 3)
        J = np.linalg.det(F)
        i = int(J.argmin())
        x = cc[i] * 1e3
        print(f"{t:>8s} {J.min():8.4f} {J.max():8.4f}   "
              f"({x[0]:7.1f} {x[1]:6.1f} {x[2]:6.1f})  cell {i}")


if __name__ == "__main__":
    main()
