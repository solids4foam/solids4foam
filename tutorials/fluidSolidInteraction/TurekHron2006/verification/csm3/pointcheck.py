#!/usr/bin/env python3
"""Compare the solidPointDisplacement value at tip A with the pointD value at
the mesh vertex (0.6, 0.2) of the written steady solution (time 1)."""
import re, sys
import numpy as np

V = r"\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)"

def points(case):
    t = open(f"{case}/constant/polyMesh/points").read()
    return np.array(re.findall(V, t[t.index("FoamFile"):].split("}", 1)[1]), float)

def pointD(case):
    t = open(f"{case}/1/pointD").read()
    j = t.index("internalField")
    return np.array(re.findall(V, t[j:t.index("boundaryField")]), float)

for case in sys.argv[1:]:
    p = points(case); d = pointD(case)
    i = np.argmin(np.linalg.norm(p[:, :2] - np.array([0.6, 0.2]), axis=1) + abs(p[:, 2]))
    fn = open(f"{case}/postProcessing/0/solidPointDisplacement_pointDisp.dat").read().split("\n")[-2].split()
    print(case.rstrip("/").split("/")[-1], "vertex", p[i], "pointD uy %.4f mm" % (1e3*d[i][1]),
          "| function uy %.4f mm" % (1e3*float(fn[2])),
          "| pointD ux %.5f function ux %.5f" % (1e3*d[i][0], 1e3*float(fn[1])))
