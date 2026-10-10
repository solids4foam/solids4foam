#!/usr/bin/env python3
"""Smallest fan-triangle area on the FSI3 plate interface (ascii polyMesh).

The AMI point-transfer weights of solids4foam (before PR #546) used
triangle::pointToBarycentric on the fan triangles (two face vertices + face
centre) of the interface faces, which returns (1/3, 1/3, 1/3) when
4 A^2 < SMALL = 1e-15 m^4 (A < 1.58e-8 m^2). This prints, for the solid and
fluid plate patches of each run directory given, the smallest face and fan
triangle areas and 4 A^2 relative to that threshold.

    python3 plate_interface_faces.py work/iqnils_mesh_1x work/iqnils_mesh_4x
"""
import re, sys, os
def body(path):
    s = open(path).read()
    s = re.sub(r"//.*", "", s); s = s[s.index("}", s.index("FoamFile"))+1:]
    return s
def points(p):
    s = body(p); i = s.index("("); n = int(s[:i].split()[-1])
    return [tuple(map(float, m)) for m in re.findall(r"\(([^()]+) ([^()]+) ([^()]+)\)", s[i:])][:n]
def faces(p, start, n):
    s = body(p); i = s.index("("); f = re.findall(r"\d+\(([\d ]+)\)", s[i+1:])
    return [list(map(int, x.split())) for x in f[start:start+n]]
def patch(p, name):
    s = body(p); m = re.search(name + r"\s*\{[^}]*nFaces\s+(\d+);[^}]*startFace\s+(\d+);", s)
    return int(m.group(2)), int(m.group(1))
def sub(a,b): return [a[k]-b[k] for k in range(3)]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def mag(a): return sum(x*x for x in a)**0.5
for case in sys.argv[1:]:
    for region in ("solid","fluid"):
        d = os.path.join(case,"constant",region,"polyMesh")
        P = points(d+"/points"); st, n = patch(d+"/boundary","plate")
        tri = []; fa = []
        for f in faces(d+"/faces", st, n):
            pts = [P[k] for k in f]; c = [sum(q[k] for q in pts)/len(pts) for k in range(3)]
            A = 0
            for j in range(len(f)):
                t = 0.5*mag(cross(sub(pts[j],c), sub(pts[(j+1)%len(f)],c))); tri.append(t); A += t
            fa.append(A)
        print(f"{os.path.basename(case):28s} {region}: nFaces {n:5d}  min face {min(fa):.3e} m2  min fan tri {min(tri):.3e} m2  4A^2 {4*min(tri)**2:.2e} m4  ratio to 1e-15: {4*min(tri)**2/1e-15:.2e}")
