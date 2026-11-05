#!/usr/bin/env python3
"""Closed STL box spanning the duct: back face at x = 0.303 m, thickness
RATIO*h (h = 0.01 m, the cell size); usage: makePlateStl.py RATIO > plate.stl"""
import sys
h = 0.01
x0, x1 = 0.303, 0.303 + float(sys.argv[1])*h
y0, y1, z0, z1 = -0.05, 0.25, -0.05, 0.06
v = [(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)]
quads = [(0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
print("solid plate")
for a, b, c, d in quads:
    for t in ((a, b, c), (a, c, d)):
        print("facet normal 0 0 0\n outer loop")
        for i in t:
            print("  vertex %.9g %.9g %.9g" % v[i])
        print(" endloop\nendfacet")
print("endsolid plate")
