#!/usr/bin/env python3
# Closed surface of the rigid cylinder of radius 0.05 m centred at (0.2, 0.2),
# extending beyond the mesh in the z direction, with 128 segments
import math, os
cx, cy, r = 0.2, 0.2, 0.05
z0, z1 = -0.015, 0.03
n = 128
bot = [(cx + r*math.cos(2*math.pi*i/n), cy + r*math.sin(2*math.pi*i/n), z0)
       for i in range(n)]
top = [(p[0], p[1], z1) for p in bot]
cb, ct = (cx, cy, z0), (cx, cy, z1)
tris = []
for i in range(n):
    j = (i + 1) % n
    tris.append((bot[i], bot[j], top[j]))
    tris.append((bot[i], top[j], top[i]))
    tris.append((cb, bot[j], bot[i]))
    tris.append((ct, top[i], top[j]))
os.makedirs('constant/fluid/triSurface', exist_ok=True)
with open('constant/fluid/triSurface/cylinder.stl', 'w') as f:
    f.write('solid cylinder\n')
    for a, b, c in tris:
        u = [b[k] - a[k] for k in range(3)]; v = [c[k] - a[k] for k in range(3)]
        nn = (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])
        m = math.sqrt(sum(x*x for x in nn))
        f.write('  facet normal %.9e %.9e %.9e\n    outer loop\n' % tuple(x/m for x in nn))
        for p in (a, b, c):
            f.write('      vertex %.9e %.9e %.9e\n' % p)
        f.write('    endloop\n  endfacet\n')
    f.write('endsolid cylinder\n')
