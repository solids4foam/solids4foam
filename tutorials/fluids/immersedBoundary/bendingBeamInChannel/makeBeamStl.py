#!/usr/bin/env python3
# Closed surface of the beam: 0.2 m wide (11.9 < x < 12.1), from below the
# bottom wall (y = -0.1 m) to its tip (y = 2.0295 m), extending beyond the mesh
# in the z direction, with 80 segments along its height so that it can bend
import math, os
x0, x1 = 11.9, 12.1
y0, y1 = -0.1, 2.0295
z0, z1 = -0.1, 0.2
ny = 80
ys = [y0 + (y1 - y0)*i/ny for i in range(ny + 1)]
tris = []
def quad(a, b, c, d):
    tris.append((a, b, c)); tris.append((a, c, d))
for i in range(ny):
    ya, yb = ys[i], ys[i + 1]
    quad((x0, ya, z0), (x0, ya, z1), (x0, yb, z1), (x0, yb, z0))  # x = x0, -x
    quad((x1, ya, z0), (x1, yb, z0), (x1, yb, z1), (x1, ya, z1))  # x = x1, +x
    quad((x0, ya, z0), (x0, yb, z0), (x1, yb, z0), (x1, ya, z0))  # z = z0, -z
    quad((x0, ya, z1), (x1, ya, z1), (x1, yb, z1), (x0, yb, z1))  # z = z1, +z
quad((x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1))      # bottom, -y
quad((x0, y1, z0), (x0, y1, z1), (x1, y1, z1), (x1, y1, z0))      # tip, +y
os.makedirs('constant/triSurface', exist_ok=True)
with open('constant/triSurface/beam.stl', 'w') as f:
    f.write('solid beam\n')
    for a, b, c in tris:
        u = [b[k] - a[k] for k in range(3)]; v = [c[k] - a[k] for k in range(3)]
        n = (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])
        m = math.sqrt(sum(x*x for x in n))
        f.write('  facet normal %g %g %g\n    outer loop\n' % tuple(x/m for x in n))
        for p in (a, b, c):
            f.write('      vertex %.8f %.8f %.8f\n' % p)
        f.write('    endloop\n  endfacet\n')
    f.write('endsolid beam\n')
