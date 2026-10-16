#!/usr/bin/env python3
# Synthetic heart valves as closed, thick surfaces (binary STL):
#   valveTube.stl    three leaflets on a cylinder of radius R from the annulus
#                    (z = 0) to z = L, for the valveSliceAxis motion
#   annulus.stl      the ring of base vertices of valveTube.stl
#   valveClosed.stl  three flat leaflets closing the orifice at z = 0, and
#   valveOpen.stl    the same leaflets rotated about their hinges, with the
#                    same vertices, for the valveMorph motion
#   plate.stl        an annular plate below the annulus, from the radius
#                    R - t/2 to 2R, which closes the duct around the valve
import math, os, struct, sys

R = 0.01          # annulus radius [m]
L = 0.016         # leaflet length [m]
t = 0.0015        # leaflet thickness [m]
gap = 0.03        # angular gap between the leaflets [rad]
nTheta, nAxial = 24, 16


def thick_sheet(mid, nu, nv, thickness):
    """Closed surface of thickness 'thickness' about the mid-surface
    mid(i, j), i in [0, nu], j in [0, nv]; returns (points, triangles)."""
    P = [[mid(i, j) for j in range(nv + 1)] for i in range(nu + 1)]

    def sub(a, b): return [a[k] - b[k] for k in range(3)]
    def cross(a, b): return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
                             a[0]*b[1] - a[1]*b[0]]
    def unit(a):
        m = math.sqrt(sum(x*x for x in a)); return [x/m for x in a]

    def normal(i, j):
        du = sub(P[min(i + 1, nu)][j], P[max(i - 1, 0)][j])
        dv = sub(P[i][min(j + 1, nv)], P[i][max(j - 1, 0)])
        return unit(cross(du, dv))

    pts, idx = [], {}
    for side, s in ((0, 0.5), (1, -0.5)):
        for i in range(nu + 1):
            for j in range(nv + 1):
                n = normal(i, j)
                idx[side, i, j] = len(pts)
                pts.append([P[i][j][k] + s*thickness*n[k] for k in range(3)])

    tris = []
    def quad(a, b, c, d): tris.extend([(a, b, c), (a, c, d)])
    for i in range(nu):
        for j in range(nv):
            quad(idx[0, i, j], idx[0, i + 1, j], idx[0, i + 1, j + 1],
                 idx[0, i, j + 1])
            quad(idx[1, i, j], idx[1, i, j + 1], idx[1, i + 1, j + 1],
                 idx[1, i + 1, j])
    for i in range(nu):
        quad(idx[0, i, 0], idx[1, i, 0], idx[1, i + 1, 0], idx[0, i + 1, 0])
        quad(idx[0, i, nv], idx[0, i + 1, nv], idx[1, i + 1, nv],
             idx[1, i, nv])
    for j in range(nv):
        quad(idx[0, 0, j], idx[0, 0, j + 1], idx[1, 0, j + 1], idx[1, 0, j])
        quad(idx[0, nu, j], idx[1, nu, j], idx[1, nu, j + 1],
             idx[0, nu, j + 1])

    # Orient the normals outwards: positive enclosed volume
    vol = sum(
        sum(pts[a][k]*cross(pts[b], pts[c])[k] for k in range(3))
        for a, b, c in tris)/6
    if vol < 0:
        tris = [(a, c, b) for a, b, c in tris]
    return pts, tris


def write_stl(name, surfaces):
    """Binary STL of the triangles of the surfaces"""
    tris = []
    for pts, faces in surfaces:
        for a, b, c in faces:
            A, B, C = pts[a], pts[b], pts[c]
            u = [B[k] - A[k] for k in range(3)]
            v = [C[k] - A[k] for k in range(3)]
            n = [u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2],
                 u[0]*v[1] - u[1]*v[0]]
            m = math.sqrt(sum(x*x for x in n)) or 1
            tris.append([x/m for x in n] + A + B + C)
    with open(name, 'wb') as f:
        f.write(os.path.basename(name).encode().ljust(80, b' '))
        f.write(struct.pack('<I', len(tris)))
        for t in tris:
            f.write(struct.pack('<12fH', *t, 0))

def sector(k):
    th0 = 2*math.pi*k/3 + 0.5*gap
    return th0, th0 + 2*math.pi/3 - gap


out = sys.argv[1] if len(sys.argv) > 1 else 'constant/triSurface'
os.makedirs(out, exist_ok=True)

# Tube leaflets: mid-surface on the cylinder of radius R
tube = []
for k in range(3):
    th0, th1 = sector(k)
    def mid(i, j, th0=th0, th1=th1):
        th = th0 + (th1 - th0)*i/nTheta
        return [R*math.cos(th), R*math.sin(th), L*j/nAxial]
    tube.append(thick_sheet(mid, nTheta, nAxial, t))
write_stl(os.path.join(out, 'valveTube.stl'), tube)

# Annulus: a thin ring through the base vertices (z = 0) of the tube
ring = []
for k in range(3):
    th0, th1 = sector(k)
    def mid(i, j, th0=th0, th1=th1):
        th = th0 + (th1 - th0)*i/nTheta
        return [R*math.cos(th), R*math.sin(th), 0.25*t*(j - 0.5)]
    ring.append(thick_sheet(mid, nTheta, 1, t))
write_stl(os.path.join(out, 'annulus.stl'), ring)


# Hinged leaflets: radial lines from the hinge at radius R towards the axis,
# rotated by beta about the tangent of the hinge (beta = 0: flat, closed)
def hinged(beta):
    surfaces = []
    for k in range(3):
        th0, th1 = sector(k)
        def mid(i, j, th0=th0, th1=th1):
            th = th0 + (th1 - th0)*i/nTheta
            s = 0.95*R*j/nAxial
            r = R - s*math.cos(beta)
            return [r*math.cos(th), r*math.sin(th), s*math.sin(beta)]
        surfaces.append(thick_sheet(mid, nTheta, nAxial, t))
    return surfaces


write_stl(os.path.join(out, 'valveClosed.stl'), hinged(0))
write_stl(os.path.join(out, 'valveOpen.stl'), hinged(math.radians(80)))


# Annular plate below the annulus: rIn < r < rOut, -tp < z < 0
def annular_slab(rIn, rOut, z0, z1, n):
    pts, tris = [], []
    for r, z in ((rIn, z0), (rIn, z1), (rOut, z0), (rOut, z1)):
        for i in range(n):
            th = 2*math.pi*i/n
            pts.append([r*math.cos(th), r*math.sin(th), z])
    ib, it, ob, ot = 0, n, 2*n, 3*n
    def quad(a, b, c, d): tris.extend([(a, b, c), (a, c, d)])
    for i in range(n):
        j = (i + 1) % n
        quad(ot + i, it + i, it + j, ot + j)   # top
        quad(ob + i, ob + j, ib + j, ib + i)   # bottom
        quad(ob + i, ot + i, ot + j, ob + j)   # outer
        quad(ib + i, ib + j, it + j, it + i)   # inner
    return pts, [(a, c, b) for a, b, c in tris]  # outward normals


write_stl(os.path.join(out, 'plate.stl'),
          [annular_slab(R - 0.5*t, 2*R, -0.002, 0, 96)])
