import math, sys
h=0.01
def tri(s,a,b,c): return "facet normal 0 0 0\n outer loop\n"+"".join("  vertex %.9g %.9g %.9g\n"%p for p in (a,b,c))+" endloop\nendfacet\n"
def quad(a,b,c,d): return tri(0,a,b,c)+tri(0,a,c,d)
def box(x0,x1,y0,y1,z0,z1):
    v=[(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)]
    q=[(0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
    return "solid plate\n"+"".join(quad(*[v[i] for i in f]) for f in q)+"endsolid plate\n"
def shell(curve,t,y0,y1,z0,z1,n=400):
    """thin shell: the curve x=curve(y), thickness t measured along the normal; closed prism in z"""
    ys=[y0+(y1-y0)*i/n for i in range(n+1)]
    F=[];B=[]
    for y in ys:
        e=1e-6; dx=(curve(y+e)-curve(y-e))/(2*e)
        nx,ny=1.0,-dx; m=math.hypot(nx,ny); nx/=m; ny/=m   # unit normal, +x side
        x=curve(y)
        F.append((x+0.5*t*nx,y+0.5*t*ny)); B.append((x-0.5*t*nx,y-0.5*t*ny))
    s="solid plate\n"
    for i in range(n):
        f0,f1,b0,b1=F[i],F[i+1],B[i],B[i+1]
        # front (+x side) surface, normal +x: orientation with y increasing and z increasing
        P=lambda p,z:(p[0],p[1],z)
        s+=quad(P(f0,z0),P(f1,z0),P(f1,z1),P(f0,z1))*0+ ""
        # outward normal of front face: tangent t=(f1-f0), z-hat; n = t x z = (ty,-tx,0)... check numerically below
        s+=quad(P(f0,z0),P(f1,z0),P(f1,z1),P(f0,z1))
        s+=quad(P(b0,z0),P(b0,z1),P(b1,z1),P(b1,z0))
        s+=quad(P(f0,z1),P(f1,z1),P(b1,z1),P(b0,z1))   # top, +z
        s+=quad(P(f0,z0),P(b0,z0),P(b1,z0),P(f1,z0))   # bottom, -z
    s+=quad((F[0][0],F[0][1],z0),(F[0][0],F[0][1],z1),(B[0][0],B[0][1],z1),(B[0][0],B[0][1],z0))
    s+=quad((F[n][0],F[n][1],z0),(B[n][0],B[n][1],z0),(B[n][0],B[n][1],z1),(F[n][0],F[n][1],z1))
    return s+"endsolid plate\n"

ang=math.radians(30)
fn={"tilt30":lambda y:0.35+(y-0.1)*math.tan(ang),"arc":lambda y:0.35-2.0*(y-0.1)**2}
open(sys.argv[2],"w").write(shell(fn[sys.argv[1]],float(sys.argv[3])*h,-0.05,0.25,-0.05,0.06))
