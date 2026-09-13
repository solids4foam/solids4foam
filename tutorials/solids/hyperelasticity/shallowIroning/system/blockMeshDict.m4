/*--------------------------------*- C++ -*----------------------------------*\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     v2.4                                                           |
| Web:         https://solids4foam.github.io                                  |
| Disclaimer:  This offering is not approved or endorsed by OpenCFD Limited,  |
|              producer and distributor of the OpenFOAM software via          |
|              www.openfoam.com, and owner of the OPENFOAM® and OpenCFD®      |
|              trade marks.                                                   |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// M4
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'printf ($1)')])
m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

// GEOMETRY
// Taken from Yastrebov book "Numerical Methods in Contact Mechanics"
m4_define(h1, 0.95)
m4_define(h2, 4)
m4_define(a1, 0.3)
m4_define(a2, 0.2)
m4_define(a3, 0.9)
m4_define(d1, 0.2)
m4_define(d2, 1.2)
m4_define(d3, 10.6)
m4_define(l, 1)
m4_define(r, 0.75)

// MESH
m4_define(foundationGrid, 120 30 1) 
m4_define(slabGrid, 12 6 1) 
m4_define(grading, 1 1 1)

// FRONT AND BACK PLANES
m4_define(zA, 0)
m4_define(zB, l)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
// Foundation: lower rectangular body ("foundation" cell zone)

    //Plane A:
    (0 0 zA) vlabel(A0)
    (calc(d1+d2+d3) 0 zA) vlabel(A1)
    (calc(d1+d2+d3) h2 zA) vlabel(A2)
    (0 h2 zA) vlabel(A3)

    //Plane B:
    (0 0 zB) vlabel(B0)
    (calc(d1+d2+d3) 0 zB) vlabel(B1)
    (calc(d1+d2+d3) h2 zB) vlabel(B2)
    (0 h2 zB) vlabel(B3)

// Block: upper rounded body ("slab" cell zone)

    //Plane A:
    (d1 calc(h2+a2+a1) zA) vlabel(C0)
    (calc(d1+d2) calc(h2+a2+a1) zA) vlabel(C1)
    (calc(d1+d2) calc(h2+a2+a3+a1) zA) vlabel(C2)
    (d1 calc(h2+a2+a3+a1) zA) vlabel(C3)

    //Plane B:
    (d1 calc(h2+a2+a1) zB) vlabel(D0)
    (calc(d1+d2) calc(h2+a2+a1) zB) vlabel(D1)
    (calc(d1+d2) calc(h2+a2+a3+a1) zB) vlabel(D2)
    (d1 calc(h2+a2+a3+a1) zB) vlabel(D3)
);

blocks
(
    hex ( A0 A1 A2 A3 B0 B1 B2 B3 ) foundation (foundationGrid) simpleGrading (grading)
    hex ( C0 C1 C2 C3 D0 D1 D2 D3 ) slab (slabGrid) simpleGrading (grading)
);

edges
(
    // Plane A
    arc  D0 D1  (calc(d1+d2*0.5) calc(h2+a2) zB)
    // Plane B
    arc  C0 C1  (calc(d1+d2*0.5) calc(h2+a2) zA)
);

boundary
(
    frontAndBack
    {
        type empty;
        faces
        (
            (A0 A1 A2 A3)
            (B1 B2 B3 B0)
            (D0 D1 D2 D3)
            (C1 C0 C3 C2)
        );
    }
    foundationContact
    {
        type patch;
        faces
        (
            (A3 B3 B2 A2)
        );
    }
    slabContact
    {
        type patch;
        faces
        (
            (D0 C0 C1 D1)
            (D1 C1 C2 D2)
        );
    }
    fixed
    {
        type patch;
        faces
        (
            (B0 A0 A1 B1)
        );
    }
    displacement
    {
        type patch;
        faces
        (
            (D3 D2 C2 C3)
        );
    }
    zeroTraction
    {
        type patch;
        faces
        (
            (B0 B3 A3 A0)
            (B2 B1 A1 A2)
            (C0 D0 D3 C3)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
