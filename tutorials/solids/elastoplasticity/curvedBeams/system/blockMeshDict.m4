/*--------------------------------*- C++ -*----------------------------------*\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     v2.0                                                           |
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
    location    "constant/polyMesh";
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

m4_define(R8, 8)
m4_define(R10, 10)
m4_define(R12, 12)
m4_define(L, 1)
m4_define(H, 17)

m4_define(xCenterUpperCyl, calc(-1*(((R10+R12)*(R10+R12)-H*H)**0.5-R12) ))

// MESH
m4_define(cylinderUpper, 5 50 1)
m4_define(cylinderLower, 5 50 1)
m4_define(grading, 1 1 1)

// FRONT AND BACK PLANES
m4_define(zA, 0)
m4_define(zB, L)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
// Slab

    //Plane A:
    (0             0 zA) vlabel(A0)
    (calc(R12-R10) 0 zA) vlabel(A1)
    (calc(R12+R10) 0 zA) vlabel(A2)
    (calc(2*R12)   0 zA) vlabel(A3)

    //Plane B:
    (0             0 zB) vlabel(B0)
    (calc(R12-R10) 0 zB) vlabel(B1)
    (calc(R12+R10) 0 zB) vlabel(B2)
    (calc(2*R12)   0 zB) vlabel(B3)

// Block

    //Plane C = Plane A:
    (calc(xCenterUpperCyl-R10) H zA) vlabel(C0)
    (calc(xCenterUpperCyl-R8)  H zA) vlabel(C1)
    (calc(xCenterUpperCyl+R8)  H zA) vlabel(C2)
    (calc(xCenterUpperCyl+R10) H zA) vlabel(C3)

    //Plane D = plane B:
    (calc(xCenterUpperCyl-R10) H zB) vlabel(D0)
    (calc(xCenterUpperCyl-R8)  H zB) vlabel(D1)
    (calc(xCenterUpperCyl+R8)  H zB) vlabel(D2)
    (calc(xCenterUpperCyl+R10) H zB) vlabel(D3)
);

blocks
(
    hex ( A0 A1 A2 A3 B0 B1 B2 B3 ) (cylinderLower) simpleGrading (grading)
    hex ( C3 C2 C1 C0 D3 D2 D1 D0 ) (cylinderUpper) simpleGrading (grading)
);

edges
(
    // Plane A
    arc  A0 A3 (R12 R12 zA)
    arc  A1 A2 (R12 R10 zA)
    arc  C1 C2 (xCenterUpperCyl calc(H-R8) zA)
    arc  C0 C3 (xCenterUpperCyl calc(H-R10) zA)
    // Plane B
    arc  B0 B3 (R12 R12 zB)
    arc  B1 B2 (R12 R10 zB)
    arc  D1 D2 (xCenterUpperCyl calc(H-R8) zB)
    arc  D0 D3 (xCenterUpperCyl calc(H-R10) zB)
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

    fixed
    {
        type patch;
        faces
        (
            (B0 A0 A1 B1)
            (A3 B3 B2 A2)
        );
    }

    displacement
    {
        type patch;
        faces
        (
            (D3 D2 C2 C3)
            (D0 C0 C1 D1)
        );
    }

    zeroTraction
    {
        type patch;
        faces
        (
            (D1 C1 C2 D2)
            (B2 B1 A1 A2)
        );
    }

    lowerCylContact
    {
        type patch;
        faces
        (
            (C0 D0 D3 C3)
        );
    }

    upperCylContact
    {
        type patch;
        faces
        (
            (B0 B3 A3 A0)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
