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
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

// Thickness of the mid-layer relative to W: 0, 0.002, 0.005 or 0.01. Allrun
// sets it with "m4 -DhOverW=<value>"
m4_ifdef([hOverW], [], [m4_define(hOverW, 0.01)])

// The first argument is kept if there is a mid-layer, the second if not
m4_define(ifMidLayer, [m4_ifelse(hOverW, 0, [$2], [$1])])

// GEOMETRY, in mm
// The origin is at the intersection of the two symmetry planes
m4_define(W, 30)
m4_define(h, calc(hOverW*W))
m4_define(xLeft, calc(-3*W - h))
m4_define(xInterface13, calc(-1.5*W - h))
m4_define(xInterface23, calc(-1.5*W))
m4_define(xRight, 0)
m4_define(yBottom, calc(-W))
m4_define(yTop, 0)

// MESH
// Materials 1 and 2 are graded towards the interfaces and the free edge with
// a ratio of 10000 between the largest and the smallest cell; the mid-layer
// is uniform in x
m4_define(nxSide, 120)
m4_define(nxMid, m4_ifelse(hOverW, 0.01, 500, 250))
m4_define(ny, 240)
m4_define(grading, 10000)
m4_define(invGrading, 0.0001)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 0.001;

vertices
(
    // Left boundary
    (xLeft yBottom 0) vlabel(L0)
    (xLeft yTop 0) vlabel(L1)
    (xLeft yBottom 1) vlabel(L2)
    (xLeft yTop 1) vlabel(L3)

    // Interface 23, which is interface 12 if there is no mid-layer
    (xInterface23 yBottom 0) vlabel(B0)
    (xInterface23 yTop 0) vlabel(B1)
    (xInterface23 yBottom 1) vlabel(B2)
    (xInterface23 yTop 1) vlabel(B3)

    // Right boundary
    (xRight yBottom 0) vlabel(R0)
    (xRight yTop 0) vlabel(R1)
    (xRight yBottom 1) vlabel(R2)
    (xRight yTop 1) vlabel(R3)
ifMidLayer([
    // Interface 13
    (xInterface13 yBottom 0) vlabel(A0)
    (xInterface13 yTop 0) vlabel(A1)
    (xInterface13 yBottom 1) vlabel(A2)
    (xInterface13 yTop 1) vlabel(A3)], [m4_define(A0, B0)m4_define(A1, B1)m4_define(A2, B2)m4_define(A3, B3)])
);

blocks
(
    hex (L0 A0 A1 L1 L2 A2 A3 L3) material1 (nxSide ny 1) simpleGrading (invGrading grading 1)
    hex (B0 R0 R1 B1 B2 R2 R3 B3) material2 (nxSide ny 1) simpleGrading (grading grading 1)ifMidLayer([
    hex (A0 B0 B1 A1 A2 B2 B3 A3) material3 (nxMid ny 1) simpleGrading (1 grading 1)])
);

edges
(
);

patches
(
    patch left
    (
        (L0 L2 L3 L1)
    )
    symmetry right
    (
        (R0 R1 R3 R2)
    )
    patch bottom
    (
        (L0 A0 A2 L2)
        (B0 R0 R2 B2)ifMidLayer([
        (A0 B0 B2 A2)])
    )
    symmetry top
    (
        (L1 L3 A3 A1)
        (B1 B3 R3 R1)ifMidLayer([
        (A1 A3 B3 B1)])
    )
    empty frontAndBack
    (
        (L0 L1 A1 A0)
        (L2 A2 A3 L3)
        (B0 B1 R1 R0)
        (B2 R2 R3 B3)ifMidLayer([
        (A0 A1 B1 B0)
        (A2 B2 B3 A3)])
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
