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

// Thickness of the mid-layer relative to W: 0.002, 0.005 or 0.01. Allrun
// sets it with "m4 -DhOverW=<value>"
m4_ifdef([hOverW], [], [m4_define(hOverW, 0.01)])

// GEOMETRY, in mm
m4_define(W, 30)
m4_define(h, calc(hOverW*W))
m4_define(xInterface23, calc(1.5*W + h))
m4_define(xRight, calc(3*W + h))

// MESH
// Materials 1 and 2 are graded towards the interfaces and the free edge with
// a ratio of 10000 between the largest and the smallest cell; the mid-layer
// is uniform in x
m4_define(nxMid, m4_ifelse(hOverW, 0.01, 500, 250))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 0.001;

vertices
(
    (0 0 0)                 // 0
    (45 0 0)                // 1
    (xInterface23 0 0)      // 2
    (xRight 0 0)            // 3
    (0 30 0)                // 4
    (45 30 0)               // 5
    (xInterface23 30 0)     // 6
    (xRight 30 0)           // 7
    (0 0 1)                 // 8
    (45 0 1)                // 9
    (xInterface23 0 1)      // 10
    (xRight 0 1)            // 11
    (0 30 1)                // 12
    (45 30 1)               // 13
    (xInterface23 30 1)     // 14
    (xRight 30 1)           // 15
);

blocks
(
    hex (0 1 5 4 8 9 13 12) material1 (120 240 1) simpleGrading (0.0001 10000 1)
    hex (2 3 7 6 10 11 15 14) material2 (120 240 1) simpleGrading (10000 10000 1)
    hex (1 2 6 5 9 10 14 13) material3 (nxMid 240 1) simpleGrading (1 10000 1)
);

edges
(
);

patches
(
    patch left
    (
        (0 8 12 4)
    )
    symmetry right
    (
        (3 7 15 11)
    )
    patch bottom
    (
        (0 1 9 8)
        (1 2 10 9)
        (2 3 11 10)
    )
    symmetry top
    (
        (4 12 13 5)
        (5 13 14 6)
        (6 14 15 7)
    )
    empty frontAndBack
    (
        (0 4 5 1)
        (1 5 6 2)
        (2 6 7 3)
        (8 9 13 12)
        (9 10 14 13)
        (10 11 15 14)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
