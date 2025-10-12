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

// Model Description
// Axisymmetric cylinder mesh with compression die

// Setup m4 stuff
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])

// define geometry in mm
define(PI, 3.14159265359)
define(l, 1) // length of cylinder section
define(ri, 10) // cylinder inner radius
define(ro, 20) // cylinder outer radius
define(theta, 1) // angle of wedge in degrees

// calculated quantities
define(thetaRad2, calc(0.5*theta*PI/180))
define(sint, calc(sin(thetaRad2)))
define(cost, calc(cos(thetaRad2)))

// define mesh density
define(mr, 10) // number of cells in radial direction

// start of blockMeshDict

scale 0.001;

vertices
(
    //- dimension in mm
    (0 calc(ri*cost) calc(-ri*sint))
    (l calc(ri*cost) calc(-ri*sint))
    (l calc(ro*cost) calc(-ro*sint))
    (0 calc(ro*cost) calc(-ro*sint))

    (0 calc(ri*cost) calc(ri*sint))
    (l calc(ri*cost) calc(ri*sint))
    (l calc(ro*cost) calc(ro*sint))
    (0 calc(ro*cost) calc(ro*sint))
);

blocks
(
    hex (0 1 2 3 4 5 6 7) billet (1 mr 1) simpleGrading (1 1 1)
);

edges
(
);

patches
(
    symmetryPlane symmPlane
    (
        (0 4 7 3)
        (1 2 6 5)
    )

    patch inner
    (
        (0 1 5 4)
    )

    patch outer
    (
        (2 3 7 6)
    )

    wedge back
    (
        (3 2 1 0)
    )

    wedge front
    (
        (4 5 6 7)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
