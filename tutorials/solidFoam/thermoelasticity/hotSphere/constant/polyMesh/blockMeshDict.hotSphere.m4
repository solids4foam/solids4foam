/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM Extend Project: Open Source CFD        |
|  \\    /   O peration     | Version:  1.6-ext                               |
|   \\  /    A nd           | Web:      www.extend-project.de                 |
|    \\/     M anipulation  |                                                 |
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
// One eighth of a hollow sphere cylinder

// Setup m4 stuff
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])

// Geometry parameters in metres
define(PI, 3.14159265359)
define(ri, 0.19) // sphere inner radius
define(ro, 0.2) // sphere outer radius

// For the curves, we need to define
define(theta, 45)
define(thetaRad, calc(theta*PI/180))
define(thetaRad2, calc(0.5*theta*PI/180))
define(sint, calc(sin(thetaRad)))
define(cost, calc(cos(thetaRad)))
define(sint2, calc(sin(thetaRad2)))
define(cost2, calc(cos(thetaRad2)))

// define mesh density: ratio mc/mr = 15 for uniform cells
define(mr, 5) // number of cells in radial direction
//define(mc, 75) // half the number of cells in circumferential direction
define(mc, 25) // half the number of cells in circumferential direction

// start of blockMeshDict

convertToMeters 1;

vertices
(
    (0 0 ri)
    (calc(ri*cost) 0 calc(ri*cost))
    (ri 0 0)
    (calc(ri*cost) calc(ri*cost) 0)
    (0 ri 0)
    (0 calc(ri*cost) calc(ri*cost))
    (calc(ri/sqrt(3)) calc(ri/sqrt(3)) calc(ri/sqrt(3)))

    (0 0 ro)
    (calc(ro*cost) 0 calc(ro*cost))
    (ro 0 0)
    (calc(ro*cost) calc(ro*cost) 0)
    (0 ro 0)
    (0 calc(ro*cost) calc(ro*cost))
    (calc(ro/sqrt(3)) calc(ro/sqrt(3)) calc(ro/sqrt(3)))
);

blocks
(
    hex (0 1 6 5 7 8 13 12) billet (mc mc mr) simpleGrading (1 1 1)
    hex (1 2 3 6 8 9 10 13) billet (mc mc mr) simpleGrading (1 1 1)
    hex (5 6 3 4 12 13 10 11) billet (mc mc mr) simpleGrading (1 1 1)
);

edges
(
    arc 0 1 (calc(ri*cost2) 0 calc(ri*sint2))
    arc 7 8 (calc(ro*cost2) 0 calc(ro*sint2))

    arc 2 1 (calc(ri*sint2) 0 calc(ri*cost2))
    arc 9 8 (calc(ro*sint2) 0 calc(ro*cost2))

    arc 2 3 (calc(ri*cost2) calc(ri*sint2) 0)
    arc 10 9 (calc(ro*cost2) calc(ro*sint2) 0)

    arc 0 5 (0 calc(ri*sint2) calc(ri*cost2))
    arc 12 7 (0 calc(ro*sint2) calc(ro*cost2))

    arc 4 3 (calc(ri*sint2) calc(ri*cost2) 0)
    arc 10 11 (calc(ro*sint2) calc(ro*cost2) 0)

    arc 4 5 (0 calc(ri*cost2) calc(ri*sint2))
    arc 12 11 (0 calc(ro*cost2) calc(ro*sint2))

    arc 5 6 (calc(ri/3.0) calc(2.0*ri/3.0) calc(2.0*ri/3.0))
    arc 13 12 (calc(ro/3.0) calc(2.0*ro/3.0) calc(2.0*ro/3.0))

    arc 3 6 (calc(2.0*ri/3.0) calc(2.0*ri/3.0) calc(ri/3.0))
    arc 13 10 (calc(2.0*ro/3.0) calc(2.0*ro/3.0) calc(ro/3.0))

    arc 1 6 (calc(2.0*ri/3.0) calc(ri/3.0) calc(2.0*ri/3.0))
    arc 13 8 (calc(2.0*ro/3.0) calc(ro/3.0) calc(2.0*ro/3.0))
);

patches
(
    symmetryPlane symmx
    (
        (0 5 12 7)
        (5 4 11 12)
    )

    symmetryPlane symmy
    (
        (7 8 1 0)
        (8 9 2 1)
    )

    symmetryPlane symmz
    (
        (9 10 3 2)
        (10 11 4 3)
    )

    patch inside
    (
        (0 1 6 5)
        (1 2 3 6)
        (6 3 4 5)
    )

    patch outside
    (
        (7 12 13 8)
        (8 13 10 9)
        (12 11 10 13)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
