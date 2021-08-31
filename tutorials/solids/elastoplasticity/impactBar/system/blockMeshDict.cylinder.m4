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
// Axisymmetric cylinder mesh with compression die

// Setup m4 stuff
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])

// define geometry in mm
define(PI, 3.14159265359)
define(l, 32.4) // full cylinder height/length
define(r, 3.2) // cylinder radius
define(theta, 1) // angle of wedge in degrees
define(dt, 1) // wall thickness
define(dw, 10) // wall radius

// calculated quantities
define(thetaRad2, calc(0.5*theta*PI/180))
define(sint, calc(sin(thetaRad2)))
define(cost, calc(cos(thetaRad2)))

// define mesh density
define(ml, 60) // number of cells in axial direction
define(mr, 6) // number of cells in radial direction
define(mr2, 1) // rigid wall

// start of blockMeshDict

convertToMeters 0.001;

vertices
(
    //- dimension in mm
    (0 0 0)
    (l 0 0)
    (l calc(r*cost) calc(-r*sint))
    (0 calc(r*cost) calc(-r*sint))
    (l calc(r*cost) calc(r*sint))
    (0 calc(r*cost) calc(r*sint))

    (l 0 0)
    (calc(l+dt) 0 0)
    (calc(l+dt) calc(dw*cost) calc(-dw*sint))
    (l calc(dw*cost) calc(-dw*sint))
    (calc(l+dt) calc(dw*cost) calc(dw*sint))
    (l calc(dw*cost) calc(dw*sint))
);

blocks
(
    hex (0 1 2 3 0 1 4 5) billet (ml mr 1) simpleGrading (1 1 1)
    hex (6 7 8 9 6 7 10 11) die (1 mr2 1) simpleGrading (1 1 1)
);

edges
(
);

patches
(
    patch left
    (
        (0 5 3 0)
    )

    patch dieContact
    (
        (6 11 9 6)
        (8 9 11 10)
    )

    patch billetContact
    (
        (1 4 2 1)
        (2 3 5 4)
    )

    patch loading
    (
        (7 10 8 7)
    )

    wedge back
    (
        (3 2 1 0)
        (9 8 7 6)
    )

    wedge front
    (
        (0 1 4 5)
        (6 7 10 11)
    )

    empty axis
    (
        (0 1 1 0)
        (6 7 7 6)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
