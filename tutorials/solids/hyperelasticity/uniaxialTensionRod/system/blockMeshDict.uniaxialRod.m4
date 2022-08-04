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
define(l, calc(0.5*53.334)) // half length of cylinder
define(r, 6.413) // cylinder radius
define(theta, 1) // angle of wedge in degrees

// calculated quantities
define(thetaRad2, calc(0.5*theta*PI/180))
define(sint, calc(sin(thetaRad2)))
define(cost, calc(cos(thetaRad2)))

// define mesh density
define(ml1, 40) // number of cells in length direction (refined region)
define(ml2, 20) // number of cells in length direction
define(mr, 10) // number of cells in radial direction
define(lx, calc(0.35*l)) // x-coordinate where refined region ends

// start of blockMeshDict

convertToMeters 0.001;

vertices
(
    //- dimension in mm
    (0 0 0)
    (lx 0 0)
    (l 0 0)
    (l calc(r*cost) calc(-r*sint))
    (lx calc(r*cost) calc(-r*sint))
    (0 calc(r*cost) calc(-r*sint))
    (l calc(r*cost) calc(r*sint))
    (lx calc(r*cost) calc(r*sint))
    (0 calc(r*cost) calc(r*sint))
);

blocks
(
    hex (0 1 4 5 0 1 7 8) billet (ml1 mr 1) simpleGrading (1 1 1)
    hex (1 2 3 4 1 2 6 7) billet (ml2 mr 1) simpleGrading (1 1 1)
);

edges
(
);

patches
(   symmetry symmPlane
    (
        (0 8 5 0)
    )

    patch loading
    (
        (2 6 3 2)
    )

    patch tracFree
    (
        (3 4 7 6)
        (4 5 8 7)
    )

    symmetry back
    (
        (5 4 1 0)
        (4 3 2 1)
    )

    symmetry front
    (
        (0 1 7 8)
        (1 2 6 7)
    )

    empty axis
    (
        (0 1 1 0) 
        (1 2 2 1) 
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
