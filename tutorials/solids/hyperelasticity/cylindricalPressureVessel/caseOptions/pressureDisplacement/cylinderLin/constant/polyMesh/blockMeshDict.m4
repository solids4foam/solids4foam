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
// M4
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'printf ($1)')])
m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

// GEOMETRY
// Inner radius
m4_define(r1, 8)
// thickness
m4_define(h, 8)
// lenght
m4_define(l, 0.5)

// MESH
// Abaqus mesh
//m4_define(BLOCKSIZE, 40 54 1)
//m4_define(BLOCKSIZE, 20 51 1) // from paper
m4_define(BLOCKSIZE, 10 26 1)
m4_define(grading, 3 1 1)

m4_define(zA, -0.5)
m4_define(zB, l)
m4_define(r2,calc(r1+h))
m4_define(angle, rad(45))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    //Plane A:
    (r1 0 zA) vlabel(A0)
    (r2 0 zA) vlabel(A1)
    (0 r2 zA) vlabel(A2)
    (0 r1 zA) vlabel(A3)

    //Plane B:
    (r1 0 zB) vlabel(B0)
    (r2 0 zB) vlabel(B1)
    (0 r2 zB) vlabel(B2)
    (0 r1 zB) vlabel(B3)
);

blocks
(
    hex ( A0 A1 A2 A3 B0 B1 B2 B3 ) (BLOCKSIZE) simpleGrading (grading)
);

edges
(
    // Plane A
    arc  A0 A3  (calc(r1*cos(angle)) calc(r1*cos(angle)) zA)
    arc  A1 A2  (calc(r2*cos(angle)) calc(r2*cos(angle)) zA)
    // Plane B
    arc  B0 B3  (calc(r1*cos(angle)) calc(r1*cos(angle)) zB)
    arc  B1 B2  (calc(r2*cos(angle)) calc(r2*cos(angle)) zB)
);

boundary
(
    inside
    {
        type patch;
        faces
        (
            (A3 B3 B0 A0)
        );
    }

    outside
    {
        type patch;
        faces
        (
            (A1 B1 B2 A2)
        );
    }

    frontAndBack
    {
        type empty;
        faces
        (
            (A0 A1 A2 A3)
            (B1 B2 B3 B0)
        );
    }
/*
    symm1
    {
       type symmetryPlane;
        faces
        (
            (A0 A1 A2 A3)
        );
    }
    symm2
    {
        type symmetryPlane;
        faces
        (
            (B1 B2 B3 B0)
        );
    }
*/
    symm3
    {
        type symmetryPlane;
        faces
        (
            (A2 B2 B3 A3)
        );
    }

    symm4
    {
        type symmetryPlane;
        faces
        (
            (B1 A1 A0 B0)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
