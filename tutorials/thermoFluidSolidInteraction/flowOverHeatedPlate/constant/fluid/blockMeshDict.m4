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

// GEOMETRY
// visina
//m4_define(H, 0.010)
// duljina
//udaljenost stijenke od izlaza
//m4_define(L, 0.05)
//udaljenost stijenke od ulaza
//m4_define(Lp, 0.025)
// pola sirine flapa
m4_define(l, 0.04)
// visina flapa
m4_define(h, 0.04)

// z koordinata
m4_define(zA, -0.001)
m4_define(zB, 0.001)

// MESH
// Abaqus mesh
//m4_define(BLOCKSIZE, 40 54 1)
//m4_define(BLOCKSIZE, 20 51 1) // from paper
//m4_define(BLOCKSIZE1,  calc(H*1000) calc((L-l)*2000) 1)

m4_define(BLOCKSIZE,  60 60 1)

m4_define(grading, 1 1 1)

//m4_define(r2,calc(r1+h))
//m4_define(angle, rad(45))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 1;

vertices
(
 //block1
    //Plane A:
    ( -l  0 zA) vlabel(A0)
    (  0  0 zA) vlabel(A1)
    (  0  h zA) vlabel(A2)
    ( -l  h zA) vlabel(A3)
    //Plane B:
    ( -l  0 zB) vlabel(B0)
    (  0  0 zB) vlabel(B1)
    (  0  h zB) vlabel(B2)
    ( -l  h zB) vlabel(B3)
 );

blocks
(
    hex (A0 A1 A2 A3 B0 B1 B2 B3) (BLOCKSIZE) simpleGrading (grading)
);

edges
(

);

boundary
(
    top
    {
        type patch;
        faces
        (
            (A2 A3 B3 B2)
        );
    }

    left
    {
        type patch;
        faces
        (
            (B0 B3 A3 A0)
        );
    }

    right
    {
        type patch;
        faces
        (
            (A1 A2 B2 B1)
        );
    }

    bottom
    {
        type patch;
        faces
        (
            (A0 A1 B1 B0)
        );
    }

    frontAndBack
    {
        type empty;
        faces
        (
            (A0 A3 A2 A1)
            (B0 B1 B2 B3)
        );
    }
);

//mergePatchPairs
//(
//);

// ************************************************************************* //
