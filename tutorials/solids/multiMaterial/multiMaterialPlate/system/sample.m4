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
    object      sample;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// M4
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'printf ($1)')])

// Thickness of the mid-layer relative to W: 0, 0.002, 0.005 or 0.01. Allrun
// sets it with "m4 -DhOverW=<value>"
m4_ifdef([hOverW], [], [m4_define(hOverW, 0.01)])

// The first argument is kept if there is a mid-layer, the second if not
m4_define(ifMidLayer, [m4_ifelse(hOverW, 0, [$2], [$1])])

// GEOMETRY, in m
m4_define(W, 0.03)
m4_define(h, calc(hOverW*W))

// The lines run along the interfaces from the free edge to the symmetry
// plane, 0.05 micrometres inside materials 1 and 2 rather than on the
// interface faces themselves: this is within the first cell next to the
// interface, which is 0.33 micrometres wide
m4_define(offset, 5e-8)
m4_define(xInterface13, calc(-1.5*W - h - offset))
m4_define(xInterface23, calc(-1.5*W + offset))
m4_define(yBottom, calc(-W))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type            sets;
libs            ("libsampling.so");

interpolationScheme cellPointFace;

setFormat       raw;

surfaceFormat   vtk;

fields
(
    sigma
);

sets
(
    ifMidLayer([interface13], [interface12])
    {
        type            face;
        axis            y;
        start           (xInterface13 yBottom 0.0005);
        end             (xInterface13 0 0.0005);
    }ifMidLayer([
    interface23
    {
        type            face;
        axis            y;
        start           (xInterface23 yBottom 0.0005);
        end             (xInterface23 0 0.0005);
    }])
);

surfaces
(
);

// ************************************************************************* //
