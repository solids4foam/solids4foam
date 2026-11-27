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

// Thickness of the mid-layer relative to W: 0.002, 0.005 or 0.01. Allrun
// sets it with "m4 -DhOverW=<value>"
m4_ifdef([hOverW], [], [m4_define(hOverW, 0.01)])

// The lines run along the interfaces from the free edge to the symmetry
// plane, 0.05 micrometres inside materials 1 and 2 rather than on the
// interface faces themselves: this is within the first cell next to the
// interface, which is 0.33 micrometres wide
m4_define(W, 0.03)
m4_define(h, calc(hOverW*W))
m4_define(xInterface23, calc(1.5*W + h + 5e-8))

// The mid-line runs across the plate at half its height. The part of it
// around the mid-layer is sampled again on its own, as the sample utility of
// foam-extend skips cells that are smaller than 1e-5 of the length of a line
m4_define(xRight, calc(3*W + h))
m4_define(xZoomStart, calc(1.5*W - h))
m4_define(xZoomEnd, calc(1.5*W + 2*h))

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
    interface13
    {
        type            face;
        axis            y;
        start           (0.04499995 0 0.0005);
        end             (0.04499995 0.03 0.0005);
    }
    interface23
    {
        type            face;
        axis            y;
        start           (xInterface23 0 0.0005);
        end             (xInterface23 0.03 0.0005);
    }
    midLine
    {
        type            face;
        axis            x;
        start           (0 0.015 0.0005);
        end             (xRight 0.015 0.0005);
    }
    midLineZoom
    {
        type            face;
        axis            x;
        start           (xZoomStart 0.015 0.0005);
        end             (xZoomEnd 0.015 0.0005);
    }
);

surfaces
(
);

// ************************************************************************* //
