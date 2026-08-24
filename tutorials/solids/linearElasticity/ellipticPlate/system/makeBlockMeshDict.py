#!/usr/bin/env python3
#*--------------------------------*- C++ -*----------------------------------*#
# solids4foam: solid mechanics and fluid-solid interaction simulations        #
# Version:     v2.3                                                           #
# Web:         https://solids4foam.github.io                                  #
# Disclaimer:  This offering is not approved or endorsed by OpenCFD Limited,  #
#              producer and distributor of the OpenFOAM software via          #
#              www.openfoam.com, and owner of the OPENFOAM® and OpenCFD®      #
#              trade marks.                                                   #
#*---------------------------------------------------------------------------*#
'''
Description:
    Generate OpenFOAM spline blockEdge entries for the ellipticPlate
    blockMeshDict from the analytical inner and outer ellipse equations.

Author:
    Ivan Batistic, UCD.
'''

import argparse
import math


def fmt(value):
    if abs(value) < 5e-13:
        value = 0.0

    text = f"{value:.10f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def ellipse_points(a, b, z, n_points, reverse=False):
    angles = [0.5*math.pi*i/(n_points - 1) for i in range(n_points)]
    points = [(a*math.cos(t), b*math.sin(t), z) for t in angles]

    if reverse:
        points.reverse()

    return points


def print_spline_edge(name, a, b, z, n_points, reverse=False):
    print(f"    {name}")
    print("    (")

    for x, y, z in ellipse_points(a, b, z, n_points, reverse):
        print(f"        ({fmt(x)} {fmt(y)} {fmt(z)})")

    print("    )")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=
        "Generate spline edge entries for the ellipticPlate OpenFOAM "
        "blockMeshDict."
    )
    parser.add_argument(
        "-n",
        "--n-points",
        type=int,
        default=33,
        help="Number of sampled points per spline edge, including endpoints.",
    )
    args = parser.parse_args()

    if args.n_points < 3:
        raise SystemExit("At least 3 points are required for each spline edge")

    print_spline_edge("spline 0 1", 2.0, 1.0, 0.0, args.n_points, True)
    print_spline_edge("spline 4 5", 2.0, 1.0, 0.6, args.n_points, True)
    print_spline_edge("spline 2 3", 3.25, 2.75, 0.0, args.n_points)
    print_spline_edge("spline 6 7", 3.25, 2.75, 0.6, args.n_points)
