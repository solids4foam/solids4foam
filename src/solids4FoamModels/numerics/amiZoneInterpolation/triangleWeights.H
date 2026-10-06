/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Function
    Foam::triangleWeights

Description
    Barycentric weights of a point in a triangle, with a degeneracy test
    relative to the size of the triangle.

    triangle::pointToBarycentric treats a triangle as degenerate, and returns
    the weights (1/3, 1/3, 1/3), when d00*d11 - d01*d01 < SMALL. That quantity
    is four times the squared area, in units of length^4, so with SMALL = 1e-15
    every triangle with an area below about 1.6e-8 m^2 is taken as degenerate:
    on interface faces finer than about 0.3 mm the point interpolation becomes
    a three-point average of the vertices and the face centre. Here the test
    is relative: the triangle is degenerate when the squared sine of the angle
    between its edges from point a is below SMALL.

SourceFiles
    triangleWeights.H

\*---------------------------------------------------------------------------*/

#ifndef triangleWeights_H
#define triangleWeights_H

#include "triPointRef.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Barycentric weights (of a, b and c) of the point p in the triangle t
inline FixedList<scalar, 3> triangleWeights
(
    const triPointRef& t,
    const point& p
)
{
    const vector v0(t.b() - t.a());
    const vector v1(t.c() - t.a());
    const vector v2(p - t.a());

    const scalar d00 = v0 & v0;
    const scalar d01 = v0 & v1;
    const scalar d11 = v1 & v1;
    const scalar d20 = v2 & v0;
    const scalar d21 = v2 & v1;

    const scalar denom = d00*d11 - d01*d01;

    FixedList<scalar, 3> w;

    if (mag(denom) <= SMALL*d00*d11)
    {
        // Degenerate triangle
        w[0] = 1.0/3.0;
        w[1] = 1.0/3.0;
        w[2] = 1.0/3.0;
        return w;
    }

    w[1] = (d11*d20 - d01*d21)/denom;
    w[2] = (d00*d21 - d01*d20)/denom;
    w[0] = 1.0 - w[1] - w[2];

    return w;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
