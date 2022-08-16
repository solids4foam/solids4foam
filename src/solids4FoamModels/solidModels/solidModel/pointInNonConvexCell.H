/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tetCellList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


label orient3d(const point& a, const point& b, const point& c, const point& d)
{
    // Computes the sign of the signed volume
    // of the tetrahedron (a,b,c,d)
    // return Sign( dot(cross(b-a,c-a),d-a) )
    return Foam::sign
    (
        ((b - a)^(c - a)) & (d -a)
    );
}


bool checkIntersectionOfLineAndTriangle(const Pair<vector>& S, const triPointRef& t)
{
    // Intersect([q1,q2],(t1,t2,t3))
    //   s1 = orient3d(q1,t1,t2,t3)
    //   s2 = orient3d(q2,t1,t2,t3)
    //   // Test whether the two extermities of the segment
    //   // are on the same side of the supporting plane of
    //   // the triangle
    //   If(s1 == s2)
    //      return false
    //   End if
    //   // Now we know that the segment 'straddles' the supporing
    //   // plane. We need to test whether the three tetrahedra formed
    //   // by the segment and the three edges of the triangle have
    //   // the same orientation
    //   s3 = orient3d(q1,q2,t1,t2)
    //   s4 = orient3d(q1,q2,t2,t3)
    //   s5 = orient3d(q1,q2,t3,t1)
    //   return (s3 == s4 AND s4 == s5)
    // End

    const point& q1 = S.first();
    const point& q2 = S.second();

    const point& t1 = t.a();
    const point& t2 = t.b();
    const point& t3 = t.c();

    const label s1 = orient3d(q1, t1, t2, t3);
    const label s2 = orient3d(q2, t1, t2, t3);

    if (s1 == s2)
    {
        return false;
    }

    const label s3 = orient3d(q1, q2, t1, t2);
    const label s4 = orient3d(q1, q2, t2, t3);
    const label s5 = orient3d(q1, q2, t3, t1);

    return ((s3 == s4) && (s4 == s5));
}


bool pointInNonConvexCell(const polyMesh& mesh, const point& p, const label cellI)
{
    // Algorithm taken from https://stackoverflow.com/questions/44513525/
    // testing-whether-a-3d-point-is-inside-a-3d-polyhedron
    // Inside(Polyhedron P, point q)
    //     Segment S = [q, q+(0,0,1e30)]
    //     count = 0
    //     For each triangle T of P
    //         If Intersect(S,T)
    //             count = count + 1
    //         End if
    //     End for
    //     return odd(count)
    // End
    // Where a Segment is a pair of points

    const Pair<vector> S(p, p + vector(0, 0, 1e30));
    label count = 0;

    const cell& curCell = mesh.cells()[cellI];
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    // Loop through all polygonal faces
    forAll(curCell, fI)
    {
        const label faceID = curCell[fI];
        const face& curFace = faces[faceID];

        // We will decompose the face into triangles using the face centre point
        // and then loop through these triangular faces
        const point ctr = Foam::average(curFace.points(points));
        point nextPoint = vector::zero;
        forAll(curFace, pI)
        {
            nextPoint = points[curFace.nextLabel(pI)];

            // Construct triangle
            triPointRef t
            (
                points[curFace[pI]],
                nextPoint,
                ctr
            );

            //vector n = t.normal();
            //scalar A = mag(n);
            //n /= A;

            // Check for intersection between t and S
            if (checkIntersectionOfLineAndTriangle(S, t))
            {
                count++;
            }
        }
    }

    return bool(count % 2);
}


bool sameSide
(
    const point& v1,
    const point& v2,
    const point& v3,
    const point& v4,
    const point& p
)
{
    const vector normal = (v2 - v1)^(v3 - v1);
    const scalar dotV4 = (normal & (v4 - v1));
    const scalar dotP = (normal & (p - v1));
    return (Foam::sign(dotV4) == Foam::sign(dotP));
}


bool pointInTetrahedron
(
    const point& v1,
    const point& v2,
    const point& v3,
    const point& v4,
    const point& p
)
{
    return
    (
        sameSide(v1, v2, v3, v4, p)
     && sameSide(v2, v3, v4, v1, p)
     && sameSide(v3, v4, v1, v2, p)
     && sameSide(v4, v1, v2, v3, p)
    );
}


bool pointInCellTetDecomp
(
    const polyMesh& mesh,
    const tetCellList& tets,
    const pointField& tetPoints,
    const point& p,
    const label cellI
)
{
    forAll(tets, tetI)
    {
        const tetPointRef curTet = tets[tetI].tet(tetPoints);

        if
        (
            pointInTetrahedron
            (
                curTet.a(), curTet.b(), curTet.c(), curTet.d(), p
            )
        )
        {
            return true;
        }
    }

    return false;
}

} // End namespace Foam

// ************************************************************************* //
