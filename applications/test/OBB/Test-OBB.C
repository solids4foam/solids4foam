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

Application
    Test-OBB

Description
    Tests oriented bounding box fitting and overlap queries.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "OBB.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

bool contains(const OBB& box, const point& p, const scalar tolerance)
{
    const vector local(box.R() & (p - box.midpoint()));

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (mag(local[cmpt]) > box.ext()[cmpt] + tolerance)
        {
            return false;
        }
    }

    return true;
}


bool check(const bool condition, const char* message)
{
    if (!condition)
    {
        Info<< "FAILED: " << message << nl;
        return false;
    }

    return true;
}

} // End anonymous namespace


int main()
{
    bool passed = true;
    const scalar tolerance = 1e-10;
    const scalar angle = 0.63;
    const scalar c = Foam::cos(angle);
    const scalar s = Foam::sin(angle);
    const tensor rotatedAxes
    (
        c,  s,  0,
       -s,  c,  0,
        0,  0,  1
    );

    const point centre(1.2, -0.7, 2.1);
    const vector halfLength(2.0, 0.8, 0.3);
    pointField points(8);
    label pointI = 0;

    for (label i = -1; i <= 1; i += 2)
    {
        for (label j = -1; j <= 1; j += 2)
        {
            for (label k = -1; k <= 1; k += 2)
            {
                const vector local
                (
                    i*halfLength.x(),
                    j*halfLength.y(),
                    k*halfLength.z()
                );
                points[pointI++] = centre + (local & rotatedAxes);
            }
        }
    }

    const OBB fittedBox(points);

    forAll(points, i)
    {
        passed = check
        (
            contains(fittedBox, points[i], tolerance),
            "PCA-fitted box does not contain every input point"
        ) && passed;
    }

    passed = check
    (
        mag(fittedBox.midpoint() - centre) < tolerance,
        "PCA-fitted box centre is incorrect"
    ) && passed;

    const OBB boxA(point::zero, vector(1, 1, 1));
    const OBB touchingBox(point(2, 0, 0), vector(1, 1, 1));
    const OBB separatedBox(point(2.01, 0, 0), vector(1, 1, 1));
    const OBB rotatedBox(point(1.2, 0, 0), vector(1, 0.2, 0.2), rotatedAxes);

    passed = check(boxA.overlaps(touchingBox), "Touching boxes must overlap")
        && passed;
    passed = check
    (
        !boxA.overlaps(separatedBox),
        "Separated axis-aligned boxes reported an overlap"
    ) && passed;
    passed = check
    (
        boxA.overlaps(rotatedBox) && rotatedBox.overlaps(boxA),
        "Rotated overlap is not detected symmetrically"
    ) && passed;

    pointField facePoints(4);
    facePoints[0] = point(0, 0, 0);
    facePoints[1] = point(4*c, 4*s, 0);
    facePoints[2] = point(4*c - 2*s, 4*s + 2*c, 0);
    facePoints[3] = point(-2*s, 2*c, 0);

    // A non-unit normal is intentional: the constructor normalises it.
    const OBB faceBox(facePoints, vector(0, 0, 3), 0.4, 0.2);
    const scalar fittedFaceArea = 4*faceBox.ext().y()*faceBox.ext().z();

    passed = check
    (
        mag(faceBox.ext().x() - 0.3) < tolerance,
        "Face normal half-length is incorrect"
    ) && passed;
    passed = check
    (
        mag(fittedFaceArea - 8.0) < tolerance,
        "Minimum-area face fit is incorrect"
    ) && passed;
    passed = check
    (
        mag(faceBox.midpoint().z() - 0.1) < tolerance,
        "Asymmetric face extrusion centre is incorrect"
    ) && passed;

    if (!passed)
    {
        return 1;
    }

    Info<< "All OBB tests passed" << nl;
    return 0;
}


// ************************************************************************* //
