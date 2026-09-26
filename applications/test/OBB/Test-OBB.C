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
    Tests oriented bounding box fitting, growth, containment and overlap.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "OBB.H"
#include "IOstreams.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include <cmath>
#include <limits>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

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

    OBB fittedBox(points);
    fittedBox.grow(tolerance);

    forAll(points, i)
    {
        passed = check
        (
            fittedBox.contains(points[i]),
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

    OBB grownBox(centre, halfLength, rotatedAxes);
    grownBox.grow(0);
    passed = check
    (
        grownBox == OBB(centre, halfLength, rotatedAxes),
        "Zero growth changed the box"
    ) && passed;
    grownBox.grow(0.2);
    passed = check
    (
        mag(grownBox.ext() - vector(2.2, 1.0, 0.5)) < tolerance
     && grownBox.midpoint() == centre
     && grownBox.R() == OBB(centre, halfLength, rotatedAxes).R(),
        "Absolute growth changed the centre, axes or expected half-lengths"
    ) && passed;
    passed = check
    (
        grownBox.contains(centre + (vector(2.1, 0, 0) & rotatedAxes))
     && !grownBox.contains(centre + (vector(2.3, 0, 0) & rotatedAxes))
     && boxA.contains(point(1, 1, 1)),
        "Containment failed for rotated points or an axis-aligned boundary"
    ) && passed;

    // Absolute growth must also create thickness in degenerate boxes.
    OBB planarBox(facePoints);
    planarBox.grow(0.2);
    const point faceCentre(0.5*(facePoints[0] + facePoints[2]));
    passed = check
    (
        planarBox.contains(faceCentre + vector(0, 0, 0.19))
     && !planarBox.contains(faceCentre + vector(0, 0, 0.21)),
        "Growth of a planar box did not create the expected thickness"
    ) && passed;
    OBB pointBox(pointField(1, centre));
    pointBox.grow(0.2);
    passed = check
    (
        pointBox.contains(centre)
     && mag(pointBox.ext() - vector(0.2, 0.2, 0.2)) < tolerance,
        "Growth of a single-point box failed"
    ) && passed;

    // Uniform coordinate scaling must preserve candidate membership.
    const scalar scale = 1000;
    OBB scaledBox(scale*centre, scale*halfLength, rotatedAxes);
    scaledBox.grow(scale*0.2);
    forAll(points, i)
    {
        passed = check
        (
            grownBox.contains(points[i])
         == scaledBox.contains(scale*points[i]),
            "Coordinate scaling changed containment"
        ) && passed;
    }

#ifdef OPENFOAM_COM
    boundBox axisAligned(centre - halfLength, centre + halfLength);
    OBB sameGrowth(axisAligned);
    axisAligned.grow(0.2);
    sameGrowth.grow(0.2);
    passed = check
    (
        mag(sameGrowth.midpoint() - axisAligned.midpoint()) < tolerance
     && mag(sameGrowth.ext() - 0.5*axisAligned.span()) < tolerance,
        "OBB growth differs from OpenFOAM.com absolute growth"
    ) && passed;
#endif

    OBB emptyBox((pointField()));
    OBB defaultBox;
    defaultBox.grow(100);
    passed = check
    (
        defaultBox == emptyBox && defaultBox.empty()
     && !defaultBox.contains(point::zero)
     && !defaultBox.overlaps(boxA) && !boxA.overlaps(defaultBox),
        "A default box must remain empty after growth and overlap nothing"
    ) && passed;
    emptyBox.grow(100);
    passed = check
    (
        emptyBox.empty() && !emptyBox.contains(point::zero)
     && !emptyBox.overlaps(boxA) && !boxA.overlaps(emptyBox),
        "An empty box must remain empty after growth and overlap nothing"
    ) && passed;

    for (label small = 0; small < 2; ++small)
    {
        const scalar lengthScale = small ? 1e-20 : 1;
        pointField line(3);
        forAll(line, i)
        {
            line[i] = lengthScale*scalar(i)*vector(1, 2, 3);
        }
        OBB lineBox(line);
        lineBox.grow(lengthScale*1e-10);
        forAll(line, i)
        {
            passed = check
            (
                lineBox.contains(line[i]),
                "A collinear or very small point set was not enclosed"
            ) && passed;
        }
        for (direction cmpt = 0; cmpt < tensor::nComponents; ++cmpt)
        {
            passed = check
            (
                std::isfinite(lineBox.R()[cmpt]),
                "Degenerate point set produced non-finite axes"
            ) && passed;
        }
    }

    tensor invalidAxes(tensor::I);
    invalidAxes.xx() = std::numeric_limits<scalar>::quiet_NaN();
    const OBB repairedBox(centre, halfLength, invalidAxes);
    passed = check
    (
        repairedBox.R() == tensor::I && repairedBox.contains(centre),
        "Non-finite axes were not repaired"
    ) && passed;

    // Exercise both stream formats, including an empty box in the same stream.
    for (label binary = 0; binary < 2; ++binary)
    {
        const IOstream::streamFormat format =
            binary ? IOstream::BINARY : IOstream::ASCII;
        OStringStream os(format);
        os.precision(16);
        os << grownBox << token::SPACE << emptyBox;
        IStringStream is(os.str(), format);
        OBB restored(is);
        OBB restoredEmpty(is);
        passed = check
        (
            mag(restored.midpoint() - grownBox.midpoint()) < tolerance
         && mag(restored.ext() - grownBox.ext()) < tolerance
         && mag(restored.R() - grownBox.R()) < tolerance
         && restoredEmpty.empty(),
            "OBB stream round-trip failed"
        ) && passed;
    }

    if (!passed)
    {
        return 1;
    }

    Info<< "All OBB tests passed" << nl;
    return 0;
}


// ************************************************************************* //
