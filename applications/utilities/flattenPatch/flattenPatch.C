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

Application
    flattenPatch

Description
    Project points on given patch to a specified plane.

    The plane given by a point and normal.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("patch name");
    argList::validArgs.append("point");
    argList::validArgs.append("normal");

#   include "setRootCase.H"

    const word patchName(args.additionalArgs()[0]);

    const vector planePoint(IStringStream(args.additionalArgs()[1])());

    vector planeNormal(IStringStream(args.additionalArgs()[2])());

    if (mag(planeNormal) < SMALL)
    {
        FatalError
            << "The magnitude of the plane normal should be greater than zero!"
            << abort(FatalError);
    }

    // Normalise the normal
    planeNormal /= mag(planeNormal);

    Info<< "Patch:" << patchName << nl
        << "point:" << planePoint << nl
        << "normal:" << planeNormal << endl;

#   include "createTime.H"
#   include "createMesh.H"


    // Get patch ID

    const label patchID = mesh.boundaryMesh().findPatchID(patchName);

    if (patchID == -1)
    {
        FatalError
            << "Cannot find patch " << patchName
            << abort(FatalError);
    }

    // Read the points field
    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Info<< "Reading points field from " << points.path() << nl << endl;

    // Patch points
    const vectorField localPoints = mesh.boundaryMesh()[patchID].localPoints();

    // Patch point addressing
    const labelList meshPoints = mesh.boundaryMesh()[patchID].meshPoints();

    // Calculate the new position of each point on the patch

    scalar maxDisp = 0.0;

    forAll(localPoints, pointI)
    {
        const point& oldPoint = localPoints[pointI];

        // Vector from the planePoint to the current point
        const vector d = oldPoint - planePoint;

        // Component of d normal to the plane
        const vector nd = sqr(planeNormal) & d;

        // Point global ID
        const label pointID = meshPoints[pointI];

        // Project the point back to the plane
        points[pointID] = oldPoint - nd;

        // Record max displacement
        maxDisp = max(maxDisp, mag(nd));
    }

    Info<< "Maximum point displacement: " << maxDisp << nl << endl;

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    // Write the points
    Info<< "Writing the new points to " << points.path() << endl;
    points.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
