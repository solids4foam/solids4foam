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

Description
    Project points of the specified patch onto the surface of a sphere.

    The sphere is specified by origin and radius, e.g.

        projectPatchToSphere leftPatch (0 0 0) 0.5

    Take care as this does not move internal points so it may cause invalid
    cells.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "pointMesh.H"
#include "pointFields.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("patch name");
    argList::validArgs.append("origin");
    argList::validArgs.append("radius");

#   include "setRootCase.H"

#ifdef OPENFOAM_NOT_EXTEND
    const word patchName(args[1]);
    #ifdef OPENFOAM_COM
    const vector origin(args.get<vector>(2));
    const scalar radius(args.get<scalar>(3));
    #else
    const vector origin(args.argRead<vector>(2));
    const scalar radius(args.argRead<scalar>(3));
    #endif
#else
    const word patchName(args.additionalArgs()[0]);
    const vector origin(IStringStream(args.additionalArgs()[1])());
    const scalar radius(readScalar(IStringStream(args.additionalArgs()[2])()));
#endif

    if (radius < SMALL)
    {
        FatalError
            << "radius must be greater than zero"
            << abort(FatalError);
    }

    Info<< "Patch:" << patchName << nl
        << "origin:" << origin << nl
        << "Radius:" << radius << endl;

#   include "createTime.H"
#   include "createMesh.H"

    const word oldInstance = mesh.pointsInstance();


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

    // Create pointMesh so we can create a pointVectorField for visualisation
    pointMesh pMesh(mesh);

    // Create pointVectorField to visualise the mesh motion
    pointVectorField pointMotion
    (
        IOobject
        (
            "pointMotion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );
    vectorField& pointMotionI = pointMotion;


    Info<< "Reading points field from " << points.path() << nl << endl;

    // Calculate point motion field

    // Patch points
    const vectorField localPoints = mesh.boundaryMesh()[patchID].localPoints();

    // Patch point addressing
    const labelList meshPoints = mesh.boundaryMesh()[patchID].meshPoints();

    // Calculate the new position of each point on the patch

    scalar maxDisp = 0.0;

    forAll(localPoints, pointI)
    {
        // Global point index
        const label pointID = meshPoints[pointI];

        // Old point
        //const point oldPoint = localPoints[pointI];
        const point oldPoint = points[pointID];

        // Vector from the origin to the patch point
        const vector r = oldPoint - origin;
        scalar magR = mag(r);

        if (magR < SMALL)
        {
            FatalError
                << "The sphere origin cannot be on the original patch"
                << abort(FatalError);
        }

        // Scale point to the correct radius
        points[pointID] = origin + radius*(r/magR);

        // Store motion
        pointMotionI[pointID] = points[pointID] - oldPoint;

        // Store max motion
        maxDisp = max(maxDisp, mag(pointMotionI[pointID]));
    }

    Info<< "Maximum point displacement: " << maxDisp << nl << endl;

    // Increase the write precision
    IOstream::defaultPrecision(16);

    // Write the points
    Info<< "Writing the new points to " << points.path() << endl;
    points.write();

    Info<< "Writing the " << pointMotion.name() << " field to "
        << runTime.timeName() << endl;
    pointMotion.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
