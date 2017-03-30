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

        projectPatchToSphere leftPatch "(0 0 0)" 0.5

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "fvMesh.H"
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
    Foam::argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"

    const word patchName(args.additionalArgs()[0]);

    const vector origin(IStringStream(args.additionalArgs()[1])());

    const scalar radius(readScalar(IStringStream(args.additionalArgs()[2])()));

    if (radius < SMALL)
    {
        FatalError
            << "radius must be greater than zero"
            << abort(FatalError);
    }

    const bool overwrite = args.optionFound("overwrite");

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

    // Calculate point motion field

    pointMesh pMesh(mesh);

    pointVectorField pointMotion
    (
        IOobject
        (
            "pointMotion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Patch points
    const vectorField localPoints = mesh.boundaryMesh()[patchID].localPoints();

    // Patch point addressing
    const labelList meshPoints = mesh.boundaryMesh()[patchID].meshPoints();

    // New mesh points
    vectorField newPoints = mesh.points();

    // Calculate the new position of each point on the patch
    forAll(localPoints, pointI)
    {
        // Vector from the origin to the patch point
        const vector r = localPoints[pointI] - origin;
        scalar magR = mag(r);

        if (magR < SMALL)
        {
            FatalError
                << "The sphere origin cannot be on the original patch"
                << abort(FatalError);
        }

        // Scale point to the correct radius
        const vector newPoint = origin + radius*(r/magR);

        const label pointID = meshPoints[pointI];

        newPoints[pointID] = newPoint;
        pointMotion[pointID] = newPoint - localPoints[pointI];
    }

    // Write mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    else
    {
        runTime++;
    }

    mesh.movePoints(newPoints);
    mesh.moving(false);
    mesh.changing(false);
    mesh.setPhi().writeOpt() = IOobject::NO_WRITE;

    Info<< "Writing mesh and pointMotion field to time "
        << runTime.value() << endl;

    mesh.write();
    pointMotion.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
