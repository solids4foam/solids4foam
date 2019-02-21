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


Description


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("thickness");

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    // Find empty patch
    label emptyPatchID = -1;
    
    forAll(mesh.boundaryMesh(), patchI)
    {
        if 
        (
            mesh.boundaryMesh()[patchI].type() 
         == emptyPolyPatch::typeName
        )
        {
            emptyPatchID = patchI;
        }
    }

    if(emptyPatchID == -1)
    {
        FatalErrorIn("...")
            << "Can not find empty patch"
                << abort(FatalError);
    }

    const labelList& meshPoints = 
        mesh.boundaryMesh()[emptyPatchID].meshPoints();

    // 2-D domain thickness
    scalar L(readScalar(IStringStream(args.additionalArgs()[0])()));
//     scalar L = 0.01;

    vectorField newPoints = mesh.points();

    boundBox box(newPoints);

    scalar zMin = box.min().z();
//     scalar zMax = box.max().z();

    forAll(meshPoints, pointI)
    {
        label curPoint = meshPoints[pointI];

        if (mag(newPoints[curPoint].z()-zMin)<SMALL)
        {
            newPoints[curPoint].z() = 0;
        }
        else
        {
            newPoints[curPoint].z() = L;
        }
    }

    twoDPointCorrector twoDCorrector(mesh);

    twoDCorrector.correctPoints(newPoints);
    
    mesh.movePoints(newPoints);

    runTime++;

    Info << "Writing scaled mesh ...";
    mesh.write();
    Info << "done" << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
