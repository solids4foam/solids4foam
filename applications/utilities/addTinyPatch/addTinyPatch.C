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
    addTinyPatch

Description
    Given a vector coordinate, find the closest boundary face on the specified
    patch and add it to a new patch.

    The utilty requires three arguements:

        addTinyPatch currentPatchName newTinyPatchName "(1 2 3)"

    where "(1 2 3)" is the coordinate used to find the closest boundary face on
    the "currentPatchName" patch.

    The new mesh overwrites the previous mesh.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "polyTopoChange.H"
#else
    #include "directTopoChange.H"
#endif
#include "polyModifyFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::validArgs.append("currentPatchName");
    argList::validArgs.append("newTinyPatchName");
    argList::validArgs.append("approximateFaceCoordinate");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    // Store points instance (i.e. where the mesh was read from)
    const word pointsInstance = mesh.pointsInstance();

    // Read arguments
#ifdef OPENFOAMESIORFOUNDATION
    const word currentPatchName(args[0]);
    const word newTinyPatchName(args[1]);
    #ifdef OPENFOAMESI
    const vector approxFaceC(args.get<vector>(2));
    #else
    const vector approxFaceC(args.argRead<vector>(2));
    #endif
#else
    const word currentPatchName(args.additionalArgs()[0]);
    const word newTinyPatchName(args.additionalArgs()[1]);
    const vector approxFaceC(IStringStream(args.additionalArgs()[2])());
#endif

    // Method:
    // Step 1: Add new patch with no faces to the end of the boundary list
    // Step 2: Find the closest boundary face and move it to the new patch


    // Step 1: Add new patch with no faces to the end of the boundary list

    // Previous boundary patches
    const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

    // Create new boundary patches
    DynamicList<polyPatch*> newPatches(oldPatches.size() + 1);

    // Copy old patches

    label startFaceI = mesh.nInternalFaces();

    forAll(oldPatches, patchI)
    {
        const polyPatch& pp = oldPatches[patchI];

        newPatches.append
        (
            pp.clone
            (
                oldPatches,
                patchI,
                pp.size(),
                startFaceI
            ).ptr()
        );

        startFaceI += pp.size();
    }

    // Append new patch to the end

    const label newPatchID = oldPatches.size();

    newPatches.append
    (
        polyPatch::New
        (
            "patch",              // patch type
            newTinyPatchName,     // name
            0,                    // size
            startFaceI,           // start face index
            newPatchID,           // patch index
            oldPatches            // polyBoundaryMesh
        ).ptr()
    );


    // Remove the old patches and add the new patches
    newPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(newPatches);


    // Step 2: Find the closest boundary face and move it to the new patch

    // Create topo changer to perform mesh changes
    // This class records what changes are to be made to the mesh
#ifdef OPENFOAMESIORFOUNDATION
    polyTopoChange meshMod(mesh);
#else
    directTopoChange meshMod(mesh);
#endif

    // Find patch ID of specified patch

    const label curPatchID = mesh.boundaryMesh().findPatchID(currentPatchName);

    if (curPatchID == -1)
    {
        FatalError
            << currentPatchName << " patch not found!" << abort(FatalError);
    }

    // Find the face on the patch closest to approxFaceC

    const vectorField& patchCf = mesh.boundaryMesh()[curPatchID].faceCentres();

    label closestFaceID = -1;
    scalar closestDist = GREAT;

    forAll(patchCf, faceI)
    {
        scalar dist = mag(approxFaceC - patchCf[faceI]);

        if (dist < closestDist)
        {
            closestDist = dist;
            closestFaceID = faceI;
        }
    }

    if (closestFaceID == -1)
    {
        FatalError
            << "No closest face found! "
            << "Are you sure the specified patch has faces?"
            << abort(FatalError);
    }

    // Change local face index to global face index
    closestFaceID += mesh.boundaryMesh()[curPatchID].start();

    // Set face to be moved to the new patch
    meshMod.modifyFace
    (
        mesh.faces()[closestFaceID],        // face
        closestFaceID,                      // face ID,
        mesh.faceOwner()[closestFaceID],    // owner
        -1,                                 // neighbour,
        false,                              // flip face
        newPatchID,                         // patchID,
        -1,                                 // zone ID
        false                               // zone flip
    );

    // Perform the mesh change
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    // this->updateMesh(map);
    // if (map().hasMotionPoints())
    // {
    //     movePoints(map().preMotionPoints());
    // }
    meshMod.changeMesh(mesh, true);


    // Write the mesh
    Info<< "Overwriting the mesh" << endl;
    mesh.setInstance(pointsInstance);
    mesh.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
