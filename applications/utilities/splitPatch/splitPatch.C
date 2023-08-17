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
    splitPatch

Description
    Splits up a patch by putting faces in the given bounding box in a new patch

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "polyTopoChange.H"
#else
    #include "directTopoChange.H"
#endif
#include "boundBox.H"


// * * * * * * * * * * * * *  Main Program * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    runTime.functionObjects().off();

    const bool overwrite = args.optionFound("overwrite");

    const word oldInstance = mesh.pointsInstance();


    // Read dictionary
    const IOdictionary dict
    (
        IOobject
        (
            "splitPatchDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    // Read patch to split
    const word patchToSplitName(dict.lookup("patchToSplitName"));
    Info<< "Reading patchToSplit: " << patchToSplitName << nl << endl;

    const label patchToSplitID =
        mesh.boundaryMesh().findPatchID(patchToSplitName);

    if (patchToSplitID == -1)
    {
        FatalError
            << "patch " << patchToSplitName << " not found"
            << abort(FatalError);
    }

    // New patch name
    const word newPatchName(dict.lookup("newPatchName"));

    // Read the bounding boxes
    // Faces within the bounding boxes will be put in a new patch
    const List<boundBox> boundBoxes = dict.lookup("boundBoxes");

    Info<< "Adding a new patch: " << newPatchName << nl << endl;
    {
        // Add a new patch
        // Copy old patches and increase size by 1

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        const label nOldPatches = patches.size();
        const label newPatchID = nOldPatches;

        DynamicList<polyPatch*> allPatches(nOldPatches + 1);

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            allPatches.append
            (
                pp.clone
                (
                    patches,
                    allPatches.size(),
                    pp.size(),
                    pp.start()
                ).ptr()
            );
        }

        // Add new patch at the end
        allPatches.append
        (
            new polyPatch
            (
                newPatchName,     // patch name
                0,                // size
                mesh.nFaces(),    // start face index
                newPatchID,       // patch index
                patches           // boundary mesh
#ifdef OPENFOAM
                ,
                "patch"           // patch type
#endif
            )
        );

        Info<< "Removing patches." << endl;
        allPatches.shrink();
        mesh.removeBoundary();
        mesh.addPatches(allPatches);
    }

    // Find faces within the bounding boxes
    Info<< "Find faces inside bounding boxes" << endl;
    const polyPatch& ppatch = mesh.boundaryMesh()[patchToSplitID];
    const vectorField& patchCf = ppatch.faceCentres();
    SLList<label> facesToSplit;
    forAll(boundBoxes, boxI)
    {
        const boundBox& bb = boundBoxes[boxI];

        forAll(ppatch, faceI)
        {
            if (bb.contains(patchCf[faceI]))
            {
                facesToSplit.insert(faceI);
            }
        }
    }



    // Add faces to the new patch
    const labelList facesToSplitList(facesToSplit);
    Info<< "    faces found: " << facesToSplitList.size() << endl;

    const label newPatchID = mesh.boundaryMesh().size() - 1;

    // Create topo changer to perform mesh changes
#ifdef OPENFOAM_NOT_EXTEND
    polyTopoChange meshMod(mesh);
#else
    directTopoChange meshMod(mesh);
#endif

    forAll(facesToSplitList, faceI)
    {
        // Add face to the new patch

        const label curFaceID = ppatch.start() + facesToSplitList[faceI];

        meshMod.modifyFace
        (
            mesh.faces()[curFaceID],         // face
            curFaceID,                       // face ID,
            mesh.faceOwner()[curFaceID],     // owner,
            -1,                              // neighbour,
            false,                           // flip face
            newPatchID,                      // patchID,
            -1,                              // zone ID
            false                            // zone flip
        );
    }

    // Perform the topo change to move the faces to the new
    // patch
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);

    // Update mesh fields
    //mesh.updateMesh(map);


    // Write the mesh

    Info<< nl << "Writing the mesh" << endl;

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    mesh.write();


    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
