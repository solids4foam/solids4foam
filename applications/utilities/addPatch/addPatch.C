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
    createCrackPatch

Description
    Adds an empty patch to the end of the boundary

    It overwrites the previous mesh

Author
    Ripu Manchanda UT
    Philip Cardiff UCD/UT

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "directTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
label findPatchID(const polyMesh& mesh, const word& name)
{
    label patchI = mesh.boundaryMesh().findPatchID(name);

    if (patchI > -1)
    {
        FatalErrorIn("findPatchID(const polyMesh&, const word&)")
            << "Patch " << name << "already exists" << endl
            << exit(FatalError);
    }
    return patchI;
}

int main(int argc, char *argv[])
{
    argList::validArgs.append("patchName");
    argList::validArgs.append("patchType");
    argList::validOptions.insert("additionalPatches", "(patch2 .. patchN)");
    argList::validOptions.insert("additionalPatchTypes", "(type2 .. typeN)");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    const word patchName(args.additionalArgs()[0]);
    const word patchType(args.additionalArgs()[1]);

    wordList allPatches(1, patchName);
    wordList allPatchTypes(1, patchType);

    // Additional patches
    if (args.optionFound("additionalPatches"))
    {
        const wordList additionalPatchNames
        (
            args.optionLookup("additionalPatches")()
        );
        allPatches.append(additionalPatchNames);

        // Additional patch types (defaults to patch)
        if (args.optionFound("additionalPatchTypes"))
        {
            allPatchTypes.append
            (
                wordList(args.optionLookup("additionalPatchTypes")())
            );
        }
        else
        {
            allPatchTypes.append
            (
                wordList(additionalPatchNames.size(), "patch")
            );
        }
    }

    const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

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

    label patchIndex = oldPatches.size();
    forAll (allPatches, patchI)
    {
        findPatchID(mesh, allPatches[patchI]);

        Info<< "Adding " << allPatches[patchI]
            << " patch to the end of the boundary patch list with patch type "
            << allPatchTypes[patchI] << endl;

        // Add new patch
        newPatches.append
        (
            polyPatch::New
            (
                allPatchTypes[patchI],           // patch type
                allPatches[patchI],              // name
                0,                    // size
                startFaceI,           // start
                patchIndex,    // patch index
                oldPatches            // polyBoundaryMesh
            ).ptr()
        );

        patchIndex++;
    }

    //// Add new patch
    //newPatches.append
    //(
    //    polyPatch::New
    //    (
    //        patchType,           // patch type
    //        patchName,              // name
    //        0,                    // size
    //        startFaceI,           // start
    //        oldPatches.size(),    // patch index
    //        oldPatches            // polyBoundaryMesh
    //    ).ptr()
    //);

    // Remove the old patches and add the new patches
    newPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(newPatches);

    // Write the mesh

    Info<< "Writing the mesh" << endl;

    mesh.setInstance(mesh.pointsInstance());
    mesh.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
