/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    abaqusMeshToFoam

Description
    Converts an Abaqus mesh from an input file (*.inp) into the FOAM format.

    The approach taken is:
        - Count the nodes
        - Read the nodes and construct the pointField
        - Count the elements
        - Read the elements and construct the faceList and cellList
            - Each cell adds its own faces so faces will be duplicated
        - Add all faces to a default boundary patch
        - Remove duplicate faces using the procedures from mergeOrSplitBaffles

    Limitations:
        - Only the following element types are supported:
            - C3D8
            - C3D8R
        - Only the first PART is used and the rest are ignored
        - Node sets, element sets and surfaces are not converted

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "IFstream.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "Time.H"
    #include "polyTopoChange.H"
#else
    #include "foamTime.H"
    #include "directTopoChange.H"
#endif
#include "localPointRegion.H"
#include "processorPolyPatch.H"
#include "faceSet.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"

using namespace Foam;

// * * * * * * * * * * *  Helper Functions * * * * * * * * * * * * * * * * * //

// Return true if the header was succesfully read
bool readHeader(IFstream& inputStream)
{
    while (!inputStream.eof())
    {
        if (!inputStream.good())
        {
            break;
        }

        token tok(inputStream);

        // Look for the word "Part,"
        if (tok.isWord())
        {
            if (tok.wordToken() == "Part,")
            {
                // Read the next tok
                while (inputStream.read(tok))
                {
                    // Read the part name
                    const word partName(tok.wordToken());

                    if (tok.pToken() == token::MULTIPLY)
                    {
                        // Read the punctuation *
                        break;
                    }
                }

                //inputStream.read(tok);

                // Read the word "Node"
                inputStream.read(tok);
                if (tok.wordToken() != "Node")
                {
                    FatalErrorIn("readHeader(...)")
                        << ": error reading the node list!"
                        << exit(FatalError);
                }

                return true;
            }
        }
    }

    return false;
}


// Count the nodes and set the threeD flag
// Returns the number of nodes
label countNodes(IFstream& inputStream, bool& threeD)
{
    label nNodes = 0;
    bool checkThreeD = true;
    while (!inputStream.eof())
    {
        if (!inputStream.good())
        {
            break;
        }

        // Node list:
        // Starts with p *, then the word "Part," then another word (name of
        // part), then p "*", and then word "Node"
        // The each line of the node list will be:
        // label scalar scalar (for 2-D)
        // label scalar scalar scalar (for 3-D)
        // where label is the node index
        // Element list:
        // Then the element list starts with p *, word = "Element," and
        // word = "type=C3D8R" or "type=C3D8"
        // The p * indicates node sets etc.

        // Node format:
        // label, x, y        for 2-D
        // label, x, y, z     for 2-D

        // The node list will then end when the element list starts,
        // where the element list starts with a *

        bool error = false;

        // First, we will count the number of nodes
        token tok(inputStream);
        // Info<< tok.labelToken() << " " << __LINE__ << endl;
        if (!tok.isLabel())
        {
            error = true;
        }

        inputStream.read(tok);
        if (!tok.isPunctuation())
        {
            error = true;
        }

        inputStream.read(tok);
        // Info<< tok.doubleScalarToken() << " " << __LINE__ << endl;
        if (!tok.isDoubleScalar())
        {
            error = true;
        }

        inputStream.read(tok);
        if (!tok.isPunctuation())
        {
            error = true;
        }

        inputStream.read(tok);
        // Info<< tok.doubleScalarToken() << " " << __LINE__ << endl;
        if (!tok.isDoubleScalar())
        {
            error = true;
        }

        // Check if 3-D
        if (threeD)
        {
            inputStream.read(tok);
            if (!tok.isPunctuation())
            {
                error = true;
            }

            inputStream.read(tok);
            if (!tok.isDoubleScalar())
            {
                error = true;
            }
        }
        else if (checkThreeD)
        {
            checkThreeD = false;

            inputStream.read(tok);
            if (tok.isPunctuation())
            {
                Info<< "Model is 3-D" << endl;
                threeD = true;

                // Read z coord
                inputStream.read(tok);
            }
            else
            {
                Info<< "Model is 2-D" << endl;

                // Put the token back
                inputStream.putBack(tok);
            }
        }

        if (error)
        {
            FatalErrorIn("countNodes(...)")
                << ": error reading the node list!"
                << exit(FatalError);
        }

        nNodes++;

        // Check if we have reached the end of the node list
        inputStream.read(tok);
        if (tok.isPunctuation())
        {
            if (tok.pToken() == token::MULTIPLY)
            {
                Info<< "End of node list" << endl;
                break;
            }
            else
            {
                FatalErrorIn("countNodes(...)")
                    << ": found " << tok.pToken() << " in the node list!"
                    << exit(FatalError);
            }
        }
        else
        {
            // Put token back
            inputStream.putBack(tok);
        }
    }

    return nNodes;
}

// Read the nodes into the points field
// It is assumed the the header has already been read from the input stream
void readNodes
(
    IFstream& inputStream, const label nNodes, pointField& points, const bool threeD
)
{
    forAll(points, pointI)
    {
        // Node format:
        // label, x, y        for 2-D
        // label, x, y, z     for 2-D

        // Read node index
        token tok(inputStream);

        // Comma
        inputStream.read(tok);

        // X coord
        inputStream.read(tok);
        points[pointI][vector::X] = tok.doubleScalarToken();

        // Comma
        inputStream.read(tok);

        // Y coord
        inputStream.read(tok);
        points[pointI][vector::Y] = tok.doubleScalarToken();

        // Check if 3-D
        if (threeD)
        {
            // Comma
            inputStream.read(tok);

            // Z coord
            inputStream.read(tok);
            points[pointI][vector::Z] = tok.doubleScalarToken();
        }
    }
}


// Count and return the elements
label countElements(IFstream& inputStream, bool& threeD)
{
    // Read element list header
    token tok;

    // Read *
    inputStream.read(tok);

    // Read "Element,"
    inputStream.read(tok);

    // Read element type
    inputStream.read(tok);
    word elemType = tok.wordToken();
    // Remove "type="
    elemType = elemType.substr(5);
    Info<< "Element type = " << elemType << endl;
    if (elemType != "C3D8R" && elemType != "C3D8")
    {
        FatalError
            << "Unknown element type: " << elemType << endl;
    }

    label nElements = 0;
    while (!inputStream.eof())
    {
        if (!inputStream.good())
        {
            break;
        }

        // Element list:
        // Starts with p *, then the word "Element," then another word like
        // type=<ELEMENT_TYPE>, where ELEMENT_TYPE is the element type.
        // For now, we only support 3-D continuum elements: C3D8R and C3D8
        // Then the elements are listed
        // Element format:
        // element_index, i1, i2, i3, i4, i5, i6, i7, i8        for 2-D
        // where i* refer to the node indices in the element

        // The element list will then end with a * which is the start of a node
        // set or element set or something else

        bool error = false;


        // Read element index and 8 node indices
        for (label i = 0; i < 9; i++)
        {
            inputStream.read(tok);
            //Info<< "label " << tok.labelToken() << " " << endl;
            if (!tok.isLabel())
            {
                error = true;
            }

            if (i == 8)
            {
                break;
            }
            else
            {
                inputStream.read(tok);
                // Info<< "P " << tok.pToken() << " " << endl;
                if (!tok.isPunctuation())
                {
                    error = true;
                }
            }
        }

        if (error)
        {
            FatalErrorIn("countElements(...)")
                << ": error reading the element list!"
                << exit(FatalError);
        }

        nElements++;

        // Check if we have reached the end of the element list
        inputStream.read(tok);
        if (tok.isPunctuation())
        {
            if (tok.pToken() == token::MULTIPLY)
            {
                Info<< "End of element list" << endl;
                break;
            }
            else
            {
                FatalErrorIn("countElements(...)")
                    << ": found " << tok.pToken() << " in the element list!"
                    << exit(FatalError);
            }
        }
        else
        {
            // Put token back
            inputStream.putBack(tok);
        }
    }

    return nElements;
}


// Read the elements and populate the faces and cells lists
void readElements
(
    IFstream& inputStream, const label nElements, faceList& faces, cellList& cells
)
{
    // Read element list header
    token tok;

    // Read "Element,"
    inputStream.read(tok);

    // Read element type
    inputStream.read(tok);
    word elemType = tok.wordToken();
    // Remove "type="
    elemType = elemType.substr(5);
    Info<< "Element type = " << elemType << endl;
    if (elemType != "C3D8R" && elemType != "C3D8")
    {
        FatalError
            << "Unknown element type: " << elemType << endl;
    }

    label faceI = 0;
    forAll(cells, cellI)
    {
        // Element format:
        // element_index, i1, i2, i3, i4, i5, i6, i7, i8        for 2-D
        // where i* refer to the node indices in the element
        // Node indices start from 1 not from 0!

        // The element list will then end with a * which is the start of a node
        // set or element set or something else

        // Read element index
        inputStream.read(tok);

        // Read comma
        inputStream.read(tok);

        // Read the nodes in the element
        labelList p(label(8), label(-1));
        forAll(p, pI)
        {
            // Read node index
            // Subtract 1 as Abaqus counts from 1 rather than 0
            inputStream.read(tok);
            p[pI] = tok.labelToken() - 1;

            if (pI != (p.size() - 1))
            {
                // Read comma
                inputStream.read(tok);
            }
        }

        // Create faces
        faces[faceI].resize(4);
        faces[faceI][0] = p[3];
        faces[faceI][1] = p[2];
        faces[faceI][2] = p[1];
        faces[faceI][3] = p[0];

        faces[faceI + 1].resize(4);
        faces[faceI + 1][0] = p[4];
        faces[faceI + 1][1] = p[5];
        faces[faceI + 1][2] = p[6];
        faces[faceI + 1][3] = p[7];

        faces[faceI + 2].resize(4);
        faces[faceI + 2][0] = p[1];
        faces[faceI + 2][1] = p[2];
        faces[faceI + 2][2] = p[6];
        faces[faceI + 2][3] = p[5];

        faces[faceI + 3].resize(4);
        faces[faceI + 3][0] = p[2];
        faces[faceI + 3][1] = p[3];
        faces[faceI + 3][2] = p[7];
        faces[faceI + 3][3] = p[6];

        faces[faceI + 4].resize(4);
        faces[faceI + 4][0] = p[3];
        faces[faceI + 4][1] = p[0];
        faces[faceI + 4][2] = p[4];
        faces[faceI + 4][3] = p[7];

        faces[faceI + 5].resize(4);
        faces[faceI + 5][0] = p[0];
        faces[faceI + 5][1] = p[1];
        faces[faceI + 5][2] = p[5];
        faces[faceI + 5][3] = p[4];

        // Create cell
        labelList c(6);
        forAll(c, fI)
        {
            c[fI] = faceI++;
        }
        cells[cellI] = cell(c);

        // Check if we have reached the end of the element list
        inputStream.read(tok);
        if (tok.isPunctuation())
        {
            if (tok.pToken() == token::MULTIPLY)
            {
                Info<< "End of element list" << endl;
                break;
            }
            else
            {
                FatalErrorIn("countElements(...)")
                    << ": found " << tok.pToken() << " in the element list!"
                    << exit(FatalError);
            }
        }
        else
        {
            // Put token back
            inputStream.putBack(tok);
        }
    }
}


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces)
{
    // Get all duplicate face labels (in boundaryFaces indices!).
    labelList duplicates = localPointRegion::findDuplicateFaces
    (
        mesh,
        boundaryFaces
    );


    // Check that none are on processor patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(duplicates, bFaceI)
    {
        if (duplicates[bFaceI] != -1)
        {
            label faceI = mesh.nInternalFaces() + bFaceI;
            label patchI = patches.whichPatch(faceI);

            if (isA<processorPolyPatch>(patches[patchI]))
            {
                FatalErrorIn("findBaffles(const polyMesh&, const labelList&)")
                    << "Duplicate face " << faceI
                    << " is on a processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << faceI
                    << " is on patch:" << patches[patchI].name()
                    << abort(FatalError);
            }
        }
    }


    // Write to faceSet for ease of postprocessing.
    {
        faceSet duplicateSet
        (
            mesh,
            "duplicateFaces",
            (mesh.nFaces() - mesh.nInternalFaces())/256
        );

        forAll(duplicates, bFaceI)
        {
            label otherFaceI = duplicates[bFaceI];

            if (otherFaceI != -1 && otherFaceI > bFaceI)
            {
                duplicateSet.insert(mesh.nInternalFaces() + bFaceI);
                duplicateSet.insert(mesh.nInternalFaces() + otherFaceI);
            }
        }

        Pout<< "Writing " << duplicateSet.size()
            << " duplicate faces to faceSet " << duplicateSet.objectPath()
            << nl << endl;
        duplicateSet.write();
    }

    return duplicates;
}


void insertDuplicateMerge
(
    const polyMesh& mesh,
    const labelList& duplicates,
#ifdef OPENFOAMESIORFOUNDATION
    polyTopoChange& meshMod
#else
    directTopoChange& meshMod
#endif
)
{
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
#ifdef OPENFOAMFOUNDATION
    const meshFaceZones& faceZones = mesh.faceZones();
#else
    const faceZoneMesh& faceZones = mesh.faceZones();
#endif

    forAll(duplicates, bFaceI)
    {
        label otherFaceI = duplicates[bFaceI];

        if (otherFaceI != -1 && otherFaceI > bFaceI)
        {
            // Two duplicate faces. Merge.

            label face0 = mesh.nInternalFaces() + bFaceI;
            label face1 = mesh.nInternalFaces() + otherFaceI;

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        own1,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }
    }
}


void mergeDuplicateFaces(polyMesh& mesh)
{
    // Collect all boundary faces
    labelList boundaryFaces(mesh.nFaces() - mesh.nInternalFaces());
    forAll(boundaryFaces, i)
    {
        boundaryFaces[i] = i + mesh.nInternalFaces();
    }

    // Mesh change engine
#ifdef OPENFOAMESIORFOUNDATION
    polyTopoChange meshMod(mesh);
#else
    directTopoChange meshMod(mesh);
#endif

    // Get all duplicate face labels (in boundaryFaces indices!)
    labelList duplicates(findBaffles(mesh, boundaryFaces));

    // Merge into internal faces.
    insertDuplicateMerge(mesh, duplicates, meshMod);

    // Change the mesh. No inflation.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("Abaqus input file");

#   include "setRootCase.H"
#   include "createTime.H"

    // Read abaqus input file name
#ifdef OPENFOAMESIORFOUNDATION
    const fileName inputFileName(args[1]);
#else
    const fileName inputFileName(args.additionalArgs()[0]);
#endif
    Info<< "Reading " << inputFileName << nl << endl;

    // Read abaqus input file
    autoPtr<IFstream> inputStreamPtr(new IFstream(inputFileName));

    if (!inputStreamPtr())
    {
        FatalErrorIn(args.executable())
            << ": file " << inputFileName << " not found"
            << exit(FatalError);
    }

    // Check the header was read
    if (!readHeader(inputStreamPtr()))
    {
        FatalErrorIn(args.executable())
            << ": did not find the file header!"
            << exit(FatalError);
    }

    // Count the nodes and set the threeD flag
    bool threeD = false;
    Info<< "Counting the nodes" << endl;
    const label nNodes = countNodes(inputStreamPtr(), threeD);
    Info<< "Number of nodes: " << nNodes << endl;

    // Re-read the input file
    inputStreamPtr.clear();
    inputStreamPtr.set(new IFstream(inputFileName));

    // Re-read the header
    readHeader(inputStreamPtr());

    // Read the nodes
    Info<< "Reading the nodes" << endl;
    pointField points(nNodes, vector::zero);
    readNodes(inputStreamPtr(), nNodes, points, threeD);

    // Count elements
    Info<< "Counting the elements" << endl;
    const label nElements = countElements(inputStreamPtr(), threeD);
    Info<< "Number of elements: " << nElements << endl;

    // Re-read the input file
    inputStreamPtr.clear();
    inputStreamPtr.set(new IFstream(inputFileName));

    // Re-read the header
    readHeader(inputStreamPtr());

    // Re-count the nodes to move the stream to the elements
    {
        bool tmp = false;
        countNodes(inputStreamPtr(), tmp);
    }

    // Read the elements
    // We will remove the duplicate faces after
    Info<< "Reading the elements" << endl;
    faceList faces(nElements*6, face(4));
    cellList cells(nElements);
    readElements(inputStreamPtr(), nElements, faces, cells);

    // Create polyMesh
    // Note that the neighbour list will have length zero as every cell defines
    // its own faces
    // These duplicate faces are fixed after the patches are added, by using
    // the method from the mergeOrSplitBaffles utility
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ
        ),
#ifdef OPENFOAMESIORFOUNDATION
        std::move(points),
        std::move(faces),
        std::move(cells),
#else
        Xfer<pointField>(points),
        Xfer<faceList>(faces),
        Xfer<cellList>(cells),
#endif
        false
    );

    // Add boundary patches

    // void addPatches
    // (
    //     const List<polyPatch*>&,
    //     const bool validBoundary = true
    // );

    // Create new boundary patches
    Info<< "Creating boundary" << endl;
    DynamicList<polyPatch*> p(1);

    // Add all faces to one patch
    p.append
    (
        polyPatch::New
        (
            "patch",              // patch type
            "defaultPatch",       // name
            mesh.nFaces(),        // size
            0,                    // start face index
            0,                    // patch index
            mesh.boundaryMesh()   // polyBoundaryMesh
        ).ptr()
    );

    p.shrink();
    mesh.removeBoundary();
    mesh.addPatches(p);

    // Optionally, we could add all sets or surfaces are zones
    // Leave for another day, or until I really need it!
    // Add mesh zones
    // void addZones
    // (
    //     const List<pointZone*>& pz,
    //     const List<faceZone*>& fz,
    //     const List<cellZone*>& cz
    // );

    // Finally, we will remove all duplicate faces
    // This will update the boundary and result in a non-zero neighbour list
    // The procedure below is taken from the mergeOrSplitBaffles utility
    Info<< "Removing duplicate faces" << endl;
    mergeDuplicateFaces(mesh);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    // Write the mesh
    Info<< "Writing mesh to " << mesh.instance() << endl;
    mesh.setInstance("constant");
    mesh.write();

    Info<< nl << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
