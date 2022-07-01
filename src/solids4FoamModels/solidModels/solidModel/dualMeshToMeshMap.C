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

\*---------------------------------------------------------------------------*/

#include "dualMeshToMeshMap.H"
#include "polyMesh.H"
#include "Random.H"
#include "pointInNonConvexCell.H"
#ifdef FOAMEXTEND
    #include "tetPolyMesh.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dualMeshToMeshMap, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualMeshToMeshMap::dualMeshToMeshMap
(
    const polyMesh& mesh,
    const polyMesh& dualMesh,
    const meshDualiser& dualiser
)
:
    mesh_(mesh),
    dualMesh_(dualMesh),
    pointToDualCells_(dualiser.pointToDualCells()),
    pointToDualCell_(pointToDualCells_.size(), -1),
    pointToDualFaces_(pointToDualCells_.size()),
    pointToDualPoint_(dualiser.pointToDualPoint()),
    cellToDualPoint_(dualiser.cellToDualPoint()),
    faceToDualPoint_(dualiser.faceToDualPoint()),
    edgeToDualPoint_(dualiser.edgeToDualPoint()),
    dualFaceToCell_(dualMesh.nFaces(), -1),
    dualCellToPoint_(dualMesh.nCells(), -1)
{
    // Create pointToDualCell from pointToDualCells
    forAll(pointToDualCells_, pointI)
    {
        if (pointToDualCells_[pointI].size() != 1)
        {
            FatalErrorIn("dualMeshToMeshMap::dualMeshToMeshMap(...)")
                << "pointToDualCells_[" << pointI << "].size() != 1"
                << " for point = " << mesh.points()[pointI]
                << abort(FatalError);
        }

        pointToDualCell_[pointI] = pointToDualCells_[pointI][0];
    }

    // Calculate the minimum edge length, which we will use as a reference
    // length
    scalar minEdgeLength = GREAT;
    const edgeList& dualEdges = dualMesh.edges();
    const pointField& dualPoints = dualMesh.points();
    forAll(dualEdges, dualEdgeI)
    {
        const edge& curEdge = dualEdges[dualEdgeI];
        const scalar curEdgeLength = curEdge.mag(dualPoints);

        if (curEdgeLength < minEdgeLength)
        {
            minEdgeLength = curEdgeLength;
        }
    }

    // Construct pointToDualFaces from pointToDualCell
    //const faceList& dualFaces = dualMesh.faces();
    const cellList& dualCells = dualMesh.cells();
    forAll(pointToDualFaces_, pointI)
    {
        labelHashSet curDualFaces;

        const label dualCellID = pointToDualCell_[pointI];
        const cell& curDualCell = dualCells[dualCellID];

        // Add all faces in the dual cells
        // This must be taken into account when using this map
        pointToDualFaces_[pointI] = labelList(curDualCell);
    }

    // Set dualFaceToCell
    // All internal dual faces should uniquely lie within one primary cell. For
    // boundary faces, each dual face will uniquely lie on the boundary of one
    // cell. Note that the dual mesh and primary mesh boundaries may not exactly
    // overlap, so we must take care when using geometric checks, e.g.
    // pointInCell.
    // To avoid an expensive global search, we will use the pointToDualCellMap:
    // we will match the dual faces within each dual cell with the point cells
    // about each point, where we will use the pointInCell function.
    // For boundary dual faces, the pointInCell function may fail as the dual
    // face centre may not lie inside any primary mesh cell. Mostly, there is a
    // one-to-one map between boundary dual faces and boundary primary points,
    // but it is also possible for boundary dual points to directly lie on
    // primary mesh points.

    // Approach
    // For all points (pointI), get the corresponding dual cell (dualCellI)
    // Loop over the dual faces (dualFaceID) of that dual cell
    // If the dualFaceToCell map is already set for this dual face: skip
    // If pointCell[pointI].size == 1, then we found the cell: continue
    // If the dual face is internal:
    //     Loop over the pointCells (cellID) of pointI
    //     Use the pointInCell function to find which cell the dual face centre
    //     lies within
    // Else if the dual face is on the boundary:
    //     Loop over all boundary pointFaces (faceID) of pointI
    //     Project the dualFace centre to find which pointFace it hits
    //     The faceCell of the pointFace is correct cell: done
    //     If the hit misses, then no map is found
    const labelListList& pointCells = mesh.pointCells();
    const labelListList& pointFaces = mesh.pointFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const vectorField& dualFaceCentres = dualMesh.faceCentres();
    const vectorField& CI = mesh.cellCentres();

    // This is only needed in foam extend: OF can use
    // polyMesh::pointInCell(CELL_TETS)
#ifdef FOAMEXTEND
    tetPolyMesh tetMesh(mesh);
    const pointField tetPoints = tetMesh.points();
#endif

    forAll(pointToDualCell_, pointI)
    {
        const label dCellID = pointToDualCell_[pointI];
        const cell& dCell = dualCells[dCellID];

        const labelList& curPointCells = pointCells[pointI];

        if (debug)
        {
            Info<< nl << "    dCell = " << dualMesh.cellCentres()[dCellID] << nl
                << "    dCell nFaces = " << dCell.size() << endl;
        }

        forAll(dCell, fI)
        {
            const label dFaceID = dCell[fI];

            if (debug)
            {
                Info<< "    dFaceID = " << dFaceID << nl
                    << "    face centre = " << dualFaceCentres[dFaceID] << endl;
            }

            if (dualFaceToCell_[dFaceID] != -1)
            {
                // Already set: skip
                if (debug)
                {
                    Info<< "    already set" << endl;
                }
                continue;
            }

            // If there is only one point cell then the current dual face must
            // map to it
            if (curPointCells.size() == 1)
            {
                const label cellID = curPointCells[0];
                dualFaceToCell_[dFaceID] = cellID;
                continue;
            }

            if (dualMesh.isInternalFace(dFaceID))
            {
                // Dual face centre
                // We apply a small perturbation as pointInCell can fail when
                // the point is exactly on the face (e.g. of a decomposed tetra)
#ifdef FOAMEXTEND
                vector perturb = Random(0).vector01();
#else
                vector perturb = Random(0).sample01<vector>();
#endif
                const vector dFaceC =
                    dualFaceCentres[dFaceID]
                  + 0.0001*minEdgeLength*perturb/mag(perturb);

                // Check if the centre of the dual face dFaceID is within one of the
                // point cells
                forAll(curPointCells, pcI)
                {
                    const label cellID = curPointCells[pcI];

                    if (debug)
                    {
                        Info<< "    point cell = " << mesh.cellCentres()[cellID] << endl;
                    }

                    if
                    (
#ifdef FOAMEXTEND
                        pointInCellTetDecomp
                        (
                            mesh,
                            tetMesh.tets(cellID),
                            tetPoints,
                            dFaceC,
                            cellID
                        )
#else
                        mesh.polyMesh::pointInCell
                        (
                            dFaceC, cellID, polyMesh::CELL_TETS
                        )
#endif
                    )
                    {
                        if (debug)
                        {
                            Info<< "    found" << endl;
                        }
                        dualFaceToCell_[dFaceID] = cellID;
                        break;
                    }
                }

                // If we didn't find it then we will pick the closest cell of
                // the point cells. This is not ideal, but how bad can it be...
                // if (dualFaceToCell_[dFaceID] == -1)
                // {
                //     label closestCellID = -1;
                //     scalar dist = GREAT;
                //     forAll(curPointCells, pcI)
                //     {
                //         const label cellID = curPointCells[pcI];
                //         const vector& curC = CI[cellID];
                //         const scalar newDist = mag(dFaceC - curC);

                //         if (newDist < dist)
                //         {
                //             dist = newDist;
                //             closestCellID = cellID;
                //         }
                //     }

                //     dualFaceToCell_[dFaceID] = closestCellID;
                // }
            }
            else // boundary dual face
            {
                // Point faces for pointI
                const labelList& curPointFaces = pointFaces[pointI];

                // Dual face centre
                const vector dFaceC = dualFaceCentres[dFaceID];

                // Dual patch index
                const label dPatchID =
                    dualMesh.boundaryMesh().whichPatch(dFaceID);

                // Dual patch
                const polyPatch& dPp = dualMesh.boundaryMesh()[dPatchID];

                // Dual local face index
                const label dFaceLocalID = dFaceID - dPp.start();

                // Dual patch normal for the local dual face
                const vector dN = dPp.faceNormals()[dFaceLocalID];

                // Loop over all point faces of pointI and find which one the
                // dualFace projects on to
                forAll(curPointFaces, pfI)
                {
                    const label faceID = curPointFaces[pfI];

                    // Only loop over boundary faces
                    if (!mesh.isInternalFace(faceID))
                    {
                        // Get face
                        const face& curFace = faces[faceID];

                        // Project dual face centre to this face using the dual
                        // face normal (we could have also used the face normal)
                        if (curFace.ray(dFaceC, dN, points).hit())
                        {
                            // The dual face projects to this face, so the
                            // dual face must map to the faceCell of this face
                            const label cellID = faceOwner[faceID];

                            dualFaceToCell_[dFaceID] = cellID;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Check all were set correctly
    if (gMin(dualFaceToCell_) == -1)
    {
        label countIntFace = 0;
        label countBndFace = 0;
        forAll(dualFaceToCell_, dfI)
        {
            if (dualFaceToCell_[dfI] == -1)
            {
                if (dualMesh.isInternalFace(dfI))
                {
                    countIntFace++;
                }
                else
                {
                    countBndFace++;
                }
            }
        }
        reduce(countIntFace, sumOp<label>());
        reduce(countBndFace, sumOp<label>());

        Info<< "Constructing dualFaceToCell:" << nl
            << "    internal faces not set: " << countIntFace << nl
            << "    boundary faces not set: " << countBndFace << nl
            << endl;

        // It is an error for internal faces to be not set
        if (countIntFace > 1)
        {
            Info<< nl << "The dualFaceToCell map was not set for at least one "
                << "internal face!" << nl
                << "Details on up to the first 10 internal faces that were"
                << " not set:" << endl;
            int i = 0;
            forAll(dualFaceToCell_, dfI)
            {
                if (dualFaceToCell_[dfI] == -1)
                {
                    if (dualMesh.isInternalFace(dfI))
                    {
                        Pout<< "dualFace " << dfI
                            << ",  " << dualMesh.faceCentres()[dfI] << endl;

                        if (i++ == 10)
                        {
                            break;
                        }
                    }
                }
            }

            FatalErrorIn("dualMeshToMeshMap::dualMeshToMeshMap(...)")
                << "Problem setting dualFaceToCell map"
                << abort(FatalError);
        }
    }


    // Set dualCellToPoint
    // This is the reverse of the pointToDualCell map

    forAll(pointToDualCell_, pointI)
    {
        const label dCellID = pointToDualCell_[pointI];

        dualCellToPoint_[dCellID] = pointI;
    }

    // Check all were set correctly
    if (gMin(dualCellToPoint_) == -1)
    {
        label count = 0;
        forAll(dualCellToPoint_, dcI)
        {
            if (dualCellToPoint_[dcI] == -1)
            {
                count++;
            }
        }

        FatalErrorIn("dualMeshToMeshMap::dualMeshToMeshMap(...)")
            << "gMin(dualCellToPoint_) == -1: problem setting dualCellToPoint map"
            << nl << "There are " << count << " cells not set"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dualMeshToMeshMap::~dualMeshToMeshMap()
{}


// ************************************************************************* //
