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

#include "domainDecomposition.H"
#include "foamTime.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "globalMeshData.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

domainDecomposition::domainDecomposition(const IOobject& io)
:
    fvMesh(io),
    decompositionDict_
    (
        IOobject
        (
            "decomposeParDict",
            time().system(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nProcs_(readInt(decompositionDict_.lookup("numberOfSubdomains"))),
    distributed_(false),
    cellToProc_(nCells()),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    nInternalProcFaces_(nProcs_),
    nLiveProcFaces_(nProcs_),
    procCellAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    globallySharedPoints_(0),
    cyclicParallel_(false)
{
    if (decompositionDict_.found("distributed"))
    {
        distributed_ = Switch(decompositionDict_.lookup("distributed"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

domainDecomposition::~domainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool domainDecomposition::writeDecomposition()
{
    Info<< "\nConstructing processor meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(2*globallySharedPoints_.size());

    forAll (globallySharedPoints_, pointi)
    {
        sharedPointLookup.insert(globallySharedPoints_[pointi], pointi);
    }


    // Mark point/faces/cells that are in zones.  Bad coding - removed
    // HJ, 31/Mar/2009

    label totProcFaces = 0;
    label maxProcPatches = 0;
    label maxProcFaces = 0;


    // Write out the meshes
    for (label procI = 0; procI < nProcs_; procI++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[procI];

        // Access list of all points in the mesh.  HJ, 27/Mar/2009
        const pointField& meshPoints = allPoints();

        labelList pointLookup(meshPoints.size(), -1);

        pointField procPoints(curPointLabels.size());

        forAll (curPointLabels, pointi)
        {
            procPoints[pointi] = meshPoints[curPointLabels[pointi]];

            pointLookup[curPointLabels[pointi]] = pointi;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[procI];

        // Access list of all faces in the mesh.  HJ, 27/Mar/2009
        const faceList& meshFaces = allFaces();

        labelList faceLookup(meshFaces.size(), -1);

        faceList procFaces(curFaceLabels.size());

        forAll (curFaceLabels, facei)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            // HJ, 5/Dec/2001
            label curF = mag(curFaceLabels[facei]) - 1;

            faceLookup[curF] = facei;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[facei] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[facei];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll (origFaceLabels, pointi)
            {
                procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
            }
        }

        // Create cell lookup
        labelList cellLookup(nCells(), -1);
        const labelList& curCellLabels = procCellAddressing_[procI];

        forAll (curCellLabels, cellI)
        {
            cellLookup[curCellLabels[cellI]] = cellI;
        }

        // Get complete owner-neighour addressing in the mesh
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        // Calculate owner and neighbour list
        // Owner list is sized to number of live faces
        // Neighbour list is sized to number of internal faces

        labelList procOwner(nLiveProcFaces_[procI]);

        // Note: loop over owner, not all faces: sizes are different
        forAll (procOwner, faceI)
        {
            // Remember to decrement the index by one (turning index)
            // HJ, 28/Mar/2009
            label curF = mag(curFaceLabels[faceI]) - 1;

            if (curFaceLabels[faceI] >= 0)
            {
                procOwner[faceI] = cellLookup[own[curF]];
            }
            else
            {
                procOwner[faceI] = cellLookup[nei[curF]];
            }
        }

        labelList procNeighbour(nInternalProcFaces_[procI]);

        // Note: loop over neighbour, not all faces: sizes are different
        forAll (procNeighbour, faceI)
        {
            // Remember to decrement the index by one (turning index)
            // HJ, 28/Mar/2009
            label curF = mag(curFaceLabels[faceI]) - 1;

            if (curFaceLabels[faceI] >= 0)
            {
                procNeighbour[faceI] = cellLookup[nei[curF]];
            }
            else
            {
                procNeighbour[faceI] = cellLookup[own[curF]];
            }
        }

        // Create processor cells.  No longer needed: using owner and neighbour
        // HJ, 28/Mar/2009
//         const cellList& meshCells = cells();

//         cellList procCells(curCellLabels.size());

//         forAll (curCellLabels, cellI)
//         {
//             const labelList& origCellLabels = meshCells[curCellLabels[cellI]];

//             cell& curCell = procCells[cellI];

//             curCell.setSize(origCellLabels.size());

//             forAll (origCellLabels, cellFaceI)
//             {
//                 curCell[cellFaceI] = faceLookup[origCellLabels[cellFaceI]];
//             }
//         }

        // Create processor mesh without a boundary

        fileName processorCasePath
        (
            time().caseName()/fileName(word("processor") + Foam::name(procI))
        );

        // make the processor directory
        mkDir(time().rootPath()/processorCasePath);

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            processorCasePath,
            "system",
            "constant",
            true
        );

        // Create the mesh
        polyMesh procMesh
        (
            IOobject
            (
                this->polyMesh::name(),  // region name of undecomposed mesh
                pointsInstance(),
                processorDb
            ),
            xferMove(procPoints),
            xferMove(procFaces),
            xferMove(procOwner),
            xferMove(procNeighbour),
            false          // Do not sync par
//   xferMove(procCells)   // Old-fashioned mesh creation using cells.
                           // Deprecated: using face owner/neighbour
                           // HJ, 30/Mar/2009
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

        const labelList& curPatchSizes = procPatchSize_[procI];

        const labelList& curPatchStarts = procPatchStartIndex_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex_[procI];

        const polyPatchList& meshPatches = boundaryMesh();

        List<polyPatch*> procPatches
        (
            curPatchSizes.size()
          + curProcessorPatchSizes.size(),
            reinterpret_cast<polyPatch*>(0)
        );

        label nPatches = 0;

        forAll (curPatchSizes, patchi)
        {
            procPatches[nPatches] =
                meshPatches[curBoundaryAddressing[patchi]].clone
                (
                    procMesh.boundaryMesh(),
                    nPatches,
                    curPatchSizes[patchi],
                    curPatchStarts[patchi]
                ).ptr();

            nPatches++;
        }

        forAll (curProcessorPatchSizes, procPatchI)
        {
            procPatches[nPatches] =
                new processorPolyPatch
                (
                    word("procBoundary") + Foam::name(procI)
                  + word("to")
                  + Foam::name(curNeighbourProcessors[procPatchI]),
                    curProcessorPatchSizes[procPatchI],
                    curProcessorPatchStarts[procPatchI],
                    nPatches,
                    procMesh.boundaryMesh(),
                    procI,
                    curNeighbourProcessors[procPatchI]
            );

            nPatches++;
        }

        // Add boundary patches
        procMesh.addPatches(procPatches);

        // Create and add zones

        // Note:
        // This coding was all wrong, as each point/face/cell may only
        //  belong to a single zone.
        // Additionally, ordering of points/faces/cells in the processor mesh
        // needs to match the ordering in global mesh zones.  Full rewrite.
        // HJ, 30/Mar/2009

        // Create zones if needed
        if
        (
            pointZones().size() > 0
         || faceZones().size() > 0
         || cellZones().size() > 0
        )
        {
            // Point zones
            List<pointZone*> procPz(pointZones().size());

            {
                const pointZoneMesh& pz = pointZones();

                // Go through all the zoned points and find out if they
                // belong to a processor.  If so, add it to the zone as
                // necessary
                forAll (pz, zoneI)
                {
                    const labelList& zonePoints = pz[zoneI];

                    labelList procZonePoints(zonePoints.size());
                    label nZonePoints = 0;

                    forAll (zonePoints, pointI)
                    {
                        const label localIndex =
                            pointLookup[zonePoints[pointI]];

                        if (localIndex >= 0)
                        {
                            // Point live on processor: add to zone
                            procZonePoints[nZonePoints] = localIndex;
                            nZonePoints++;
                        }
                    }

                    // Add the zone
                    procZonePoints.setSize(nZonePoints);

                    procPz[zoneI] = new pointZone
                    (
                        pz[zoneI].name(),
                        procZonePoints,
                        zoneI,
                        procMesh.pointZones()
                    );
                }
            }


            // Face zones
            List<faceZone*> procFz(faceZones().size());

            {
                const faceZoneMesh& fz = faceZones();

                forAll (fz, zoneI)
                {
                    const labelList& zoneFaces = fz[zoneI];
                    const boolList& flipMap = fz[zoneI].flipMap();

                    // Go through all the zoned faces and find out if they
                    // belong to a processor.  If so, add it to the zone as
                    // necessary

                    labelList procZoneFaces(zoneFaces.size());
                    boolList procZoneFaceFlips(zoneFaces.size());
                    label nZoneFaces = 0;

                    forAll (zoneFaces, faceI)
                    {
                        const label localIndex = faceLookup[zoneFaces[faceI]];

                        if (localIndex >= 0)
                        {
                            // Face is present on the processor

                            // Add the face to the zone
                            procZoneFaces[nZoneFaces] = localIndex;

                            // Grab the flip
                            bool flip = flipMap[faceI];

                            if (curFaceLabels[localIndex] < 0)
                            {
                                flip = !flip;
                            }

                            procZoneFaceFlips[nZoneFaces] = flip;
                            nZoneFaces++;
                        }
                    }

                    // Add the zone
                    procZoneFaces.setSize(nZoneFaces);
                    procZoneFaceFlips.setSize(nZoneFaces);

                    procFz[zoneI] = new faceZone
                    (
                        fz[zoneI].name(),
                        procZoneFaces,
                        procZoneFaceFlips,
                        zoneI,
                        procMesh.faceZones()
                    );
                }
            }

            // Cell zones
            List<cellZone*> procCz(cellZones().size());

            {
                const cellZoneMesh& cz = cellZones();

                // Go through all the zoned cells and find out if they
                // belong to a processor.  If so, add it to the zone as
                // necessary

                forAll (cz, zoneI)
                {
                    const labelList& zoneCells = cz[zoneI];

                    labelList procZoneCells(zoneCells.size());
                    label nZoneCells = 0;

                    forAll (zoneCells, cellI)
                    {
                        const label localIndex = cellLookup[zoneCells[cellI]];

                        if (localIndex >= 0)
                        {
                            procZoneCells[nZoneCells] = localIndex;
                            nZoneCells++;
                        }
                    }

                    // Add the zone
                    procZoneCells.setSize(nZoneCells);

                    procCz[zoneI] = new cellZone
                    (
                        cz[zoneI].name(),
                        procZoneCells,
                        zoneI,
                        procMesh.cellZones()
                    );
                }
            }

            // Add zones
            procMesh.addZones(procPz, procFz, procCz);
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        procMesh.write();

        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << procMesh.nCells()
            << endl;

        label nBoundaryFaces = 0;
        label nProcPatches = 0;
        label nProcFaces = 0;

        forAll (procMesh.boundaryMesh(), patchi)
        {
            if
            (
                procMesh.boundaryMesh()[patchi].type()
             == processorPolyPatch::typeName
            )
            {
                const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>
                (
                    procMesh.boundaryMesh()[patchi]
                );

                Info<< "    Number of faces shared with processor "
                    << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                nProcPatches++;
                nProcFaces += ppp.size();
            }
            else
            {
                nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
            }
        }

        Info<< "    Number of processor patches = " << nProcPatches << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = " << nBoundaryFaces << endl;

        totProcFaces += nProcFaces;
        maxProcPatches = max(maxProcPatches, nProcPatches);
        maxProcFaces = max(maxProcFaces, nProcFaces);

        // create and write the addressing information
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procPointAddressing_[procI]
        );
        pointProcAddressing.write();

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procFaceAddressing_[procI]
        );
        faceProcAddressing.write();

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procCellAddressing_[procI]
        );
        cellProcAddressing.write();

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procBoundaryAddressing_[procI]
        );
        boundaryProcAddressing.write();
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of processor patches = " << maxProcPatches << nl
        << "Max number of faces between processors = " << maxProcFaces
        << endl;

    return true;
}


// ************************************************************************* //
