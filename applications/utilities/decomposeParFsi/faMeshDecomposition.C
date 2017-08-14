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

#include "faMeshDecomposition.H"
#include "foamTime.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorFaPatch.H"
#include "faMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "globalMeshData.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void faMeshDecomposition::distributeFaces()
{
    Info<< "\nCalculating distribution of faces" << endl;

    cpuTime decompositionTime;

    for (label procI = 0; procI < nProcs(); procI++)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/fileName(word("processor") + Foam::name(procI))
        );

        fvMesh procMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                processorDb.timeName(),
                processorDb
            )
        );

        labelHashSet faceProcAddressingHash
        (
            labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    "constant",
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        forAll (faceLabels(), faceI)
        {
            if (faceProcAddressingHash.found(faceLabels()[faceI] + 1))
            {
                faceToProc_[faceI] = procI;
            }
        }
    }

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
faMeshDecomposition::faMeshDecomposition(const fvMesh& mesh)
:
    faMesh(mesh),
    decompositionDict_
    (
        IOobject
        (
            "decomposeParDict",
            time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nProcs_(readInt(decompositionDict_.lookup("numberOfSubdomains"))),
    distributed_(false),
    faceToProc_(nFaces()),
    procFaceLabels_(nProcs_),
    procMeshEdgesMap_(nProcs_),
    procNInternalEdges_(nProcs_, 0),
    procPatchEdgeLabels_(nProcs_),
    procPatchPointAddressing_(nProcs_),
    procPatchEdgeAddressing_(nProcs_),
    procEdgeAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
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

faMeshDecomposition::~faMeshDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void faMeshDecomposition::decomposeMesh(const bool filterEmptyPatches)
{
    // Decide which cell goes to which processor
    distributeFaces();

    Info<< "\nDistributing faces to processors" << endl;

    // Memory management
    {
        List<SLList<label> > procFaceList(nProcs());

        forAll (faceToProc_, faceI)
        {
            if (faceToProc_[faceI] >= nProcs())
            {
                FatalErrorIn("Finite area mesh decomposition")
                    << "Impossible processor label " << faceToProc_[faceI]
                    << "for face " << faceI
                    << abort(FatalError);
            }
            else
            {
                procFaceList[faceToProc_[faceI]].append(faceI);
            }
        }

        // Convert linked lists into normal lists
        forAll (procFaceList, procI)
        {
            procFaceAddressing_[procI] = procFaceList[procI];
        }
    }


    // Find processor mesh faceLabels and ...

    for (label procI = 0; procI < nProcs(); procI++)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/fileName(word("processor") + Foam::name(procI))
        );

        fvMesh procFvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                processorDb.timeName(),
                processorDb
            )
        );

        labelIOList fvPointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                "constant",
                procFvMesh.meshSubDir,
                procFvMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        HashTable<label, label, Hash<label> > fvFaceProcAddressingHash;

        {
            labelIOList fvFaceProcAddressing
            (
                IOobject
                (
                    "faceProcAddressing",
                    "constant",
                    procFvMesh.meshSubDir,
                    procFvMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            forAll(fvFaceProcAddressing, faceI)
            {
                 fvFaceProcAddressingHash.insert
                 (
                     fvFaceProcAddressing[faceI], faceI
                 );
            }
        };

        const labelList& curProcFaceAddressing = procFaceAddressing_[procI];

        labelList& curFaceLabels = procFaceLabels_[procI];

        curFaceLabels = labelList(curProcFaceAddressing.size(), -1);

        forAll(curProcFaceAddressing, faceI)
        {
            curFaceLabels[faceI] =
                fvFaceProcAddressingHash.find
                (
                    faceLabels()[curProcFaceAddressing[faceI]] + 1
                )();
        }

        // create processor finite area mesh
        faMesh procMesh
        (
            procFvMesh,
            procFaceLabels_[procI]
        );

        const indirectPrimitivePatch& patch = this->patch();
        const Map<label>& map = patch.meshPointMap();

        HashTable<label, edge, Hash<edge> > edgesHash;

        label edgeI = -1;

        label nIntEdges = patch.nInternalEdges();

        for (label curEdge = 0; curEdge < nIntEdges; curEdge++)
        {
            edgesHash.insert(patch.edges()[curEdge], ++edgeI);
        }

        forAll (boundary(), patchI)
        {
            // Include emptyFaPatch

            label size = boundary()[patchI].labelList::size();

            for(int eI=0; eI<size; eI++)
            {
                edgesHash.insert(patch.edges()[boundary()[patchI][eI]], ++edgeI);
            }
        }


        const indirectPrimitivePatch& procPatch = procMesh.patch();
        const vectorField& procPoints = procPatch.localPoints();
        const labelList& procMeshPoints = procPatch.meshPoints();
        const edgeList& procEdges = procPatch.edges();

        labelList& curPatchPointAddressing = procPatchPointAddressing_[procI];
        curPatchPointAddressing.setSize(procPoints.size(), -1);

        forAll(procPoints, pointI)
        {
            curPatchPointAddressing[pointI] =
                map[fvPointProcAddressing[procMeshPoints[pointI]]];
        }

        labelList& curPatchEdgeAddressing = procPatchEdgeAddressing_[procI];
        curPatchEdgeAddressing.setSize(procEdges.size(), -1);

        forAll(procEdges, edgeI)
        {
            edge curGlobalEdge = procEdges[edgeI];
            curGlobalEdge[0] = curPatchPointAddressing[curGlobalEdge[0]];
            curGlobalEdge[1] = curPatchPointAddressing[curGlobalEdge[1]];

            curPatchEdgeAddressing[edgeI] = edgesHash.find(curGlobalEdge)();
        }

        Map<label>& curMap = procMeshEdgesMap_[procI];

        curMap = Map<label>(2*procEdges.size());

        forAll(curPatchEdgeAddressing, edgeI)
        {
            curMap.insert(curPatchEdgeAddressing[edgeI], edgeI);
        }

        procNInternalEdges_[procI] = procPatch.nInternalEdges();
    }


    Info << "\nDistributing edges to processors" << endl;

    // Loop through all internal edges and decide which processor they
    // belong to. First visit all internal edges.

    // set references to the original mesh
    const faBoundaryMesh& patches = boundary();
    const edgeList& edges = this->edges();
    const labelList& owner = edgeOwner();
    const labelList& neighbour = edgeNeighbour();

    // Memory management
    {
        List<SLList<label> > procEdgeList(nProcs());

        forAll(procEdgeList, procI)
        {
            for(label i=0; i<procNInternalEdges_[procI]; i++)
            {
                procEdgeList[procI].append(procPatchEdgeAddressing_[procI][i]);
            }
        }


        // Detect inter-processor boundaries

        // Neighbour processor for each subdomain
        List<SLList<label> > interProcBoundaries(nProcs());

        // Edge labels belonging to each inter-processor boundary
        List<SLList<SLList<label> > > interProcBEdges(nProcs());

        List<SLList<label> > procPatchIndex(nProcs());

        forAll (neighbour, edgeI)
        {
            if (faceToProc_[owner[edgeI]] != faceToProc_[neighbour[edgeI]])
            {
                // inter - processor patch edge found. Go through the list of
                // inside boundaries for the owner processor and try to find
                // this inter-processor patch.

                label ownerProc = faceToProc_[owner[edgeI]];
                label neighbourProc = faceToProc_[neighbour[edgeI]];

                SLList<label>::iterator curInterProcBdrsOwnIter =
                    interProcBoundaries[ownerProc].begin();

                SLList<SLList<label> >::iterator curInterProcBEdgesOwnIter =
                    interProcBEdges[ownerProc].begin();

                bool interProcBouFound = false;

                // WARNING: Synchronous SLList iterators

                for
                (
                    ;
                    curInterProcBdrsOwnIter
                 != interProcBoundaries[ownerProc].end()
                 && curInterProcBEdgesOwnIter
                 != interProcBEdges[ownerProc].end();
                    ++curInterProcBdrsOwnIter, ++curInterProcBEdgesOwnIter
                )
                {
                    if (curInterProcBdrsOwnIter() == neighbourProc)
                    {
                        // the inter - processor boundary exists. Add the face
                        interProcBouFound = true;

                        curInterProcBEdgesOwnIter().append(edgeI);

                        SLList<label>::iterator curInterProcBdrsNeiIter =
                            interProcBoundaries[neighbourProc].begin();

                        SLList<SLList<label> >::iterator
                            curInterProcBEdgesNeiIter =
                            interProcBEdges[neighbourProc].begin();

                        bool neighbourFound = false;

                        // WARNING: Synchronous SLList iterators

                        for
                        (
                            ;
                            curInterProcBdrsNeiIter !=
                            interProcBoundaries[neighbourProc].end()
                         && curInterProcBEdgesNeiIter !=
                            interProcBEdges[neighbourProc].end();
                            ++curInterProcBdrsNeiIter,
                            ++curInterProcBEdgesNeiIter
                        )
                        {
                            if (curInterProcBdrsNeiIter() == ownerProc)
                            {
                                // boundary found. Add the face
                                neighbourFound = true;

                                curInterProcBEdgesNeiIter().append(edgeI);
                            }

                            if (neighbourFound) break;
                        }

                        if (interProcBouFound && !neighbourFound)
                        {
                            FatalErrorIn
                                ("faDomainDecomposition::decomposeMesh()")
                                << "Inconsistency in inter - "
                                << "processor boundary lists for processors "
                                << ownerProc << " and " << neighbourProc
                                << abort(FatalError);
                        }
                    }

                    if (interProcBouFound) break;
                }

                if (!interProcBouFound)
                {
                    // inter - processor boundaries do not exist and need to
                    // be created

                    // set the new addressing information

                    // owner
                    interProcBoundaries[ownerProc].append(neighbourProc);
                    interProcBEdges[ownerProc].append(SLList<label>(edgeI));

                    // neighbour
                    interProcBoundaries[neighbourProc].append(ownerProc);
                    interProcBEdges[neighbourProc].append
                        (
                            SLList<label>(edgeI)
                        );
                }
            }
        }


        // Loop through patches. For cyclic boundaries detect inter-processor
        // edges; for all other, add edges to the edge list and remember start
        // and size of all patches.

        // for all processors, set the size of start index and patch size
        // lists to the number of patches in the mesh
        forAll (procPatchSize_, procI)
        {
            procPatchSize_[procI].setSize(patches.size());
            procPatchStartIndex_[procI].setSize(patches.size());
        }

        forAll (patches, patchI)
        {
            // Reset size and start index for all processors
            forAll (procPatchSize_, procI)
            {
                procPatchSize_[procI][patchI] = 0;
                procPatchStartIndex_[procI][patchI] =
                    procEdgeList[procI].size();
            }

            const label patchStart = patches[patchI].start();

//             if (typeid(patches[patchI]) != typeid(cyclicFaPatch))
            if (true)
            {
                // Normal patch. Add edges to processor where the face
                // next to the edge lives

                const labelListList& eF = patch().edgeFaces();

                label size = patches[patchI].labelList::size();

                labelList patchEdgeFaces(size, -1);

                for(int eI=0; eI<size; eI++)
                {
                    patchEdgeFaces[eI] = eF[patches[patchI][eI]][0];
                }

                forAll (patchEdgeFaces, edgeI)
                {
                    const label curProc = faceToProc_[patchEdgeFaces[edgeI]];

                    // add the face
                    procEdgeList[curProc].append(patchStart + edgeI);

                    // increment the number of edges for this patch
                    procPatchSize_[curProc][patchI]++;
                }
            }
            else
            {
                // Cyclic patch special treatment

                const faPatch& cPatch = patches[patchI];

                const label cycOffset = cPatch.size()/2;

                // Set reference to faceCells for both patches
                const labelList::subList firstEdgeFaces
                (
                    cPatch.edgeFaces(),
                    cycOffset
                );

                const labelList::subList secondEdgeFaces
                (
                    cPatch.edgeFaces(),
                    cycOffset,
                    cycOffset
                );

                forAll (firstEdgeFaces, edgeI)
                {
                    if
                    (
                        faceToProc_[firstEdgeFaces[edgeI]]
                     != faceToProc_[secondEdgeFaces[edgeI]]
                    )
                    {
                        // This edge becomes an inter-processor boundary edge
                        // inter - processor patch edge found. Go through
                        // the list of inside boundaries for the owner
                        // processor and try to find this inter-processor
                        // patch.

                        cyclicParallel_ = true;

                        label ownerProc = faceToProc_[firstEdgeFaces[edgeI]];
                        label neighbourProc =
                            faceToProc_[secondEdgeFaces[edgeI]];

                        SLList<label>::iterator curInterProcBdrsOwnIter =
                            interProcBoundaries[ownerProc].begin();

                        SLList<SLList<label> >::iterator
                            curInterProcBEdgesOwnIter =
                            interProcBEdges[ownerProc].begin();

                        bool interProcBouFound = false;

                        // WARNING: Synchronous SLList iterators

                        for
                        (
                            ;
                            curInterProcBdrsOwnIter !=
                            interProcBoundaries[ownerProc].end()
                         && curInterProcBEdgesOwnIter !=
                            interProcBEdges[ownerProc].end();
                            ++curInterProcBdrsOwnIter,
                            ++curInterProcBEdgesOwnIter
                        )
                        {
                            if (curInterProcBdrsOwnIter() == neighbourProc)
                            {
                                // the inter - processor boundary exists.
                                // Add the face
                                interProcBouFound = true;

                                curInterProcBEdgesOwnIter().append
                                    (patchStart + edgeI);

                                SLList<label>::iterator curInterProcBdrsNeiIter
                                   = interProcBoundaries[neighbourProc].begin();

                                SLList<SLList<label> >::iterator
                                    curInterProcBEdgesNeiIter =
                                    interProcBEdges[neighbourProc].begin();

                                bool neighbourFound = false;

                                // WARNING: Synchronous SLList iterators

                                for
                                (
                                    ;
                                    curInterProcBdrsNeiIter
                                   != interProcBoundaries[neighbourProc].end()
                                 && curInterProcBEdgesNeiIter
                                   != interProcBEdges[neighbourProc].end();
                                    ++curInterProcBdrsNeiIter,
                                    ++curInterProcBEdgesNeiIter
                                )
                                {
                                    if (curInterProcBdrsNeiIter() == ownerProc)
                                    {
                                        // boundary found. Add the face
                                        neighbourFound = true;

                                        curInterProcBEdgesNeiIter()
                                           .append
                                            (
                                                patchStart
                                              + cycOffset
                                              + edgeI
                                            );
                                    }

                                    if (neighbourFound) break;
                                }

                                if (interProcBouFound && !neighbourFound)
                                {
                                    FatalErrorIn
                                    (
                                        "faDomainDecomposition::decomposeMesh()"
                                    )   << "Inconsistency in inter-processor "
                                        << "boundary lists for processors "
                                        << ownerProc << " and " << neighbourProc
                                        << " in cyclic boundary matching"
                                        << abort(FatalError);
                                }
                            }

                            if (interProcBouFound) break;
                        }

                        if (!interProcBouFound)
                        {
                            // inter - processor boundaries do not exist
                            // and need to be created

                            // set the new addressing information

                            // owner
                            interProcBoundaries[ownerProc]
                                .append(neighbourProc);
                            interProcBEdges[ownerProc]
                                .append(SLList<label>(patchStart + edgeI));

                            // neighbour
                            interProcBoundaries[neighbourProc]
                               .append(ownerProc);
                            interProcBEdges[neighbourProc]
                               .append
                                (
                                    SLList<label>
                                    (
                                        patchStart
                                      + cycOffset
                                      + edgeI
                                    )
                                );
                        }
                    }
                    else
                    {
                        // This cyclic edge remains on the processor
                        label ownerProc = faceToProc_[firstEdgeFaces[edgeI]];

                        // add the edge
                        procEdgeList[ownerProc].append(patchStart + edgeI);

                        // increment the number of edges for this patch
                        procPatchSize_[ownerProc][patchI]++;

                        // Note: I cannot add the other side of the cyclic
                        // boundary here because this would violate the order.
                        // They will be added in a separate loop below
                    }
                }

                // Ordering in cyclic boundaries is important.
                // Add the other half of cyclic edges for cyclic boundaries
                // that remain on the processor
                forAll (secondEdgeFaces, edgeI)
                {
                    if
                    (
                        faceToProc_[firstEdgeFaces[edgeI]]
                     == faceToProc_[secondEdgeFaces[edgeI]]
                    )
                    {
                        // This cyclic edge remains on the processor
                        label ownerProc = faceToProc_[firstEdgeFaces[edgeI]];

                        // add the second edge
                        procEdgeList[ownerProc].append
                            (patchStart + cycOffset + edgeI);

                        // increment the number of edges for this patch
                        procPatchSize_[ownerProc][patchI]++;
                    }
                }
            }
        }

        // Convert linked lists into normal lists
        // Add inter-processor boundaries and remember start indices
        forAll (procEdgeList, procI)
        {
            // Get internal and regular boundary processor faces
            SLList<label>& curProcEdges = procEdgeList[procI];

            // Get reference to processor edge addressing
            labelList& curProcEdgeAddressing = procEdgeAddressing_[procI];

            labelList& curProcNeighbourProcessors =
                procNeighbourProcessors_[procI];

            labelList& curProcProcessorPatchSize =
                procProcessorPatchSize_[procI];

            labelList& curProcProcessorPatchStartIndex =
                procProcessorPatchStartIndex_[procI];

            // calculate the size
            label nEdgesOnProcessor = curProcEdges.size();

            for
            (
                SLList<SLList<label> >::iterator curInterProcBEdgesIter =
                    interProcBEdges[procI].begin();
                curInterProcBEdgesIter != interProcBEdges[procI].end();
                ++curInterProcBEdgesIter
            )
            {
                nEdgesOnProcessor += curInterProcBEdgesIter().size();
            }

            curProcEdgeAddressing.setSize(nEdgesOnProcessor);

            // Fill in the list. Calculate turning index.
            // Turning index will be -1 only for some edges on processor
            // boundaries, i.e. the ones where the current processor ID
            // is in the face which is a edge neighbour.
            // Turning index is stored as the sign of the edge addressing list

            label nEdges = 0;

            // Add internal and boundary edges
            // Remember to increment the index by one such that the
            // turning index works properly.
            for
            (
                SLList<label>::iterator curProcEdgeIter = curProcEdges.begin();
                curProcEdgeIter != curProcEdges.end();
                ++curProcEdgeIter
            )
            {
                curProcEdgeAddressing[nEdges] = curProcEdgeIter();
//                 curProcEdgeAddressing[nEdges] = curProcEdgeIter() + 1;
                nEdges++;
            }

            // Add inter-processor boundary edges. At the beginning of each
            // patch, grab the patch start index and size

            curProcNeighbourProcessors.setSize
            (
                interProcBoundaries[procI].size()
            );

            curProcProcessorPatchSize.setSize
            (
                interProcBoundaries[procI].size()
            );

            curProcProcessorPatchStartIndex.setSize
            (
                interProcBoundaries[procI].size()
            );

            label nProcPatches = 0;

            SLList<label>::iterator curInterProcBdrsIter =
                interProcBoundaries[procI].begin();

            SLList<SLList<label> >::iterator curInterProcBEdgesIter =
                interProcBEdges[procI].begin();

            for
            (
                ;
                curInterProcBdrsIter != interProcBoundaries[procI].end()
             && curInterProcBEdgesIter != interProcBEdges[procI].end();
                ++curInterProcBdrsIter, ++curInterProcBEdgesIter
            )
            {
                curProcNeighbourProcessors[nProcPatches] =
                    curInterProcBdrsIter();

                // Get start index for processor patch
                curProcProcessorPatchStartIndex[nProcPatches] = nEdges;

                label& curSize =
                    curProcProcessorPatchSize[nProcPatches];

                curSize = 0;

                // add faces for this processor boundary

                for
                (
                    SLList<label>::iterator curEdgesIter =
                        curInterProcBEdgesIter().begin();
                    curEdgesIter != curInterProcBEdgesIter().end();
                    ++curEdgesIter
                )
                {
                    // add the edges

                    // Remember to increment the index by one such that the
                    // turning index works properly.
                    if (faceToProc_[owner[curEdgesIter()]] == procI)
                    {
                        curProcEdgeAddressing[nEdges] = curEdgesIter();
//                         curProcEdgeAddressing[nEdges] = curEdgesIter() + 1;
                    }
                    else
                    {
                        // turning edge
                        curProcEdgeAddressing[nEdges] = curEdgesIter();
//                         curProcEdgeAddressing[nEdges] = -(curEdgesIter() + 1);
                    }

                    // increment the size
                    curSize++;

                    nEdges++;
                }

                nProcPatches++;
            }
        }
    }

    Info << "\nCalculating processor boundary addressing" << endl;
    // For every patch of processor boundary, find the index of the original
    // patch. Mis-alignment is caused by the fact that patches with zero size
    // are omitted. For processor patches, set index to -1.
    // At the same time, filter the procPatchSize_ and procPatchStartIndex_
    // lists to exclude zero-size patches
    forAll (procPatchSize_, procI)
    {
        // Make a local copy of old lists
        const labelList oldPatchSizes = procPatchSize_[procI];

        const labelList oldPatchStarts = procPatchStartIndex_[procI];

        labelList& curPatchSizes = procPatchSize_[procI];

        labelList& curPatchStarts = procPatchStartIndex_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

        curBoundaryAddressing.setSize
        (
            oldPatchSizes.size()
          + curProcessorPatchSizes.size()
        );

        label nPatches = 0;

        forAll (oldPatchSizes, patchI)
        {
            if (!filterEmptyPatches || oldPatchSizes[patchI] > 0)
            {
                curBoundaryAddressing[nPatches] = patchI;

                curPatchSizes[nPatches] = oldPatchSizes[patchI];

                curPatchStarts[nPatches] = oldPatchStarts[patchI];

                nPatches++;
            }
        }

        // reset to the size of live patches
        curPatchSizes.setSize(nPatches);
        curPatchStarts.setSize(nPatches);

        forAll (curProcessorPatchSizes, procPatchI)
        {
            curBoundaryAddressing[nPatches] = -1;

            nPatches++;
        }

        curBoundaryAddressing.setSize(nPatches);
    }


    // Gather data about globally shared points

    labelList globallySharedPoints_(0);

    // Memory management
    {
        labelList pointsUsage(nPoints(), 0);

        // Globally shared points are the ones used by more than 2 processors
        // Size the list approximately and gather the points
        labelHashSet gSharedPoints
        (
            min(100, nPoints()/1000)
        );

        // Loop through all the processors and mark up points used by
        // processor boundaries.  When a point is used twice, it is a
        // globally shared point

        for (label procI = 0; procI < nProcs(); procI++)
        {
            // Get list of edge labels
            const labelList& curEdgeLabels = procEdgeAddressing_[procI];

            // Get start of processor faces
            const labelList& curProcessorPatchStarts =
                procProcessorPatchStartIndex_[procI];

            const labelList& curProcessorPatchSizes =
                procProcessorPatchSize_[procI];

            // Reset the lookup list
            pointsUsage = 0;

            forAll (curProcessorPatchStarts, patchI)
            {
                const label curStart = curProcessorPatchStarts[patchI];
                const label curEnd = curStart + curProcessorPatchSizes[patchI];

                for
                (
                    label edgeI = curStart;
                    edgeI < curEnd;
                    edgeI++
                )
                {
                    // Mark the original edge as used
                    // Remember to decrement the index by one (turning index)
                    const label curE = curEdgeLabels[edgeI];

                    const edge& e = edges[curE];

                    forAll (e, pointI)
                    {
                        if (pointsUsage[e[pointI]] == 0)
                        {
                            // Point not previously used
                            pointsUsage[e[pointI]] = patchI + 1;
                        }
                        else if (pointsUsage[e[pointI]] != patchI + 1)
                        {
                            // Point used by some other patch = global point!
                            gSharedPoints.insert(e[pointI]);
                        }
                    }
                }
            }
        }

        // Grab the result from the hash list
        globallySharedPoints_ = gSharedPoints.toc();
        sort(globallySharedPoints_);
    }


    // Edge label for faPatches

    for (label procI = 0; procI < nProcs(); procI++)
    {
        fileName processorCasePath
        (
            time().caseName()/fileName(word("processor")
          + Foam::name(procI))
        );

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            processorCasePath
        );


        // read finite volume mesh
        fvMesh procFvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                processorDb.timeName(),
                processorDb
            )
        );

        // create finite area mesh
        faMesh procMesh
        (
            procFvMesh,
            procFaceLabels_[procI]
        );


        const labelList& curEdgeAddressing = procEdgeAddressing_[procI];

        const labelList& curPatchStartIndex = procPatchStartIndex_[procI];
        const labelList& curPatchSize = procPatchSize_[procI];

        const labelList& curProcessorPatchStartIndex =
            procProcessorPatchStartIndex_[procI];

        const labelList& curProcessorPatchSize =
            procProcessorPatchSize_[procI];

        labelListList& curPatchEdgeLabels = procPatchEdgeLabels_[procI];
        curPatchEdgeLabels =
            labelListList
            (
                curPatchSize.size()
              + curProcessorPatchSize.size()
            );

        forAll(curPatchSize, patchI)
        {
            labelList& curEdgeLabels = curPatchEdgeLabels[patchI];
            curEdgeLabels.setSize(curPatchSize[patchI], -1);

            label edgeI = 0;

            for
            (
                int i=curPatchStartIndex[patchI];
                i<(curPatchStartIndex[patchI]+curPatchSize[patchI]);
                i++
            )
            {
                curEdgeLabels[edgeI] =
                    procMeshEdgesMap_[procI][curEdgeAddressing[i]];
                edgeI++;
            }
        }

        forAll(curProcessorPatchSize, patchI)
        {
            labelList& curEdgeLabels  =
                curPatchEdgeLabels[curPatchSize.size() + patchI];
            curEdgeLabels.setSize(curProcessorPatchSize[patchI], -1);

            label edgeI = 0;

            for
            (
                int i=curProcessorPatchStartIndex[patchI];
                i<(curProcessorPatchStartIndex[patchI]
                +curProcessorPatchSize[patchI]);
                i++
            )
            {
                curEdgeLabels[edgeI] =
                    procMeshEdgesMap_[procI][curEdgeAddressing[i]];
                edgeI++;
            }
        }
    }
}


bool faMeshDecomposition::writeDecomposition()
{
    Info<< "\nConstructing processor FA meshes" << endl;


    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(2*globallySharedPoints_.size());

    forAll (globallySharedPoints_, pointi)
    {
        sharedPointLookup.insert(globallySharedPoints_[pointi], pointi);
    }

    label totProcEdges = 0;
    label maxProcPatches = 0;
    label maxProcEdges = 0;

    // Write out the meshes
    for (label procI = 0; procI < nProcs(); procI++)
    {
        // Create processor mesh without a boundary

        fileName processorCasePath
        (
            time().caseName()/fileName(word("processor") + Foam::name(procI))
        );

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            processorCasePath
        );

        // read finite volume mesh
        fvMesh procFvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                processorDb.timeName(),
                processorDb
            )
        );

        labelIOList fvBoundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                "constant",
                procFvMesh.meshSubDir,
                procFvMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );


        // create finite area mesh
        faMesh procMesh
        (
            procFvMesh,
            procFaceLabels_[procI]
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing =
            procBoundaryAddressing_[procI];

        const labelList& curPatchSizes = procPatchSize_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        const labelListList& curPatchEdgeLabels =
            procPatchEdgeLabels_[procI];

        const faPatchList& meshPatches = boundary();

        List<faPatch*> procPatches
        (
            curPatchSizes.size()
          + curProcessorPatchSizes.size(),
            reinterpret_cast<faPatch*>(NULL)
        );

        label nPatches = 0;

        forAll (curPatchSizes, patchi)
        {
            const labelList& curEdgeLabels = curPatchEdgeLabels[nPatches];

            label ngbPolyPatchIndex =
                findIndex
                (
                    fvBoundaryProcAddressing,
                    meshPatches[curBoundaryAddressing[patchi]]
                   .ngbPolyPatchIndex()
                );

            procPatches[nPatches] =
                meshPatches[curBoundaryAddressing[patchi]].clone
                (
                    procMesh.boundary(),
                    curEdgeLabels,
                    nPatches,
                    ngbPolyPatchIndex
                ).ptr();

            nPatches++;
        }

        forAll (curProcessorPatchSizes, procPatchI)
        {
            const labelList& curEdgeLabels = curPatchEdgeLabels[nPatches];

            procPatches[nPatches] =
                new processorFaPatch
                (
                    word("procBoundary") + Foam::name(procI)
                  + word("to")
                  + Foam::name(curNeighbourProcessors[procPatchI]),
                    curEdgeLabels,
                    nPatches,
                    procMesh.boundary(),
                    -1,
                    procI,
                    curNeighbourProcessors[procPatchI]
                );

            nPatches++;
        }

        // Add boundary patches
        procMesh.addFaPatches(procPatches);

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        procMesh.write();

        Info<< endl
            << "Processor " << procI << nl
            << "    Number of faces = " << procMesh.nFaces()
            << endl;

        label nBoundaryEdges = 0;
        label nProcPatches = 0;
        label nProcEdges = 0;

        forAll (procMesh.boundary(), patchi)
        {
            if
            (
                procMesh.boundary()[patchi].type()
             == processorFaPatch::typeName
            )
            {
                const processorFaPatch& ppp =
                refCast<const processorFaPatch>
                (
                    procMesh.boundary()[patchi]
                );

                Info<< "    Number of edges shared with processor "
                    << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                nProcPatches++;
                nProcEdges += ppp.size();
            }
            else
            {
                nBoundaryEdges += procMesh.boundary()[patchi].size();
            }
        }

        Info<< "    Number of processor patches = " << nProcPatches << nl
            << "    Number of processor edges = " << nProcEdges << nl
            << "    Number of boundary edges = " << nBoundaryEdges << endl;

        totProcEdges += nProcEdges;
        maxProcPatches = max(maxProcPatches, nProcPatches);
        maxProcEdges = max(maxProcEdges, nProcEdges);

        // create and write the addressing information
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procFvMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procPatchPointAddressing_[procI]
        );
        pointProcAddressing.write();

        labelIOList edgeProcAddressing
        (
            IOobject
            (
                "edgeProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procFvMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procEdgeAddressing_[procI]
        );
        edgeProcAddressing.write();

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procFvMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procFaceAddressing_[procI]
        );
        faceProcAddressing.write();

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procFvMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procBoundaryAddressing_[procI]
        );
        boundaryProcAddressing.write();
    }

    Info<< nl
        << "Number of processor edges = " << totProcEdges/2 << nl
        << "Max number of processor patches = " << maxProcPatches << nl
        << "Max number of faces between processors = " << maxProcEdges
        << endl;

    return true;
}
