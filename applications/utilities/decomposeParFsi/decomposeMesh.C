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

InClass
    domainDecomposition

Description
    Private member of domainDecomposition.
    Decomposes the mesh into bits

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "IOstreams.H"
#include "SLPtrList.H"
#include "boolList.H"
#include "cellList.H"
#include "primitiveMesh.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void domainDecomposition::decomposeMesh(const bool filterEmptyPatches)
{
    // Decide which cell goes to which processor
    distributeCells();

    // Distribute the cells according to the given processor label

    // calculate the addressing information for the original mesh
    Info<< "\nCalculating original mesh data" << endl;

    // set references to the original mesh
    const polyBoundaryMesh& patches = boundaryMesh();

    // Access all faces to grab the zones
    const faceList& fcs = allFaces();
    const labelList& owner = faceOwner();
    const labelList& neighbour = faceNeighbour();

    // loop through the list of processor labels for the cell and add the
    // cell shape to the list of cells for the appropriate processor

    Info<< "\nDistributing cells to processors" << endl;

    // Memory management
    {
        List<SLList<label> > procCellList(nProcs_);

        forAll (cellToProc_, celli)
        {
            if (cellToProc_[celli] >= nProcs_)
            {
                FatalErrorIn("domainDecomposition::decomposeMesh()")
                    << "Impossible processor label " << cellToProc_[celli]
                    << "for cell " << celli
                    << abort(FatalError);
            }
            else
            {
                procCellList[cellToProc_[celli]].append(celli);
            }
        }

        // Convert linked lists into normal lists
        forAll (procCellList, procI)
        {
            procCellAddressing_[procI] = procCellList[procI];
        }
    }

    Info << "\nDistributing faces to processors" << endl;

    // Loop through internal faces and decide which processor they belong to
    // First visit all internal faces. If cells at both sides belong to the
    // same processor, the face is an internal face. If they are different,
    // it belongs to both processors.

    // Memory management
    {
        List<SLList<label> > procFaceList(nProcs_);

        forAll (neighbour, facei)
        {
            if (cellToProc_[owner[facei]] == cellToProc_[neighbour[facei]])
            {
                // Face internal to processor
                procFaceList[cellToProc_[owner[facei]]].append(facei);
            }
        }

        // Record number of internal faces on each processor
        forAll (procFaceList, procI)
        {
            nInternalProcFaces_[procI] = procFaceList[procI].size();
        }

        // Detect inter-processor boundaries

        // Neighbour processor for each subdomain
        List<SLList<label> > interProcBoundaries(nProcs_);

        // Face labels belonging to each inter-processor boundary
        List<SLList<SLList<label> > > interProcBFaces(nProcs_);

        List<SLList<label> > procPatchIndex(nProcs_);

        forAll (neighbour, facei)
        {
            if (cellToProc_[owner[facei]] != cellToProc_[neighbour[facei]])
            {
                // inter - processor patch face found. Go through the list of
                // inside boundaries for the owner processor and try to find
                // this inter-processor patch.

                label ownerProc = cellToProc_[owner[facei]];
                label neighbourProc = cellToProc_[neighbour[facei]];

                SLList<label>::iterator curInterProcBdrsOwnIter =
                    interProcBoundaries[ownerProc].begin();

                SLList<SLList<label> >::iterator curInterProcBFacesOwnIter =
                    interProcBFaces[ownerProc].begin();

                bool interProcBouFound = false;

                // WARNING: Synchronous SLList iterators

                for
                (
                    ;
                    curInterProcBdrsOwnIter
                 != interProcBoundaries[ownerProc].end()
                 && curInterProcBFacesOwnIter
                 != interProcBFaces[ownerProc].end();
                    ++curInterProcBdrsOwnIter, ++curInterProcBFacesOwnIter
                )
                {
                    if (curInterProcBdrsOwnIter() == neighbourProc)
                    {
                        // the inter - processor boundary exists. Add the face
                        interProcBouFound = true;

                        curInterProcBFacesOwnIter().append(facei);

                        SLList<label>::iterator curInterProcBdrsNeiIter =
                            interProcBoundaries[neighbourProc].begin();

                        SLList<SLList<label> >::iterator
                            curInterProcBFacesNeiIter =
                            interProcBFaces[neighbourProc].begin();

                        bool neighbourFound = false;

                        // WARNING: Synchronous SLList iterators

                        for
                        (
                            ;
                            curInterProcBdrsNeiIter !=
                            interProcBoundaries[neighbourProc].end()
                         && curInterProcBFacesNeiIter !=
                            interProcBFaces[neighbourProc].end();
                            ++curInterProcBdrsNeiIter,
                            ++curInterProcBFacesNeiIter
                        )
                        {
                            if (curInterProcBdrsNeiIter() == ownerProc)
                            {
                                // boundary found. Add the face
                                neighbourFound = true;

                                curInterProcBFacesNeiIter().append(facei);
                            }

                            if (neighbourFound) break;
                        }

                        if (interProcBouFound && !neighbourFound)
                        {
                            FatalErrorIn
                            (
                                "domainDecomposition::decomposeMesh()"
                            )   << "Inconsistency in inter - "
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
                    interProcBFaces[ownerProc].append(SLList<label>(facei));

                    // neighbour
                    interProcBoundaries[neighbourProc].append(ownerProc);
                    interProcBFaces[neighbourProc].append
                    (
                        SLList<label>(facei)
                    );
                }
            }
        }

        // Loop through patches. For cyclic boundaries detect inter-processor
        // faces; for all other, add faces to the face list and remember start
        // and size of all patches.

        // for all processors, set the size of start index and patch size
        // lists to the number of patches in the mesh
        forAll (procPatchSize_, procI)
        {
            procPatchSize_[procI].setSize(patches.size());
            procPatchStartIndex_[procI].setSize(patches.size());
        }

        forAll (patches, patchi)
        {
            // Reset size and start index for all processors
            forAll (procPatchSize_, procI)
            {
                procPatchSize_[procI][patchi] = 0;
                procPatchStartIndex_[procI][patchi] =
                    procFaceList[procI].size();
            }

            const label patchStart = patches[patchi].start();

            if (!isA<cyclicPolyPatch>(patches[patchi]))
            {
                // Normal patch. Add faces to processor where the cell
                // next to the face lives

                const unallocLabelList& patchFaceCells =
                    patches[patchi].faceCells();

                forAll (patchFaceCells, facei)
                {
                    const label curProc = cellToProc_[patchFaceCells[facei]];

                    // add the face
                    procFaceList[curProc].append(patchStart + facei);

                    // increment the number of faces for this patch
                    procPatchSize_[curProc][patchi]++;
                }
            }
            else
            {
                // Cyclic patch special treatment

                const polyPatch& cPatch = patches[patchi];

                const label cycOffset = cPatch.size()/2;

                // Set reference to faceCells for both patches
                const labelList::subList firstFaceCells
                (
                    cPatch.faceCells(),
                    cycOffset
                );

                const labelList::subList secondFaceCells
                (
                    cPatch.faceCells(),
                    cycOffset,
                    cycOffset
                );

                forAll (firstFaceCells, facei)
                {
                    if
                    (
                        cellToProc_[firstFaceCells[facei]]
                     != cellToProc_[secondFaceCells[facei]]
                    )
                    {
                        // This face becomes an inter-processor boundary face
                        // inter - processor patch face found. Go through
                        // the list of inside boundaries for the owner
                        // processor and try to find this inter-processor
                        // patch.

                        cyclicParallel_ = true;

                        label ownerProc = cellToProc_[firstFaceCells[facei]];
                        label neighbourProc =
                            cellToProc_[secondFaceCells[facei]];

                        SLList<label>::iterator curInterProcBdrsOwnIter =
                            interProcBoundaries[ownerProc].begin();

                        SLList<SLList<label> >::iterator
                            curInterProcBFacesOwnIter =
                            interProcBFaces[ownerProc].begin();

                        bool interProcBouFound = false;

                        // WARNING: Synchronous SLList iterators

                        for
                        (
                            ;
                            curInterProcBdrsOwnIter !=
                            interProcBoundaries[ownerProc].end()
                         && curInterProcBFacesOwnIter !=
                            interProcBFaces[ownerProc].end();
                            ++curInterProcBdrsOwnIter,
                            ++curInterProcBFacesOwnIter
                        )
                        {
                            if (curInterProcBdrsOwnIter() == neighbourProc)
                            {
                                // the inter - processor boundary exists.
                                // Add the face
                                interProcBouFound = true;

                                curInterProcBFacesOwnIter().append
                                    (patchStart + facei);

                                SLList<label>::iterator curInterProcBdrsNeiIter
                                   = interProcBoundaries[neighbourProc].begin();

                                SLList<SLList<label> >::iterator
                                    curInterProcBFacesNeiIter =
                                    interProcBFaces[neighbourProc].begin();

                                bool neighbourFound = false;

                                // WARNING: Synchronous SLList iterators

                                for
                                (
                                    ;
                                    curInterProcBdrsNeiIter
                                   != interProcBoundaries[neighbourProc].end()
                                 && curInterProcBFacesNeiIter
                                   != interProcBFaces[neighbourProc].end();
                                    ++curInterProcBdrsNeiIter,
                                    ++curInterProcBFacesNeiIter
                                )
                                {
                                    if (curInterProcBdrsNeiIter() == ownerProc)
                                    {
                                        // boundary found. Add the face
                                        neighbourFound = true;

                                        curInterProcBFacesNeiIter()
                                            .append
                                            (
                                                patchStart
                                              + cycOffset
                                              + facei
                                            );
                                    }

                                    if (neighbourFound) break;
                                }

                                if (interProcBouFound && !neighbourFound)
                                {
                                    FatalErrorIn
                                    (
                                        "domainDecomposition::decomposeMesh()"
                                    )   << "Inconsistency in inter-processor "
                                        << "boundary lists for processors "
                                        << ownerProc << " and "
                                        << neighbourProc
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
                            interProcBFaces[ownerProc]
                                .append(SLList<label>(patchStart + facei));

                            // neighbour
                            interProcBoundaries[neighbourProc]
                                .append(ownerProc);
                            interProcBFaces[neighbourProc]
                                .append
                                (
                                    SLList<label>
                                    (
                                        patchStart
                                      + cycOffset
                                      + facei
                                    )
                                );
                        }
                    }
                    else
                    {
                        // This cyclic face remains on the processor
                        label ownerProc = cellToProc_[firstFaceCells[facei]];

                        // add the face
                        procFaceList[ownerProc].append(patchStart + facei);

                        // increment the number of faces for this patch
                        procPatchSize_[ownerProc][patchi]++;

                        // Note: I cannot add the other side of the cyclic
                        // boundary here because this would violate the order.
                        // They will be added in a separate loop below
                        // HJ, 15/Jan/2001
                    }
                }

                // Ordering in cyclic boundaries is important.
                // Add the other half of cyclic faces for cyclic boundaries
                // that remain on the processor
                forAll (secondFaceCells, facei)
                {
                    if
                    (
                        cellToProc_[firstFaceCells[facei]]
                     == cellToProc_[secondFaceCells[facei]]
                    )
                    {
                        // This cyclic face remains on the processor
                        label ownerProc = cellToProc_[firstFaceCells[facei]];

                        // add the second face
                        procFaceList[ownerProc].append
                            (patchStart + cycOffset + facei);

                        // increment the number of faces for this patch
                        procPatchSize_[ownerProc][patchi]++;
                    }
                }
            }
        }

        // Face zone treatment.  HJ, 27/Mar/2009
        // Face zones identified as global will be present on all CPUs
        List<SLList<label> > procZoneFaceList(nProcs_);

        if (decompositionDict_.found("globalFaceZones"))
        {
            wordList fzNames(decompositionDict_.lookup("globalFaceZones"));

            const faceZoneMesh& fz = faceZones();

            forAll (fzNames, nameI)
            {
                const label zoneID = fz.findZoneID(fzNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn("domainDecomposition::decomposeMesh()")
                        << "Unknown global face zone " << fzNames[nameI]
                        << nl << "Valid face zones are" << fz.names()
                        << exit(FatalError);
                }

                Info<< "Preserving global face zone " << fzNames[nameI]
                    << endl;

                const faceZone& curFz = fz[zoneID];

                // Go through all the faces in the zone.  If the owner of the
                // face equals to current processor, it has already been added;
                // otherwise, add the face to all processor face lists
                forAll (curFz, faceI)
                {
                    const label curFaceID = curFz[faceI];

                    if (curFaceID < owner.size())
                    {
                        // This is an active mesh face, and it already belongs
                        // to one CPU.  Find out which and add it to the others

                        const label curProc = cellToProc_[owner[curFaceID]];

                        forAll (procZoneFaceList, procI)
                        {
                            if (procI != curProc)
                            {
                                procZoneFaceList[procI].append(curFaceID);
                            }
                        }
                    }
                    else
                    {
                        // This is a stand-alone face, add it to all processors
                        forAll (procFaceList, procI)
                        {
                            procZoneFaceList[procI].append(curFaceID);
                        }
                    }
                }
            }
        }


        // Convert linked lists into normal lists
        // Add inter-processor boundaries and remember start indices

        forAll (procFaceList, procI)
        {
            // Get internal and regular boundary processor faces
            const SLList<label>& curProcFaces = procFaceList[procI];

            // Get reference to processor face addressing
            labelList& curProcFaceAddressing = procFaceAddressing_[procI];

            labelList& curProcNeighbourProcessors =
                procNeighbourProcessors_[procI];

            labelList& curProcProcessorPatchSize =
                procProcessorPatchSize_[procI];

            labelList& curProcProcessorPatchStartIndex =
                procProcessorPatchStartIndex_[procI];

            // calculate the size
            label nFacesOnProcessor = curProcFaces.size();

            for
            (
                SLList<SLList<label> >::const_iterator curInterProcBFacesIter =
                    interProcBFaces[procI].begin();
                curInterProcBFacesIter != interProcBFaces[procI].end();
                ++curInterProcBFacesIter
            )
            {
                nFacesOnProcessor += curInterProcBFacesIter().size();
            }

            // Add stand-alone global zone faces
            nFacesOnProcessor += procZoneFaceList[procI].size();

            curProcFaceAddressing.setSize(nFacesOnProcessor);

            // Fill in the list. Calculate turning index.
            // Turning index will be -1 only for some faces on processor
            // boundaries, i.e. the ones where the current processor ID
            // is in the cell which is a face neighbour.
            // Turning index is stored as the sign of the face addressing list

            label nFaces = 0;

            // Add internal and boundary faces
            // Remember to increment the index by one such that the
            // turning index works properly.  HJ, 5/Dec/2001
            for
            (
                SLList<label>::const_iterator curProcFacesIter =
                    curProcFaces.begin();
                curProcFacesIter != curProcFaces.end();
                ++curProcFacesIter
            )
            {
                curProcFaceAddressing[nFaces] = curProcFacesIter() + 1;
                nFaces++;
            }

            // Add inter-processor boundary faces. At the beginning of each
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

            SLList<SLList<label> >::iterator curInterProcBFacesIter =
                interProcBFaces[procI].begin();

            for
            (
                ;
                curInterProcBdrsIter != interProcBoundaries[procI].end()
             && curInterProcBFacesIter != interProcBFaces[procI].end();
                ++curInterProcBdrsIter, ++curInterProcBFacesIter
            )
            {
                curProcNeighbourProcessors[nProcPatches] =
                    curInterProcBdrsIter();

                // Get start index for processor patch
                curProcProcessorPatchStartIndex[nProcPatches] = nFaces;

                label& curSize =
                    curProcProcessorPatchSize[nProcPatches];

                curSize = 0;

                // Add faces for this processor boundary

                for
                (
                    SLList<label>::iterator curFacesIter =
                        curInterProcBFacesIter().begin();
                    curFacesIter != curInterProcBFacesIter().end();
                    ++curFacesIter
                )
                {
                    // Add the face

                    // Remember to increment the index by one such that the
                    // turning index works properly.  HJ, 5/Dec/2001
                    if (cellToProc_[owner[curFacesIter()]] == procI)
                    {
                        curProcFaceAddressing[nFaces] = curFacesIter() + 1;
                    }
                    else
                    {
                        // Turning face
                        curProcFaceAddressing[nFaces] = -(curFacesIter() + 1);
                    }

                    // increment the size
                    curSize++;

                    nFaces++;
                }

                nProcPatches++;
            }

            // Record number of live faces
            nLiveProcFaces_[procI] = nFaces;

            // Add stand-alone face zone faces
            const SLList<label>& curProcZoneFaces = procZoneFaceList[procI];

            for
            (
                SLList<label>::const_iterator curProcZoneFacesIter =
                    curProcZoneFaces.begin();
                curProcZoneFacesIter != curProcZoneFaces.end();
                ++curProcZoneFacesIter
            )
            {
                curProcFaceAddressing[nFaces] = curProcZoneFacesIter() + 1;
                nFaces++;
            }
        } // End for all processors
    } // End of memory management

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

        forAll (oldPatchSizes, patchi)
        {
            if (!filterEmptyPatches || oldPatchSizes[patchi] > 0)
            {
                curBoundaryAddressing[nPatches] = patchi;

                curPatchSizes[nPatches] = oldPatchSizes[patchi];

                curPatchStarts[nPatches] = oldPatchStarts[patchi];

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

    Info << "\nDistributing points to processors" << endl;
    // For every processor, loop through the list of faces for the processor.
    // For every face, loop through the list of points and mark the point as
    // used for the processor. Collect the list of used points for the
    // processor.

    // Record number of live points on each processor
    labelList nLivePoints(nProcs_, 0);

    forAll (procPointAddressing_, procI)
    {
        // Dimension list to all points in the mesh.  HJ, 27/Mar/2009
        boolList pointLabels(allPoints().size(), false);

        // Get reference to list of used faces
        const labelList& procFaceLabels = procFaceAddressing_[procI];

        // Collect the used points
        labelList& procPointLabels = procPointAddressing_[procI];

        procPointLabels.setSize(pointLabels.size());

        // Two-pass algorithm:
        // First loop through live faces and record them in the points list
        // Second, visit all inactive zone faces and record the points

        label nUsedPoints = 0;

        // First pass: live faces

        for (label faceI = 0; faceI < nLiveProcFaces_[procI]; faceI++)
        {
            // Because of the turning index, some labels may be negative
            const labelList& facePoints = fcs[mag(procFaceLabels[faceI]) - 1];

            forAll (facePoints, pointI)
            {
                // Mark the point as used
                pointLabels[facePoints[pointI]] = true;
            }
        }

        forAll (pointLabels, pointI)
        {
            if (pointLabels[pointI])
            {
                procPointLabels[nUsedPoints] = pointI;

                nUsedPoints++;
            }
        }

        // Record number of live points
        nLivePoints[procI] = nUsedPoints;

        // Second pass: zone faces

        // Reset point usage list
        boolList pointLabelsSecondPass(allPoints().size(), false);

        for
        (
            label faceI = nLiveProcFaces_[procI];
            faceI < procFaceLabels.size();
            faceI++
        )
        {
            // Because of the turning index, some labels may be negative
            const labelList& facePoints = fcs[mag(procFaceLabels[faceI]) - 1];

            forAll (facePoints, pointI)
            {
                // Mark the point as used
                if (!pointLabels[facePoints[pointI]])
                {
                    pointLabelsSecondPass[facePoints[pointI]] = true;
                }
            }
        }

        forAll (pointLabelsSecondPass, pointI)
        {
            if (pointLabelsSecondPass[pointI])
            {
                procPointLabels[nUsedPoints] = pointI;

                nUsedPoints++;
            }
        }

        // Reset the size of used points
        procPointLabels.setSize(nUsedPoints);
    }

    // Gather data about globally shared points

    // Memory management
    {
        // Dimension list to all points in the mesh.  HJ, 27/Mar/2009
        labelList pointsUsage(allPoints().size(), 0);

        // Globally shared points are the ones used by more than 2 processors
        // Size the list approximately and gather the points
        labelHashSet gSharedPoints
        (
            min(100, nPoints()/1000)
        );

        // Loop through all the processors and mark up points used by
        // processor boundaries.  When a point is used twice, it is a
        // globally shared point

        for (label procI = 0; procI < nProcs_; procI++)
        {
            // Get list of face labels
            const labelList& curFaceLabels = procFaceAddressing_[procI];

            // Get start of processor faces
            const labelList& curProcessorPatchStarts =
                procProcessorPatchStartIndex_[procI];

            const labelList& curProcessorPatchSizes =
                procProcessorPatchSize_[procI];

            // Reset the lookup list
            pointsUsage = 0;

            forAll (curProcessorPatchStarts, patchi)
            {
                const label curStart = curProcessorPatchStarts[patchi];
                const label curEnd = curStart + curProcessorPatchSizes[patchi];

                for
                (
                    label faceI = curStart;
                    faceI < curEnd;
                    faceI++
                )
                {
                    // Mark the original face as used
                    // Remember to decrement the index by one (turning index)
                    // HJ, 5/Dec/2001
                    const label curF = mag(curFaceLabels[faceI]) - 1;

                    const face& f = fcs[curF];

                    forAll (f, pointI)
                    {
                        if (pointsUsage[f[pointI]] == 0)
                        {
                            // Point not previously used
                            pointsUsage[f[pointI]] = patchi + 1;
                        }
                        else if (pointsUsage[f[pointI]] != patchi + 1)
                        {
                            // Point used by some other patch = global point!
                            gSharedPoints.insert(f[pointI]);
                        }
                    }
                }
            }
        }

        // Grab the result from the hash list
        globallySharedPoints_ = gSharedPoints.toc();
        sort(globallySharedPoints_);
    }
}
