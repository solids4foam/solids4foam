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

\*---------------------------------------------------------------------------*/

#include "advancingFrontAMIS4F.H"
#include "mergePoints.H"
#include "mapDistribute.H"
#include "AABBTree.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::advancingFrontAMIS4F::calcOverlappingProcs
(
    const List<treeBoundBoxList>& procBb,
    const treeBoundBox& bb,
    boolList& overlaps
) const
{
    overlaps.setSize(procBb.size());
    overlaps = false;

    label nOverlaps = 0;

    forAll(procBb, proci)
    {
        const treeBoundBoxList& bbp = procBb[proci];

        for (const treeBoundBox& tbb: bbp)
        {
            if (tbb.overlaps(bb))
            {
                overlaps[proci] = true;
                ++nOverlaps;
                break;
            }
        }
    }

    return nOverlaps;
}


void Foam::advancingFrontAMIS4F::distributePatches
(
    const mapDistribute& map,
    const standAlonePatch& pp,
    const globalIndex& gi,
    List<faceList>& faces,
    List<pointField>& points,
    List<labelList>& faceIDs
) const
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (const int domain : Pstream::allProcs())
    {
        const labelList& sendElems = map.subMap()[domain];

        if (domain != Pstream::myProcNo() && sendElems.size())
        {
            faceList subFaces(UIndirectList<face>(pp, sendElems));
            standAlonePatch subPatch
            (
                SubList<face>(subFaces, subFaces.size()),
                pp.points()
            );

            if (debug & 2)
            {
                Pout<< "distributePatches: to processor " << domain
                    << " sending faces " << subPatch.faceCentres() << endl;
            }

            UOPstream toDomain(domain, pBufs);
            toDomain
                << subPatch.localFaces() << subPatch.localPoints()
                << gi.toGlobal(sendElems);
        }
    }

    // Start receiving
    pBufs.finishedSends();

    faces.setSize(Pstream::nProcs());
    points.setSize(Pstream::nProcs());
    faceIDs.setSize(Pstream::nProcs());

    {
        // Set up 'send' to myself
        const labelList& sendElems = map.subMap()[Pstream::myProcNo()];
        faceList subFaces(UIndirectList<face>(pp, sendElems));
        standAlonePatch subPatch
        (
            SubList<face>(subFaces, subFaces.size()),
            pp.points()
        );

        // Receive
        if (debug & 2)
        {
            Pout<< "distributePatches: to processor " << Pstream::myProcNo()
                << " sending faces " << subPatch.faceCentres() << endl;
        }

        faces[Pstream::myProcNo()] = subPatch.localFaces();
        points[Pstream::myProcNo()] = subPatch.localPoints();
        faceIDs[Pstream::myProcNo()] = gi.toGlobal(sendElems);
    }

    // Consume
    for (const int domain : Pstream::allProcs())
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);

            str >> faces[domain]
                >> points[domain]
                >> faceIDs[domain];
        }
    }
}


void Foam::advancingFrontAMIS4F::distributeAndMergePatches
(
    const mapDistribute& map,
    const standAlonePatch& tgtPatch,
    const globalIndex& gi,
    faceList& tgtFaces,
    pointField& tgtPoints,
    labelList& tgtFaceIDs
) const
{
    // Exchange per-processor data
    List<faceList> allFaces;
    List<pointField> allPoints;
    List<labelList> allTgtFaceIDs;
    distributePatches(map, tgtPatch, gi, allFaces, allPoints, allTgtFaceIDs);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, proci)
    {
        nFaces += allFaces[proci].size();
        nPoints += allPoints[proci].size();
    }

    tgtFaces.setSize(nFaces);
    tgtPoints.setSize(nPoints);
    tgtFaceIDs.setSize(nFaces);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const labelList& faceIDs = allTgtFaceIDs[Pstream::myProcNo()];
        SubList<label>(tgtFaceIDs, faceIDs.size()) = faceIDs;

        const faceList& fcs = allFaces[Pstream::myProcNo()];
        for (const face& f : fcs)
        {
            face& newF = tgtFaces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        for (const point& pt: pts)
        {
            tgtPoints[nPoints++] = pt;
        }
    }


    // Other proc data follows
    forAll(allFaces, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const labelList& faceIDs = allTgtFaceIDs[proci];
            SubList<label>(tgtFaceIDs, faceIDs.size(), nFaces) = faceIDs;

            const faceList& fcs = allFaces[proci];
            for (const face& f : fcs)
            {
                face& newF = tgtFaces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[proci];
            for (const point& pt: pts)
            {
                tgtPoints[nPoints++] = pt;
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newTgtPoints;
    bool hasMerged = mergePoints
    (
        tgtPoints,
        SMALL,
        false,
        oldToNew,
        newTgtPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << tgtPoints.size()
                << " down to " << newTgtPoints.size() << " points" << endl;
        }

        tgtPoints.transfer(newTgtPoints);
        for (face& f : tgtFaces)
        {
            inplaceRenumber(oldToNew, f);
        }
    }
}


Foam::autoPtr<Foam::mapDistribute> Foam::advancingFrontAMIS4F::calcProcMap
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch
) const
{
    // Get decomposition of patch
    List<treeBoundBoxList> procBb(Pstream::nProcs());

    if (srcPatch.size())
    {
        procBb[Pstream::myProcNo()] =
            AABBTree<face>
            (
                srcPatch.localFaces(),
                srcPatch.localPoints(),
                false
            ).boundBoxes();
    }
    else
    {
        procBb[Pstream::myProcNo()] = treeBoundBoxList();
    }

    Pstream::gatherList(procBb);
    Pstream::scatterList(procBb);

    if (debug)
    {
        Info<< "Determining extent of srcPatch per processor:" << nl
            << "\tproc\tbb" << endl;
        forAll(procBb, proci)
        {
            Info<< '\t' << proci << '\t' << procBb[proci] << endl;
        }
    }

    // Determine which faces of tgtPatch overlaps srcPatch per proc
    const faceList& faces = tgtPatch.localFaces();
    const pointField& points = tgtPatch.localPoints();

    labelListList sendMap;

    {
        // Per processor indices into all segments to send
        List<DynamicList<label>> dynSendMap(Pstream::nProcs());

        // Work array - whether processor bb overlaps the face bounds
        boolList procBbOverlaps(Pstream::nProcs());

        forAll(faces, facei)
        {
            if (faces[facei].size())
            {
                treeBoundBox faceBb(points, faces[facei]);

                // Find the processor this face overlaps
                calcOverlappingProcs(procBb, faceBb, procBbOverlaps);

                forAll(procBbOverlaps, proci)
                {
                    if (procBbOverlaps[proci])
                    {
                        dynSendMap[proci].append(facei);
                    }
                }
            }
        }

        // Convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
    }

    // Debug printing
    if (debug)
    {
        Pout<< "Of my " << faces.size() << " I need to send to:" << nl
            << "\tproc\tfaces" << endl;
        forAll(sendMap, proci)
        {
            Pout<< '\t' << proci << '\t' << sendMap[proci].size() << endl;
        }
    }


    // Send over how many faces I need to receive
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        sendSizes[Pstream::myProcNo()][proci] = sendMap[proci].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // My local segment first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me
            label nRecv = sendSizes[proci][Pstream::myProcNo()];
            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; ++i)
            {
                constructMap[proci][i] = segmentI++;
            }
        }
    }

    return autoPtr<mapDistribute>::New
    (
        segmentI, // size after construction
        std::move(sendMap),
        std::move(constructMap)
    );
}


// ************************************************************************* //
