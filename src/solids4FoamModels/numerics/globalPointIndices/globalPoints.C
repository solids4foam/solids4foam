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

\*---------------------------------------------------------------------------*/

#include "globalPoints.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::globalPoints, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Total number of points on processor patches. Is upper limit for number
// of shared points
Foam::label Foam::globalPoints::countPatchPoints
(
    const polyBoundaryMesh& patches
)
{
    label nTotPoints = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            (Pstream::parRun() && isA<processorPolyPatch>(pp))
         || isA<cyclicPolyPatch>(pp)
        )
        {
            nTotPoints += pp.nPoints();
        }
    }

    return nTotPoints;
}


// Collect all topological information about a point on a patch.
// (this information is the patch faces using the point and the relative
// position of the point in the face)
void Foam::globalPoints::addToSend
(
    const primitivePatch& pp,
    const label patchPointI,
    const procPointList& knownInfo,

    dynamicLabelList& patchFaces,
    dynamicLabelList& indexInFace,
    DynamicList<procPointList>& allInfo
)
{
    label meshPointI = pp.meshPoints()[patchPointI];

    // Add all faces using the point so we are sure we find it on the
    // other side.
    const labelList& pFaces = pp.pointFaces()[patchPointI];

    forAll(pFaces, i)
    {
        label patchFaceI = pFaces[i];

        const face& f = pp[patchFaceI];

        patchFaces.append(patchFaceI);
        indexInFace.append(findIndex(f, meshPointI));
        allInfo.append(knownInfo);
    }
}


// Add nbrInfo to myInfo. Return true if anything changed.
// nbrInfo is for a point a list of all the processors using it (and the
// meshPoint label on that processor)
bool Foam::globalPoints::mergeInfo
(
    const procPointList& nbrInfo,
    procPointList& myInfo
)
{
    // Indices of entries in nbrInfo not yet in myInfo.
    dynamicLabelList newInfo(nbrInfo.size());

    forAll(nbrInfo, i)
    {
        const procPoint& info = nbrInfo[i];

        // Check if info already in myInfo.
        label index = -1;

        forAll(myInfo, j)
        {
            if (myInfo[j] == info)
            {
                // Already have information for processor/point combination
                // in my list so skip.
                index = j;

                break;
            }
        }

        if (index == -1)
        {
            // Mark this information as being not yet in myInfo
            newInfo.append(i);
        }
    }

    newInfo.shrink();

    // Append all nbrInfos referenced in newInfo to myInfo.

    label index = myInfo.size();

    myInfo.setSize(index + newInfo.size());

    forAll(newInfo, i)
    {
        myInfo[index++] = nbrInfo[newInfo[i]];
    }

    // Did anything change?
    return newInfo.size() > 0;
}


// Updates database of current information on meshpoints with nbrInfo.
// Uses mergeInfo above. Returns true if data kept for meshPointI changed.
bool Foam::globalPoints::storeInfo
(
    const procPointList& nbrInfo,
    const label meshPointI
)
{
    label infoChanged = false;

    // Get the index into the procPoints list.
    Map<label>::iterator iter = meshToProcPoint_.find(meshPointI);

    if (iter != meshToProcPoint_.end())
    {
        procPointList& knownInfo = procPoints_[iter()];

        if (mergeInfo(nbrInfo, knownInfo))
        {
            infoChanged = true;
        }
    }
    else
    {
        procPointList knownInfo(1);
        knownInfo[0][0] = Pstream::myProcNo();
        knownInfo[0][1] = meshPointI;

        if (mergeInfo(nbrInfo, knownInfo))
        {
            // Update addressing from into procPoints
            meshToProcPoint_.insert(meshPointI, procPoints_.size());
            // Insert into list of equivalences.
            procPoints_.append(knownInfo);

            infoChanged = true;
        }
    }
    return infoChanged;
}


// Insert my own points into structure and mark as changed.
void Foam::globalPoints::initOwnPoints
(
    const bool allPoints,
    labelHashSet& changedPoints
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            (Pstream::parRun() && isA<processorPolyPatch>(pp))
         || isA<cyclicPolyPatch>(pp)
        )
        {
            const labelList& meshPoints = pp.meshPoints();

            if (allPoints)
            {
                // All points on patch
                forAll(meshPoints, i)
                {
                    label meshPointI = meshPoints[i];

                    procPointList knownInfo(1);
                    knownInfo[0][0] = Pstream::myProcNo();
                    knownInfo[0][1] = meshPointI;

                    // Update addressing from meshpoint to index in procPoints
                    meshToProcPoint_.insert(meshPointI, procPoints_.size());
                    // Insert into list of equivalences.
                    procPoints_.append(knownInfo);

                    // Update changedpoints info.
                    changedPoints.insert(meshPointI);
                }
            }
            else
            {
                // Boundary points only
                const labelList& boundaryPoints = pp.boundaryPoints();

                forAll(boundaryPoints, i)
                {
                    label meshPointI = meshPoints[boundaryPoints[i]];

                    procPointList knownInfo(1);
                    knownInfo[0][0] = Pstream::myProcNo();
                    knownInfo[0][1] = meshPointI;

                    // Update addressing from meshpoint to index in procPoints
                    meshToProcPoint_.insert(meshPointI, procPoints_.size());
                    // Insert into list of equivalences.
                    procPoints_.append(knownInfo);

                    // Update changedpoints info.
                    changedPoints.insert(meshPointI);
                }
            }
        }
    }
}


// Send all my info on changedPoints_ to my neighbours.
void Foam::globalPoints::sendPatchPoints(const labelHashSet& changedPoints)
 const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            // Information to send:
            // patch face
            dynamicLabelList patchFaces(pp.nPoints());
            // index in patch face
            dynamicLabelList indexInFace(pp.nPoints());
            // all information I currently hold about this patchPoint
            DynamicList<procPointList> allInfo(pp.nPoints());


            // Now collect information on all mesh points mentioned in
            // changedPoints. Note that these points only should occur on
            // processorPatches (or rather this is a limitation!).

            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, patchPointI)
            {
                label meshPointI = meshPoints[patchPointI];

                if (changedPoints.found(meshPointI))
                {
                    label index = meshToProcPoint_[meshPointI];

                    const procPointList& knownInfo = procPoints_[index];

                    // Add my information about meshPointI to the send buffers
                    addToSend
                    (
                        pp,
                        patchPointI,
                        knownInfo,

                        patchFaces,
                        indexInFace,
                        allInfo
                    );
                }
            }
            patchFaces.shrink();
            indexInFace.shrink();
            allInfo.shrink();

            // Send to neighbour
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                if (debug)
                {
                    Pout<< " Sending to "
                        << procPatch.neighbProcNo() << "   point information:"
                        << patchFaces.size() << endl;
                }

                OPstream toNeighbour
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );

                toNeighbour << patchFaces << indexInFace << allInfo;
            }
        }
    }
}


// Receive all my neighbours' information and merge with mine.
// After finishing will have updated
// - procPoints_ : all neighbour information merged in.
// - meshToProcPoint_
// - changedPoints: all meshPoints for which something changed.
void Foam::globalPoints::receivePatchPoints(labelHashSet& changedPoints)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Reset changed points
    changedPoints.clear();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            labelList patchFaces;
            labelList indexInFace;
            List<procPointList> nbrInfo;

            {
                IPstream fromNeighbour
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );
                fromNeighbour >> patchFaces >> indexInFace >> nbrInfo;
            }

            if (debug)
            {
                Pout<< " Received from "
                    << procPatch.neighbProcNo() << "   point information:"
                    << patchFaces.size() << endl;
            }

            forAll(patchFaces, i)
            {
                const face& f = pp[patchFaces[i]];

                // Get index in this face from index on face on other side.
                label index = (f.size() - indexInFace[i]) % f.size();

                // Get the meshpoint on my side
                label meshPointI = f[index];

                //const procPointList& info = nbrInfo[i];
                //Pout << "Received for my coord "
                //    << mesh_.points()[meshPointI];
                //
                //forAll(info, j)
                //{
                //    Pout<< ' ' <<info[j];
                //}
                //Pout<< endl;

                if (storeInfo(nbrInfo[i], meshPointI))
                {
                    changedPoints.insert(meshPointI);
                }
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            // Handle cyclics: send lower half to upper half and vice versa.
            // Or since they both are in memory just do it point by point.

            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            const labelList& meshPoints = pp.meshPoints();

            //const edgeList& connections = cycPatch.coupledPoints();
            const edgeList connections(coupledPoints(cycPatch));

            forAll(connections, i)
            {
                const edge& e = connections[i];

                label meshPointA = meshPoints[e[0]];
                label meshPointB = meshPoints[e[1]];

                // Do we have information on meshPointA?
                Map<label>::iterator procPointA =
                    meshToProcPoint_.find(meshPointA);

                if (procPointA != meshToProcPoint_.end())
                {
                    // Store A info onto pointB
                    if (storeInfo(procPoints_[procPointA()], meshPointB))
                    {
                        changedPoints.insert(meshPointB);
                    }
                }

                // Same for info on pointB
                Map<label>::iterator procPointB =
                    meshToProcPoint_.find(meshPointB);

                if (procPointB != meshToProcPoint_.end())
                {
                    // Store B info onto pointA
                    if (storeInfo(procPoints_[procPointB()], meshPointA))
                    {
                        changedPoints.insert(meshPointA);
                    }
                }
            }
        }
    }
}


// Remove entries which are handled by normal face-face communication. I.e.
// those points where the equivalence list is only me and my (face)neighbour
void Foam::globalPoints::remove(const Map<label>& directNeighbours)
{
    // Save old ones.
    Map<label> oldMeshToProcPoint(meshToProcPoint_);
    meshToProcPoint_.clear();

    List<procPointList> oldProcPoints;
    oldProcPoints.transfer(procPoints_);

    // Go through all equivalences
    for
    (
        Map<label>::const_iterator iter = oldMeshToProcPoint.begin();
        iter != oldMeshToProcPoint.end();
        ++iter
    )
    {
        label meshPointI = iter.key();
        const procPointList& pointInfo = oldProcPoints[iter()];

        if (pointInfo.size() == 2)
        {
            // I will be in this equivalence list.
            // Check whether my direct (=face) neighbour
            // is in it. This would be an ordinary connection and can be
            // handled by normal face-face connectivity.

            const procPoint& a = pointInfo[0];
            const procPoint& b = pointInfo[1];

            if
            (
                (a[0] == Pstream::myProcNo() && directNeighbours.found(a[1]))
             || (b[0] == Pstream::myProcNo() && directNeighbours.found(b[1]))
            )
            {
                // Normal faceNeighbours
                if (a[0] == Pstream::myProcNo())
                {
                    //Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()[a[1]]
                    //    << endl;
                }
                else if (b[0] == Pstream::myProcNo())
                {
                    //Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()[b[1]]
                    //    << endl;
                }
            }
            else
            {
                // This condition will be very rare: points are used by
                // two processors which are not face-face connected.
                // e.g.
                // +------+------+
                // | wall |  B   |
                // +------+------+
                // |   A  | wall |
                // +------+------+
                // Processor A and B share a point. Note that this only will
                // be found if the two domains are face connected at all
                // (not shown in the picture)

                meshToProcPoint_.insert(meshPointI, procPoints_.size());
                procPoints_.append(pointInfo);
            }
        }
        else if (pointInfo.size() == 1)
        {
            // This happens for 'wedge' like cyclics where the two halves
            // come together in the same point so share the same meshPoint.
            // So this meshPoint will have info of size one only.
            if
            (
                pointInfo[0][0] != Pstream::myProcNo()
             || !directNeighbours.found(pointInfo[0][1])
            )
            {
                meshToProcPoint_.insert(meshPointI, procPoints_.size());
                procPoints_.append(pointInfo);
            }
        }
        else
        {
            meshToProcPoint_.insert(meshPointI, procPoints_.size());
            procPoints_.append(pointInfo);
        }
    }

    procPoints_.shrink();
}


// Get (indices of) points for which I am master (= lowest numbered point on
// lowest numbered processor).
// (equivalence lists should be complete by now)
Foam::labelList Foam::globalPoints::getMasterPoints() const
{
    labelList masterPoints(nPatchPoints_);
    label nMaster = 0;

    // Go through all equivalences and determine meshPoints where I am master.
    for
    (
        Map<label>::const_iterator iter = meshToProcPoint_.begin();
        iter != meshToProcPoint_.end();
        ++iter
    )
    {
        label meshPointI = iter.key();
        const procPointList& pointInfo = procPoints_[iter()];

        if (pointInfo.size() < 2)
        {
            // Points should have an equivalence list >= 2 since otherwise
            // they would be direct neighbours and have been removed in the
            // call to 'remove'.
            Pout<< "MeshPoint:" << meshPointI
                << " coord:" << mesh_.points()[meshPointI]
                << " has no corresponding point on a neighbouring processor"
                << endl;
            FatalErrorIn("globalPoints::getMasterPoints()")
                << '[' << Pstream::myProcNo() << ']'
                << " MeshPoint:" << meshPointI
                << " coord:" << mesh_.points()[meshPointI]
                << " has no corresponding point on a neighbouring processor"
                << abort(FatalError);
        }
        else
        {
            // Find lowest processor and lowest mesh point on this processor.
            label lowestProcI = labelMax;
            label lowestPointI = labelMax;

            forAll(pointInfo, i)
            {
                label proc = pointInfo[i][0];

                if (proc < lowestProcI)
                {
                    lowestProcI = proc;
                    lowestPointI = pointInfo[i][1];
                }
                else if (proc == lowestProcI)
                {
                    lowestPointI = min(lowestPointI, pointInfo[i][1]);
                }
            }

            if
            (
                lowestProcI == Pstream::myProcNo()
             && lowestPointI == meshPointI
            )
            {
                // I am lowest numbered processor and point. Add to my list.
                masterPoints[nMaster++] = meshPointI;
            }
        }
    }

    masterPoints.setSize(nMaster);

    return masterPoints;
}


// Send subset of lists
void Foam::globalPoints::sendSharedPoints(const labelList& changedIndices) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            OPstream toNeighbour(Pstream::blocking, procPatch.neighbProcNo());

            if (debug)
            {
                Pout<< "Sending to " << procPatch.neighbProcNo()
                    << "  changed sharedPoints info:"
                    << changedIndices.size() << endl;
            }

            toNeighbour
                << UIndirectList<label>(sharedPointAddr_, changedIndices)()
                << UIndirectList<label>(sharedPointLabels_, changedIndices)();
        }
    }
}


// Receive shared point indices for all my shared points. Note that since
// there are only a few here we can build a reverse map using the meshpoint
// instead of doing all this relative point indexing (patch face + index in
// face) as in send/receivePatchPoints
void Foam::globalPoints::receiveSharedPoints(labelList& changedIndices)
{
    changedIndices.setSize(sharedPointAddr_.size());
    label nChanged = 0;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Receive and set shared points
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Map from neighbouring meshPoint to sharedPoint)
            Map<label> nbrSharedPoints(sharedPointAddr_.size());

            {
                // Receive meshPoints on neighbour and sharedPoints and build
                // map from it. Note that we could have built the map on the
                // neighbour and sent it over.
                labelList nbrSharedPointAddr;
                labelList nbrSharedPointLabels;

                {
                    IPstream fromNeighbour
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNeighbour >> nbrSharedPointAddr >> nbrSharedPointLabels;
                }

                // Insert into to map
                forAll(nbrSharedPointLabels, i)
                {
                    nbrSharedPoints.insert
                    (
                        nbrSharedPointLabels[i], // meshpoint on neighbour
                        nbrSharedPointAddr[i]    // sharedPoint label
                    );
                }
            }


            // Merge into whatever information I hold.
            for
            (
                Map<label>::const_iterator iter = meshToProcPoint_.begin();
                iter != meshToProcPoint_.end();
                ++iter
            )
            {
                label meshPointI = iter.key();
                label index = iter();

                if (sharedPointAddr_[index] == -1)
                {
                    // No shared point known yet for this meshPoint.
                    // See if was received from neighbour.
                    const procPointList& knownInfo = procPoints_[index];

                    // Check through the whole equivalence list for any
                    // point from the neighbour.
                    forAll(knownInfo, j)
                    {
                        const procPoint& info = knownInfo[j];

                        if
                        (
                            (info[0] == procPatch.neighbProcNo())
                         && nbrSharedPoints.found(info[1])
                        )
                        {
                            // So this knownInfo contains the neighbour point
                            label sharedPointI = nbrSharedPoints[info[1]];

                            sharedPointAddr_[index] = sharedPointI;
                            sharedPointLabels_[index] = meshPointI;
                            changedIndices[nChanged++] = index;

                            break;
                        }
                    }
                }
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            // Build map from meshPoint to sharedPoint
            Map<label> meshToSharedPoint(sharedPointAddr_.size());
            forAll(sharedPointLabels_, i)
            {
                label meshPointI = sharedPointLabels_[i];

                meshToSharedPoint.insert(meshPointI, sharedPointAddr_[i]);
            }

            // Sync all info.
            //const edgeList& connections = cycPatch.coupledPoints();
            const edgeList connections(coupledPoints(cycPatch));

            forAll(connections, i)
            {
                const edge& e = connections[i];

                label meshPointA = pp.meshPoints()[e[0]];
                label meshPointB = pp.meshPoints()[e[1]];

                // Do we already have shared point for meshPointA?
                Map<label>::iterator fndA = meshToSharedPoint.find(meshPointA);
                Map<label>::iterator fndB = meshToSharedPoint.find(meshPointB);

                if (fndA != meshToSharedPoint.end())
                {
                    if (fndB != meshToSharedPoint.end())
                    {
                        if (fndA() != fndB())
                        {
                            FatalErrorIn
                            (
                                "globalPoints::receiveSharedPoints"
                                "(labelList&)"
                            )   << "On patch " << pp.name()
                                << " connected points " << meshPointA
                                << ' ' << mesh_.points()[meshPointA]
                                << " and " << meshPointB
                                << ' ' << mesh_.points()[meshPointB]
                                << " are mapped to different shared points: "
                                << fndA() << " and " << fndB()
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // No shared point yet for B.
                        label sharedPointI = fndA();

                        // Store shared point for meshPointB
                        label index = meshToProcPoint_[meshPointB];

                        sharedPointAddr_[index] = sharedPointI;
                        sharedPointLabels_[index] = meshPointB;
                        changedIndices[nChanged++] = index;
                    }
                }
                else
                {
                    // No shared point yet for A.
                    if (fndB != meshToSharedPoint.end())
                    {
                        label sharedPointI = fndB();

                        // Store shared point for meshPointA
                        label index = meshToProcPoint_[meshPointA];

                        sharedPointAddr_[index] = sharedPointI;
                        sharedPointLabels_[index] = meshPointA;
                        changedIndices[nChanged++] = index;
                    }
                }
            }
        }
    }

    changedIndices.setSize(nChanged);
}


Foam::edgeList Foam::globalPoints::coupledPoints(const cyclicPolyPatch& pp)
{
    // Look at cyclic patch as two halves, A and B.
    // Now all we know is that relative face index in halfA is same
    // as coupled face in halfB and also that the 0th vertex
    // corresponds.

    // From halfA point to halfB or -1.
    labelList coupledPoint(pp.nPoints(), -1);

    for (label patchFaceA = 0; patchFaceA < pp.size()/2; patchFaceA++)
    {
        const face& fA = pp.localFaces()[patchFaceA];

        forAll(fA, indexA)
        {
            label patchPointA = fA[indexA];

            if (coupledPoint[patchPointA] == -1)
            {
                const face& fB = pp.localFaces()[patchFaceA + pp.size()/2];

                label indexB = (fB.size() - indexA) % fB.size();

                coupledPoint[patchPointA] = fB[indexB];
            }
        }
    }

    edgeList connected(pp.nPoints());

    // Extract coupled points.
    label connectedI = 0;

    forAll(coupledPoint, i)
    {
        if (coupledPoint[i] != -1)
        {
            connected[connectedI++] = edge(i, coupledPoint[i]);
        }
    }

    connected.setSize(connectedI);

    return connected;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::globalPoints::globalPoints(const polyMesh& mesh)
:
    mesh_(mesh),
    nPatchPoints_(countPatchPoints(mesh.boundaryMesh())),
    procPoints_(nPatchPoints_),
    meshToProcPoint_(nPatchPoints_),
    sharedPointAddr_(0),
    sharedPointLabels_(0),
    nGlobalPoints_(0)
{
    if (debug)
    {
        Pout<< "globalPoints::globalPoints(const polyMesh&) : "
            << "doing processor to processor communication to get sharedPoints"
            << endl;
    }

    labelHashSet changedPoints(nPatchPoints_);

    // Initialize procPoints with my patch points. Keep track of points
    // inserted (in changedPoints)
    // There are two possible forms of this:
    // - initialize with all patch points (allPoints = true). This causes all
    //   patch points to be exchanged so a lot of information gets stored and
    //   transferred. This all gets filtered out later when removing the
    //   equivalence lists of size 2.
    // - initialize with boundary points of patches only (allPoints = false).
    //   This should work for all decompositions except extreme ones where a
    //   shared point is not on the boundary of any processor patches using it.
    //   This would happen if a domain was pinched such that two patches share
    //   a point or edge.
    initOwnPoints(true, changedPoints);

    // Do one exchange iteration to get neighbour points.
    sendPatchPoints(changedPoints);
    receivePatchPoints(changedPoints);


    // Save neighbours reachable through face-face communication.
    Map<label> neighbourList(meshToProcPoint_);


    // Exchange until nothing changes on all processors.
    bool changed = false;

    do
    {
        sendPatchPoints(changedPoints);
        receivePatchPoints(changedPoints);

        changed = changedPoints.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


    // Remove direct neighbours from point equivalences.
    remove(neighbourList);


    //Pout<< "After removing locally connected points:" << endl;
    //for
    //(
    //    Map<label>::const_iterator iter = meshToProcPoint_.begin();
    //    iter != meshToProcPoint_.end();
    //    ++iter
    //)
    //{
    //    label meshPointI = iter.key();
    //    const procPointList& pointInfo = procPoints_[iter()];
    //
    //    forAll(pointInfo, i)
    //    {
    //        Pout<< "    pointI:" << meshPointI << ' '
    //            << mesh.points()[meshPointI]
    //            << " connected to proc " << pointInfo[i][0]
    //            << " point:" << pointInfo[i][1]
    //        << endl;
    //    }
    //}


    // We now have - in procPoints_ - a list of points which are shared between
    // multiple processors. These are the ones for which are sharedPoint
    // needs to be determined. This is done by having the lowest numbered
    // processor in the equivalence list 'ask' for a sharedPoint number
    // and then distribute it across processor patches to the non-master
    // processors. Note: below piece of coding is not very efficient. Uses
    // a Map where possibly it shouldn't

    // Initialize sharedPoint addressing. Is for every entry in procPoints_
    // the sharedPoint.
    sharedPointAddr_.setSize(meshToProcPoint_.size());
    sharedPointAddr_ = -1;
    sharedPointLabels_.setSize(meshToProcPoint_.size());
    sharedPointLabels_ = -1;


    // Get points for which I am master (lowest numbered proc)
    labelList masterPoints(getMasterPoints());


    // Determine number of master points on all processors.
    labelList sharedPointSizes(Pstream::nProcs());
    sharedPointSizes[Pstream::myProcNo()] = masterPoints.size();

    Pstream::gatherList(sharedPointSizes);
    Pstream::scatterList(sharedPointSizes);

    if (debug)
    {
        Pout<< "sharedPointSizes:" << sharedPointSizes << endl;
    }

    // Get total number of shared points
    nGlobalPoints_ = 0;
    forAll(sharedPointSizes, procI)
    {
        nGlobalPoints_ += sharedPointSizes[procI];
    }

    // Assign sharedPoint labels. Every processor gets assigned consecutive
    // numbers for its master points.
    // These are assigned in processor order so processor0 gets
    // 0..masterPoints.size()-1 etc.

    // My point labels start after those of lower numbered processors
    label sharedPointI = 0;
    for (label procI = 0; procI < Pstream::myProcNo(); procI++)
    {
        sharedPointI += sharedPointSizes[procI];
    }

    forAll(masterPoints, i)
    {
        label meshPointI = masterPoints[i];

        label index = meshToProcPoint_[meshPointI];

        sharedPointLabels_[index] = meshPointI;
        sharedPointAddr_[index] = sharedPointI++;
    }


    // Now we have a sharedPointLabel for some of the entries in procPoints.
    // Send this information to neighbours. Receive their information.
    // Loop until nothing changes.

    // Initial subset to send is points for which I have sharedPoints
    labelList changedIndices(sharedPointAddr_.size());
    label nChanged = 0;

    forAll(sharedPointAddr_, i)
    {
        if (sharedPointAddr_[i] != -1)
        {
            changedIndices[nChanged++] = i;
        }
    }
    changedIndices.setSize(nChanged);

    changed = false;

    do
    {
        if (debug)
        {
            Pout<< "Determined " << changedIndices.size() << " shared points."
                << " Exchanging them" << endl;
        }
        sendSharedPoints(changedIndices);
        receiveSharedPoints(changedIndices);

        changed = changedIndices.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


    forAll(sharedPointLabels_, i)
    {
        if (sharedPointLabels_[i] == -1)
        {
            FatalErrorIn("globalPoints::globalPoints(const polyMesh& mesh)")
                << "Problem: shared point on processor " << Pstream::myProcNo()
                << " not set at index " << sharedPointLabels_[i] << endl
                << "This might mean the individual processor domains are not"
                << " connected and the overall domain consists of multiple"
                << " regions. You can check this with checkMesh"
                << abort(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "globalPoints::globalPoints(const polyMesh&) : "
            << "Finished global points" << endl;
    }
}


// ************************************************************************* //
