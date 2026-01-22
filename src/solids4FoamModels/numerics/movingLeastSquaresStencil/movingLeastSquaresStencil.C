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

#include "movingLeastSquaresStencil.H"
#include "fvMesh.H"
#include "PstreamBuffers.H"
#include "processorPolyPatch.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(movingLeastSquaresStencil, 0);


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

List<scalar> movingLeastSquaresStencil::calcFirstHaloDepth() const
{
    // 1. Step: get local depth
    const scalarField& V = mesh_.V();

    scalar procPatchArea = 0.0;
    scalar haloVol = 0.0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            const vectorField& Sf = pp.faceAreas();
            const labelUList& faceCells = pp.faceCells();

            forAll(faceCells, i)
            {
                const label cellID = faceCells[i];

                haloVol += V[cellID];
                procPatchArea += mag(Sf[i]);
            }
        }
    }

    const scalar localHaloDepth =
        (procPatchArea > 0 ? haloVol/procPatchArea : 0.0);

    // 2. Step: get list of neighbouring processors
    labelHashSet nbrSet;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>(pp);

            nbrSet.insert(ppp.neighbProcNo());
        }
    }

    const labelList nbrs =  nbrSet.toc();

    // 3. Step: point-to-point exchange of local depth with other processors
    List<scalar> neighbourHaloDepth(Pstream::nProcs(), 0.0);

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    forAll(nbrs, i)
    {
        UOPstream os(nbrs[i], pBufs);
        os << localHaloDepth;
    }
    pBufs.finishedSends();

    forAll(nbrs, i)
    {
        const label from = nbrs[i];

#ifdef OPENFOAM_COM
        if (pBufs.recvDataCount(from))
#endif
        {
            UIPstream is(from, pBufs);

            scalar lenN;
            is >> lenN;
            neighbourHaloDepth[from] = lenN;
        }
    }

    // Return averaged first halo depth per processor
    return neighbourHaloDepth;
}

labelList movingLeastSquaresStencil::checkProcessorOverlap
(
    const List<treeBoundBox>& allOwnedCellsBox,
    const treeBoundBox& ownedFacesBox
) const
{
    DynamicList<label> overlappingProcessor;

    for (label proc = 0; proc < Pstream::nProcs(); ++proc)
    {
        if (proc == Pstream::myProcNo())
        {
            continue;
        }
        if (allOwnedCellsBox[proc].overlaps(ownedFacesBox))
        {
            overlappingProcessor.append(proc);
        }
    }

    return overlappingProcessor.shrink();
}

treeBoundBox movingLeastSquaresStencil::calcOwnedCellsBox() const
{
    vector minPt(GREAT, GREAT, GREAT);
    vector maxPt(-GREAT, -GREAT, -GREAT);

    const pointField& C = mesh_.C().primitiveField();

    forAll(C, cellI)
    {
        minPt = min(minPt, C[cellI]);
        maxPt = max(maxPt, C[cellI]);
    }

    return treeBoundBox(minPt, maxPt);
}


treeBoundBox movingLeastSquaresStencil::calcOwnedFacesBox() const
{
    vector minPt(GREAT, GREAT, GREAT);
    vector maxPt(-GREAT, -GREAT, -GREAT);

    const pointField& faceCentres = mesh_.faceCentres();

    forAll(faceCentres, faceI)
    {
        minPt = min(minPt, faceCentres[faceI]);
        maxPt = max(maxPt, faceCentres[faceI]);
    }

    return treeBoundBox(minPt, maxPt);
}

List<labelList> movingLeastSquaresStencil::remoteCandidates
(
    const treeBoundBox& ownedFacesBox,
    const labelList& procToQuery
) const
{
    const vectorField& C = mesh_.C().primitiveField();

    List<labelList> remoteCandidatesPerProc;

    remoteCandidatesPerProc.setSize(Pstream::nProcs());

    forAll(remoteCandidatesPerProc, procI)
    {
        remoteCandidatesPerProc[procI].clear();
    }

    // Phase 1: Exchange ownedFacesBox between processors
    Map<treeBoundBox> incomingBoxesFromProc;
    {
        PstreamBuffers sBufs(Pstream::commsTypes::nonBlocking);

        forAll(procToQuery, i)
        {
            const label toProc = procToQuery[i];
            UOPstream os(toProc, sBufs);
            os << ownedFacesBox;
        }

        sBufs.finishedSends();

        forAll(procToQuery, i)
        {
            const label from = procToQuery[i];
#ifdef OPENFOAM_COM
            if (sBufs.recvDataCount(from) && from == Pstream::myProcNo())
#endif
            {
                continue;
            }
            UIPstream is(from, sBufs);
            treeBoundBox qb;
            is >> qb;

            incomingBoxesFromProc.insert(from, qb);
        }
    }

    // Phase 2: Mark local cells and send back global IDs
    {
        PstreamBuffers rBufs(Pstream::commsTypes::nonBlocking);

        forAllConstIter(Map<treeBoundBox>, incomingBoxesFromProc, it)
        {
            const label sender = it.key();
            const treeBoundBox& qb = it();

            labelHashSet usedBySender;
            forAll(C, cellI)
            {
                if (qb.contains(C[cellI]))
                {
                    usedBySender.insert(cellI);
                }
            }

            // Send back as a list
            labelList markedCells(usedBySender.size());
            label i = 0;
            forAllConstIter(labelHashSet, usedBySender, iter)
            {
                const label localCell = iter.key();
                markedCells[i++] = globalCells_.toGlobal(localCell);
            }

            UOPstream os(sender, rBufs);
            os << markedCells;
        }

        rBufs.finishedSends();


        // Phase 3: Recieve marked cells from  other processors
        forAll(procToQuery, i)
        {
            const label fromProc = procToQuery[i];
#ifdef OPENFOAM_COM
            if (rBufs.recvDataCount(fromProc))
#endif
            {
                continue;
            }
            UIPstream is(fromProc, rBufs);
            labelList lst;
            is >> lst;

            remoteCandidatesPerProc[fromProc].transfer(lst);
        }
    }

    return remoteCandidatesPerProc;
}


List<vectorField> movingLeastSquaresStencil::remoteCandidatesCellCentres
(
    const List<labelList>& remoteCandidates,
    const labelList& procToQuery
) const
{
    const vectorField& C = mesh_.C().primitiveField();

    List<vectorField> remoteCellCentres(Pstream::nProcs());

    forAll(remoteCellCentres, procI)
    {
        remoteCellCentres[procI].clear();
    }

    if (!Pstream::parRun())
    {
        return remoteCellCentres;
    }

    PstreamBuffers reqBufs(Pstream::commsTypes::nonBlocking);

    // Phase 1: send cell global IDs to each processor in contact
    forAll(procToQuery, i)
    {
        const label procI = procToQuery[i];

        if (!remoteCandidates[procI].empty())
        {
            UOPstream os(procI, reqBufs);
            os << remoteCandidates[procI];
        }
    }

    reqBufs.finishedSends();

    // Phase 2: respond with cell centre coordinates
    PstreamBuffers repBufs(Pstream::commsTypes::nonBlocking);

    forAll (procToQuery, i)
    {
        const label fromProc = procToQuery[i];

        if
        (
#ifdef OPENFOAM_COM
            reqBufs.recvDataCount(fromProc) ||
#endif
            fromProc == Pstream::myProcNo()
        )
        {
            continue;
        }

        UIPstream is(fromProc, reqBufs);

        labelList processorsCellID;
        is >> processorsCellID;

        vectorField centres(processorsCellID.size());

        forAll(processorsCellID, i)
        {
            const label globalCellID = processorsCellID[i];

            const label localCell = globalCells_.toLocal(globalCellID);

            if (localCell < 0 || localCell >= mesh_.nCells())
            {
                FatalErrorInFunction
                    << "Invalid global->local mapping: globaCellID="
                    << globalCellID << " localCell=" << localCell
                    << " on proc " << Pstream::myProcNo() << nl
                    << exit(FatalError);
            }

            centres[i] = C[localCell];
         }

         UOPstream os(fromProc, repBufs);
         os << centres;
    }

    repBufs.finishedSends();

    // Phase 3: Send back populated cell centres lists

    forAll(procToQuery, i)
    {
        const label p = procToQuery[i];

        if (remoteCandidates[p].empty())
        {
            continue;
        }

#ifdef OPENFOAM_COM
        if (repBufs.recvDataCount(p))
#endif
        {
            FatalErrorInFunction
                << "Did not receive centres from proc " << p << nl
                << exit(FatalError);
        }

        UIPstream is(p, repBufs);

        vectorField centres;
        is >> centres;

        if (centres.size() != remoteCandidates[p].size())
        {
            FatalErrorInFunction
                << "Centres reply size mismatch from proc " << p
                << ": got " << centres.size()
                << " expected " << remoteCandidates[p].size() << nl
                << exit(FatalError);
        }

        remoteCellCentres[p].transfer(centres);
    }

    return remoteCellCentres;
}

void  movingLeastSquaresStencil::calcProcessorCells() const
{
    procCellsPtr_.reset(new boolList(mesh_.nCells(), false));
    boolList& procCells = procCellsPtr_();

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            const labelUList& faceCells = pp.faceCells();
            forAll(faceCells, i)
            {
                procCells[faceCells[i]] = true;
            }
        }
    }
}

labelList movingLeastSquaresStencil::buildFacesStencil
(
    const label faceI,
    const List<labelList>& remoteCells,
    const List<vectorField>& remoteCellsCentres,
    const scalar relTol,
    const scalar overSampleFactor,
    const scalar searchExpansionFactor
) const
{
    // When this number of cell is reached, layer approach stops
    const label nbOfCandidates =
        label(std::ceil(overSampleFactor * scalar(N_)));

    // When building stencils, stencil is larger than N if there are cells with
    // equall distance. We use relTol to detect such cases. This results in
    // symmetric stencils on regular meshes.

    // Get potential candidates using complete layers looping.
    // Cell ID is in local indexing format.
    labelHashSet localCandidates;

    localCandidates.clear();

    labelHashSet currentLayer;
    labelHashSet nextLayer;

    const label ownerCell = mesh_.faceOwner()[faceI];
    currentLayer.insert(ownerCell);
    localCandidates.insert(ownerCell);

    if (faceI < mesh_.nInternalFaces())
    {
        const label nei = mesh_.faceNeighbour()[faceI];
        currentLayer.insert(nei);
        localCandidates.insert(nei);
    }

    const labelListList& cellCells = mesh_.cellCells();

    while
    (
        !currentLayer.empty()
     && localCandidates.size() < nbOfCandidates
    )
    {
        nextLayer.clear();

        forAllConstIter(labelHashSet, currentLayer, it)
        {
            const label c = it.key();
            const labelList& neighbours = cellCells[c];

            forAll(neighbours, i)
            {
                const label cellID = neighbours[i];
                if (localCandidates.insert(cellID))
                {
                    nextLayer.insert(cellID);
                }
            }
        }

        currentLayer.transfer(nextLayer);
    }

    // Above algorithm works only on owned faces since for remote cells we did
    // not transfered cellCells structure.

    // Algorithm: 1) Sort current candidates (distance based)
    //            2) Loop over existing candidates and check if some cell is
    //               at processor boundary. Looping is done up to N position.
    //            3) Get radius of the cell that is at position 1.3*N
    //            4) Check if any remote cell fits into this radius, in case it
    //               is inside this radius add it to the candidates list.
    //            5) Sort again and take first N cells into stencil

    // Phase 1: Build distance list for local candidates
    //          Using squared distance for efficiency

    const vector faceCentre = mesh_.faceCentres()[faceI];
    const vectorField& C = mesh_.C().primitiveField();;

    labelList localList(localCandidates.size());
    {
        label k = 0;
        forAllConstIter(labelHashSet, localCandidates, it)
        {
            localList[k++] = it.key();
        }
    }

    List<Tuple2<label, scalar>> localDist(localList.size());

    forAll(localList, i)
    {
        const label cI = localList[i];
        localDist[i] = Tuple2<label, scalar>(cI, magSqr(C[cI] - faceCentre));
    }

    Foam::sort
    (
        localDist,
        [](auto& A, auto& B)
        {
            return A.second() < B.second();
        }
    );

    // Phase 2: Check first N local candidates if internal cells is at processor
    //          boundary to activate looping over remote cells

    bool localStencil = true;
    const label nCheck = min(N_, localDist.size());

    for (label pos = 0; pos < nCheck; ++pos)
    {
        const label& cellID = localDist[pos].first();
        if (procCells()[cellID])
        {
            localStencil = false;
            break;
        }
    }

    // If stencil is local we can already return the stencil
    if (localStencil)
    {
        if ( localList.size() < N_ )
        {
            FatalErrorInFunction
                << "Local candidates for stencil have size of: "
                << localList.size()
                << " but required minimum  stencil size is larger: " << N_ << nl
                << "Increase the mesh size!"
                << exit(FatalError);
        }

        // When building stencils, stencil is larger than N if there are cells
        // with equall distance. We use relTol to detect such cases.
        // This results in symmetric stencils on regular meshes.
        const scalar cut = localDist[N_-1].second();
        const scalar cutTol = cut*(1.0 + relTol);

        label nPick = N_;
        while (nPick < localDist.size() && localDist[nPick].second() <= cutTol)
        {
            ++nPick;
        }

        labelList stencil(nPick);
        for (label i = 0; i < nPick; ++i)
        {
            const label localCellID = localDist[i].first();
            stencil[i] = globalCells_.toGlobal(localCellID);
        }

        return stencil;
    }

    // Phase 3: Radius from position posR = ceil(searchExpansionFactor*N)

    label posR = label(std::ceil(searchExpansionFactor * scalar(N_)));
    posR = min(max(posR, 0), localDist.size()-1);

    const scalar R2 = localDist[posR].second();

    DynamicList<Tuple2<label, scalar>> allDist;
    allDist.reserve((posR + 1) + 100);

    // Local cells within extended stencil radius squared R2
    for (label i = 0; i <= posR; ++i)
    {
        const label cellID = localDist[i].first();
        const label globalCellID = globalCells_.toGlobal(cellID);
        allDist.append
        (
            Tuple2<label, scalar>(globalCellID, localDist[i].second())
        );
    }

    // Phase 4: Add remote cells within radius squared R2
    forAll(remoteCells, p)
    {
        const labelList& procCells = remoteCells[p];
        const vectorField& procCellCentres = remoteCellsCentres[p];

        if (procCells.empty())
        {
            continue;
        }

        forAll(procCells, i)
        {
            const scalar d2 = magSqr(procCellCentres[i] - faceCentre);
            if (d2 <= R2)
            {
                allDist.append(Tuple2<label, scalar>(procCells[i], d2));
            }
        }
    }

    // Phase 5: Sort again local and remote cells together
    Foam::sort
    (
        allDist,
        [](const Tuple2<label, scalar>& A, const Tuple2<label, scalar>& B)
        {
            return A.second() < B.second();
        }
    );

    // Check minimum stencil size
    if ( allDist.size() < N_ )
    {
        FatalErrorInFunction
            << "Candidates for stencil have size of: " << allDist.size()
            << " but required minimum stencil size is larger: " << N_ << nl
            << "Increase the mesh size!"
            << exit(FatalError);
    }

    // Enlarge stencil to include cells with the same distance
    const scalar cut = allDist[N_-1].second();
    const scalar cutTol = cut*(1.0 + relTol);

    label nPick = N_;
    while (nPick < allDist.size() && allDist[nPick].second() <= cutTol)
    {
        ++nPick;
    }

    labelList faceStencil(nPick);
    for (label i = 0; i < nPick; ++i)
    {
        faceStencil[i] = allDist[i].first();
    }

    return faceStencil;
}

void movingLeastSquaresStencil::filterUnusedCandidates
(
    const List<labelList>& remoteCellsPerProc,
    const List<vectorField>& remoteCellsCentresPerProc
) const
{
    // Initialise storage
    remoteCellsPerProcPtr_.reset(new List<labelList>(Pstream::nProcs()));
    remoteCentresPerProcPtr_.reset(new List<vectorField>(Pstream::nProcs()));

    labelHashSet usedRemoteCells;

#ifdef OPENFOAM_COM
    label totalRemote = 0;
    forAll(remoteCellsPerProc, p)
    {
        totalRemote += remoteCellsPerProc[p].size();
    }
    usedRemoteCells.reserve(totalRemote);
#endif

    const CompactListList<label>& facesStencil = this->facesStencil();

    forAll(facesStencil, faceI)
    {
        const labelUList faceStencil = facesStencil[faceI];

        forAll(faceStencil, cellI)
        {
            const label gCellID = faceStencil[cellI];
            if (!globalCells_.isLocal(gCellID))
            {
                usedRemoteCells.insert(gCellID);
            }
        }
    }

    forAll(remoteCellsPerProc, procI)
    {
        const labelList& globalIDs = remoteCellsPerProc[procI];
        const vectorField& cellCentres = remoteCellsCentresPerProc[procI];

        if (globalIDs.empty())
        {
            continue;
        }

        if (globalIDs.size() != cellCentres.size())
        {
            FatalErrorInFunction
                << "Size for remote lists of cell centres and global IDs "
                << "are in mismatch for processor " << procI
                << exit(FatalError);
        }

        DynamicList<label> filteredGlobalIDs;
        DynamicList<vector> filteredCellCentres;
        filteredGlobalIDs.reserve(globalIDs.size());
        filteredCellCentres.reserve(globalIDs.size());

        forAll(globalIDs, i)
        {
            const label globalID = globalIDs[i];
            if (usedRemoteCells.found(globalID))
            {
                filteredGlobalIDs.append(globalID);
                filteredCellCentres.append(cellCentres[i]);
            }
        }

        remoteCellsPerProcPtr_()[procI].transfer(filteredGlobalIDs.shrink());
        remoteCentresPerProcPtr_()[procI].transfer(filteredCellCentres.shrink());
    }
}


void movingLeastSquaresStencil::calcFacesStencil() const
{
    // Get average first halo depth per processor boundary
    const List<scalar> neighbourHaloDepth = calcFirstHaloDepth();

    // We will take max halo depth and scale it to get multi halo depth
    const scalar scaledHaloDepth = max(neighbourHaloDepth) * haloDepthScale_;

    // Box for faces, augmented with expected multi-halo depth
    treeBoundBox ownedFacesBox = calcOwnedFacesBox();
#ifdef OPENFOAM_COM
    ownedFacesBox.grow(scaledHaloDepth);
#endif
#ifdef OPENFOAM_ORG
    ownedFacesBox.inflate(scaledHaloDepth);
#endif


    // Box for remote processors cells
    List<treeBoundBox> allOwnedCellsBox(Pstream::nProcs());
    allOwnedCellsBox[Pstream::myProcNo()] = calcOwnedCellsBox();

#ifdef OPENFOAM_COM
    Pstream::allGatherList(allOwnedCellsBox);
#endif
#ifdef OPENFOAM_ORG
    Pstream::gatherList(allOwnedCellsBox);
    Pstream::scatterList(allOwnedCellsBox);
#endif

    // Get processors that may have stencil cells, detected using overlap test
    labelList procToQuery =
       checkProcessorOverlap(allOwnedCellsBox, ownedFacesBox);

    // List of remote cells (per processor) written using global indexing
    List<labelList> remoteCells =
        remoteCandidates(ownedFacesBox, procToQuery);

    // Get cell centre for remote candidates
    List<vectorField> remoteCellsCentres =
        remoteCandidatesCellCentres(remoteCells, procToQuery);

    // Build face stencils using local and remote cells
    List<labelList> facesStencil(mesh_.nFaces());
    forAll(facesStencil, faceI)
    {
        facesStencil[faceI].clear();
    }

    // Build stencil for internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
    {
        facesStencil[faceI] =
            buildFacesStencil
            (
                faceI,
                remoteCells,
                remoteCellsCentres
            );
    }

    // Build stencil for boundary faces
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        // Skip empty faces
        if (isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        // Skip non-owned processor faces to avoid double counting
        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>(pp);

            if (!ppp.owner())
            {
                continue;
            }
        }

        forAll(pp, i)
        {
            const label faceI = pp.start() + i;

            facesStencil[faceI] =
                buildFacesStencil
                (
                    faceI,
                    remoteCells,
                    remoteCellsCentres
                );
        }
    }

    // Allocate compact storage  and store into it
    labelList sizes(facesStencil.size());
    forAll(facesStencil, faceI)
    {
        sizes[faceI] = facesStencil[faceI].size();
    }

    facesStencilPtr_.reset(new CompactListList<label>(sizes));
    CompactListList<label>& facesStencilRef = facesStencilPtr_();

    forAll(facesStencil, faceI)
    {
        const labelList& src = facesStencil[faceI];
        labelUList dst = facesStencilRef[faceI];

        forAll(src, j)
        {
            dst[j] = src[j];
        }
    }

    // Filter and store remote cells and cell centres
    filterUnusedCandidates(remoteCells, remoteCellsCentres);
}


void movingLeastSquaresStencil::calcCellsStencil(const scalar relTol) const
{
    // Prerequisites
    const CompactListList<label>& faceStencils = this->facesStencil();

    const List<labelList>& remoteCellsPerProc = this->remoteCellsPerProc();

    const List<vectorField>& remoteCentresPerProc =
        this->remoteCentresPerProc();

    // Collect remote cells and their location in one hash table
    label totalRemote = 0;
    forAll(remoteCellsPerProc, p)
    {
        totalRemote += remoteCellsPerProc[p].size();
    }

    HashTable<vector, label> remoteCellsData;
    remoteCellsData.resize(2*totalRemote + 1);

    forAll(remoteCellsPerProc, p)
    {
        const labelList& globalIDs = remoteCellsPerProc[p];
        const vectorField& cellCentres = remoteCentresPerProc[p];

        forAll(globalIDs, i)
        {
            remoteCellsData.insert(globalIDs[i], cellCentres[i]);
        }
    }

    // Loop over cells, collect cells from face stencils, sort them by distance
    // and construct cell stencil by taking first Nc_ cells
    const vectorField& C = mesh_.C();
    const cellList& cellFaces = mesh_.cells();

    List<labelList> cellStencils(mesh_.nCells());

    DynamicList<Tuple2<label, scalar>> dist;
    dist.reserve(label(std::ceil(1.5 * scalar(N_))));

    labelHashSet candidates;
    for (label cellI = 0; cellI < mesh_.nCells(); ++cellI)
    {
        candidates.clear();

        const labelList& faces = cellFaces[cellI];

        // Loop over cell faces and include stencil cells to candidates
        forAll(faces, fI)
        {
            const label faceID = faces[fI];
            const labelUList faceStencil = faceStencils[faceID];

            forAll(faceStencil, j)
            {
                candidates.insert(faceStencil[j]);
            }
        }

        // Build distance list (squared distance to cell centre)
        dist.clear();
        if (dist.capacity() < candidates.size())
        {
            // Do reserve only if needed
            dist.reserve(candidates.size());
        }

        const vector& cellCentre = C[cellI];

        forAllConstIter(labelHashSet, candidates, iter)
        {
            scalar d2 = GREAT;

            const label globalID = iter.key();

            if (globalCells_.isLocal(globalID))
            {
                const label localID = globalCells_.toLocal(globalID);
                d2 = magSqr(C[localID] - cellCentre);
            }
            else
            {
                const auto it = remoteCellsData.find(globalID);
                d2 = magSqr(it() - cellCentre);
            }

            dist.append(Tuple2<label, scalar>(globalID, d2));
        }

        // Sort cells based on squared distance
        Foam::sort
        (
            dist,
            [](const Tuple2<label, scalar>& A, const Tuple2<label, scalar>& B)
            {
                return A.second() < B.second();
            }
        );

        if (dist.size() < Nc_)
        {
            FatalErrorInFunction
                << "Cell " << cellI << " has only " << dist.size()
                << " candidates but requred minimum number is " << Nc_
                << exit(FatalError);
        }

        // Enlarge stencil to include cells with the same distance
        const scalar cut = dist[Nc_-1].second();
        const scalar cutTol = cut*(1.0 + relTol);

        label nPick = Nc_;
        while (nPick < dist.size() && dist[nPick].second() <= cutTol)
        {
            ++nPick;
        }

        labelList st(nPick);
        for (label i = 0; i < nPick; ++i)
        {
            st[i] = dist[i].first();
        }

        cellStencils[cellI].transfer(st);
    }

    // Allocate compact storage  and store into it
    labelList sizes(cellStencils.size());
    forAll(cellStencils, i)
    {
        sizes[i] = cellStencils[i].size();
    }

    cellsStencilPtr_.reset(new CompactListList<label>(sizes));

    CompactListList<label>& cellsStencilRef = cellsStencilPtr_();

    forAll(cellStencils, c)
    {
        const labelList& src = cellStencils[c];
        labelUList dst = cellsStencilRef[c];

        forAll(src, j)
        {
            dst[j] = src[j];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


movingLeastSquaresStencil::movingLeastSquaresStencil
(
    const fvMesh& mesh,
    const scalar haloDepthScale,
    const label N,
    const label Nc
)
:
    mesh_(mesh),
    globalCells_(mesh.nCells()),
    haloDepthScale_(haloDepthScale),
    N_(N),
    Nc_(Nc),
    procCellsPtr_(),
    facesStencilPtr_(),
    cellsStencilPtr_()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const CompactListList<label>& movingLeastSquaresStencil::facesStencil() const
{
    if (!facesStencilPtr_.valid())
    {
         calcFacesStencil();
    }

    return facesStencilPtr_();
}


const CompactListList<label>& movingLeastSquaresStencil::cellsStencil() const
{
    // Cells stencils are constructed from faces stencils
    if (!facesStencilPtr_.valid())
    {
         calcFacesStencil();
    }

    if (!cellsStencilPtr_.valid())
    {
        calcCellsStencil();
    }

    return cellsStencilPtr_();
}


const List<labelList>& movingLeastSquaresStencil::remoteCellsPerProc() const
{
    if (!remoteCellsPerProcPtr_.valid())
    {
         calcFacesStencil();
    }

    return remoteCellsPerProcPtr_();
}


const List<vectorField>& movingLeastSquaresStencil::remoteCentresPerProc() const
{
    if (!remoteCentresPerProcPtr_.valid())
    {
        calcFacesStencil();
    }

    return remoteCentresPerProcPtr_();
}


const boolList& movingLeastSquaresStencil::procCells() const
{
    if (!procCellsPtr_.valid())
    {
        calcProcessorCells();
    }

    return procCellsPtr_();
}


void Foam::movingLeastSquaresStencil::clear() const
{
    facesStencilPtr_.clear();
    cellsStencilPtr_.clear();
    remoteCellsPerProcPtr_.clear();
    remoteCentresPerProcPtr_.clear();
    procCellsPtr_.clear();
}

} // End namespace Foam

// ************************************************************************* //
