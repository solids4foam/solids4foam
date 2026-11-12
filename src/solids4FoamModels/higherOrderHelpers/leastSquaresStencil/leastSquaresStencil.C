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

#include "leastSquaresStencil.H"
#include "fvMesh.H"
#ifdef FOAMEXTEND
    #include "Pstream.H"
#else
    #include "PstreamBuffers.H"
#endif
#include "processorPolyPatch.H"
#include "volFields.H"
#include "Tuple2.H"
#include "cellZoneInterface.H"
#include "regionSplit.H"
#include <cmath>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(leastSquaresStencil, 0);


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

List<scalar> leastSquaresStencil::calcFirstHaloDepth() const
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

    // 3. Step: exchange local depth with neighbouring processors
    List<scalar> neighbourHaloDepth(Pstream::nProcs(), 0.0);

    if (!Pstream::parRun())
    {
        return neighbourHaloDepth;
    }

#ifdef FOAMEXTEND
    List<scalarField> sendHaloDepth(Pstream::nProcs());
    forAll(nbrs, i)
    {
        sendHaloDepth[nbrs[i]].setSize(1, localHaloDepth);
    }

    List<scalarField> receivedHaloDepth;
    labelListList sizes;
    Pstream::exchange<scalarField, scalar>
    (
        sendHaloDepth,
        receivedHaloDepth,
        sizes
    );

    forAll(nbrs, i)
    {
        const label from = nbrs[i];

        if (receivedHaloDepth[from].size() != 1)
        {
            FatalErrorInFunction
                << "Expected one halo depth from processor " << from
                << " but received " << receivedHaloDepth[from].size()
                << abort(FatalError);
        }

        neighbourHaloDepth[from] = receivedHaloDepth[from][0];
    }
#else
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
        if (!pBufs.recvDataCount(from))
        {
            continue;
        }
#endif
        UIPstream is(from, pBufs);

        scalar lenN;
        is >> lenN;
        neighbourHaloDepth[from] = lenN;
    }
#endif

    // Return averaged first halo depth per processor
    return neighbourHaloDepth;
}

labelList leastSquaresStencil::checkProcessorOverlap
(
    const List<orientedBoundBox>& allOwnedCellsBox,
    const List<orientedBoundBox>& allOwnedFacesBox
) const
{
    DynamicList<label> overlappingProcessor;
    const label myProcNo = Pstream::myProcNo();

    for (label proc = 0; proc < Pstream::nProcs(); ++proc)
    {
        if (proc == myProcNo)
        {
            continue;
        }

        // Both ranks must evaluate the same ordered SAT operations so that
        // round-off cannot produce different communication partner lists.
        const label first = min(myProcNo, proc);
        const label second = max(myProcNo, proc);

        if
        (
            allOwnedCellsBox[first].overlaps(allOwnedFacesBox[second])
         || allOwnedFacesBox[first].overlaps(allOwnedCellsBox[second])
        )
        {
            overlappingProcessor.append(proc);
        }
    }

    return overlappingProcessor.shrink();
}

orientedBoundBox leastSquaresStencil::calcOwnedCellsBox() const
{
    return orientedBoundBox(primitiveField(mesh_.C()));
}


orientedBoundBox leastSquaresStencil::calcOwnedFacesBox() const
{
    return orientedBoundBox(mesh_.faceCentres());
}

List<labelList> leastSquaresStencil::remoteCandidates
(
    const List<orientedBoundBox>& allOwnedFacesBox,
    const labelList& procToQuery
) const
{
    const vectorField& C = primitiveField(mesh_.C());

    List<labelList> remoteCandidatesPerProc;

    remoteCandidatesPerProc.setSize(Pstream::nProcs());

    forAll(remoteCandidatesPerProc, procI)
    {
        remoteCandidatesPerProc[procI].clear();
    }

    if (!Pstream::parRun())
    {
        return remoteCandidatesPerProc;
    }

    // Mark local cells inside each partner's gathered face box. Partner lists
    // are symmetric by construction, so the boxes need not be exchanged again
    List<labelList> sendCandidates(Pstream::nProcs());
    forAll(procToQuery, i)
    {
        const label sender = procToQuery[i];
        const orientedBoundBox& queryBox = allOwnedFacesBox[sender];
        DynamicList<label> markedCells;

        forAll(C, cellI)
        {
            if (queryBox.contains(C[cellI]))
            {
                markedCells.append(globalCells_.toGlobal(cellI));
            }
        }

        sendCandidates[sender].transfer(markedCells.shrink());
    }

#ifdef FOAMEXTEND
    labelListList candidateSizes;
    Pstream::exchange<labelList, label>
    (
        sendCandidates,
        remoteCandidatesPerProc,
        candidateSizes
    );
#else
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(procToQuery, i)
    {
        const label toProc = procToQuery[i];
        UOPstream os(toProc, pBufs);
        os << sendCandidates[toProc];
    }

    pBufs.finishedSends();

    forAll(procToQuery, i)
    {
        const label fromProc = procToQuery[i];
#ifdef OPENFOAM_COM
        if (!pBufs.recvDataCount(fromProc))
        {
            continue;
        }
#endif
        UIPstream is(fromProc, pBufs);
        is >> remoteCandidatesPerProc[fromProc];
    }
#endif

    return remoteCandidatesPerProc;
}


List<vectorField> leastSquaresStencil::remoteCandidatesCellCentres
(
    const List<labelList>& remoteCandidates,
    const labelList& procToQuery
) const
{
    return remoteCandidatesValues
    (
        remoteCandidates,
        procToQuery,
        primitiveField(mesh_.C())
    );
}


void leastSquaresStencil::calcCellRegion() const
{
    const labelList material(cellMaterialID(mesh_));

    bool multiMaterial = false;
    forAll(material, cellI)
    {
        if (material[cellI] != -1)
        {
            multiMaterial = true;
            break;
        }
    }
    reduce(multiMaterial, orOp<bool>());

    if (!multiMaterial)
    {
        // Single material: nothing is filtered, stencils are unchanged
        cellRegionPtr_.reset(new labelList());
        return;
    }

    // A material interface on a processor boundary would need each side's
    // reconstruction on the other processor; refuse it
    const List<boolList> coupledInterface(cellZoneInterfaceCoupled(mesh_));

    label nCoupledInterfaces = 0;
    forAll(coupledInterface, patchI)
    {
        forAll(coupledInterface[patchI], i)
        {
            if (coupledInterface[patchI][i])
            {
                ++nCoupledInterfaces;
            }
        }
    }
    reduce(nCoupledInterfaces, sumOp<label>());

    if (nCoupledInterfaces > 0)
    {
        FatalErrorInFunction
            << nCoupledInterfaces << " faces between two materials lie on a "
            << "processor boundary." << nl
            << "    High-order reconstruction with more than one material "
            << "needs every material interface inside one processor:"
            << " constrain the decomposition (singleProcessorFaceSets) so that"
            << " the cells on both sides of each interface share a processor."
            << exit(FatalError);
    }

    // The connected regions of each material: cells are connected through
    // every face but those between two materials, which are all internal
    boolList blockedFace(mesh_.nFaces(), false);
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
    {
        blockedFace[faceI] = (material[own[faceI]] != material[nei[faceI]]);
    }

    const regionSplit regions(mesh_, blockedFace);

    cellRegionPtr_.reset(new labelList(regions));
}


const labelList& leastSquaresStencil::cellRegion() const
{
    if (!cellRegionPtr_.valid())
    {
        calcCellRegion();
    }

    return cellRegionPtr_();
}


void leastSquaresStencil::calcRemoteRegions
(
    const List<labelList>& remoteCells,
    const labelList& procToQuery
) const
{
    remoteRegionPtr_.reset(new Map<label>());

    // Collective: every rank takes this branch, or none does
    if (cellRegion().empty())
    {
        return;
    }

    const List<Field<label>> regions =
        remoteCandidatesValues(remoteCells, procToQuery, cellRegion());

    Map<label>& remoteRegion = autoPtrRef(remoteRegionPtr_);

    forAll(remoteCells, procI)
    {
        forAll(remoteCells[procI], i)
        {
            remoteRegion.set(remoteCells[procI][i], regions[procI][i]);
        }
    }
}


label leastSquaresStencil::regionOfGlobalCell(const label globalCellID) const
{
    const labelList& region = cellRegion();

    if (region.empty())
    {
        return -1;
    }

    if (globalCells_.isLocal(globalCellID))
    {
        return region[globalCells_.toLocal(globalCellID)];
    }

    if (remoteRegionPtr_.valid())
    {
        const Map<label>& remoteRegion = remoteRegionPtr_();
        const Map<label>::const_iterator iter =
            remoteRegion.find(globalCellID);

        if (iter != remoteRegion.end())
        {
            return iter();
        }
    }

    FatalErrorInFunction
        << "No material region for remote cell " << globalCellID
        << abort(FatalError);

    return -1;
}


void leastSquaresStencil::calcFaceRegion() const
{
    faceRegionPtr_.reset(new labelList(mesh_.nFaces(), -1));
    labelList& faceRegion = autoPtrRef(faceRegionPtr_);

    const labelList& region = cellRegion();

    if (region.empty())
    {
        return;
    }

    // Every face: its owner's region. Inside one material that is the face's
    // region; a face between two materials is reconstructed from each side
    // separately - from the owner's region here, and from the neighbour's in
    // interfaceNeighbourFacesStencil - so that no reconstruction mixes two
    // materials and neither material is preferred
    const labelUList& own = mesh_.faceOwner();

    forAll(faceRegion, faceI)
    {
        faceRegion[faceI] = region[own[faceI]];
    }
}


void  leastSquaresStencil::calcProcessorCells() const
{
    procCellsPtr_.reset(new boolList(mesh_.nCells(), false));
    boolList& procCells = autoPtrRef(procCellsPtr_);

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


labelList leastSquaresStencil::buildFacesStencil
(
    const label faceI,
    const List<labelList>& remoteCells,
    const List<vectorField>& remoteCellsCentres,
    const scalar relTol,
    const scalar overSampleFactor,
    const scalar searchExpansionFactor,
    const label stencilRegionIn
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

    // With more than one material only cells of the face's material region
    // enter the stencil (-1: single material, no filtering)
    if (!faceRegionPtr_.valid())
    {
        calcFaceRegion();
    }
    const label stencilRegion =
        stencilRegionIn >= 0
      ? stencilRegionIn
      : faceRegionPtr_()[faceI];
    const auto sameMaterial = [&](const label localCellID)
    {
        return
            stencilRegion == -1
         || cellRegion()[localCellID] == stencilRegion;
    };

    const label ownerCell = mesh_.faceOwner()[faceI];
    if (sameMaterial(ownerCell))
    {
        currentLayer.insert(ownerCell);
        localCandidates.insert(ownerCell);
    }

    if (faceI < mesh_.nInternalFaces())
    {
        const label nei = mesh_.faceNeighbour()[faceI];
        if (sameMaterial(nei))
        {
            currentLayer.insert(nei);
            localCandidates.insert(nei);
        }
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
                if (!sameMaterial(cellID))
                {
                    continue;
                }
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
    const vectorField& C = primitiveField(mesh_.C());

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
        []
        (
            const Tuple2<label, scalar>& A,
            const Tuple2<label, scalar>& B
        )
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

    const label unboundedPosR =
        label(std::ceil(searchExpansionFactor * scalar(N_)));

    label posR = min(max(unboundedPosR, 0), localDist.size()-1);

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
            if
            (
                stencilRegion != -1
             && regionOfGlobalCell(procCells[i]) != stencilRegion
            )
            {
                continue;
            }
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

    // The local radius is under-estimated if the processor partition contains
    // too few local cells; use all exchanged remote candidates before failing.
    if (allDist.size() < N_ && unboundedPosR >= localDist.size())
    {
        allDist.clear();
        allDist.reserve(localDist.size() + 100);

        forAll(localDist, i)
        {
            const label cellID = localDist[i].first();
            const label globalCellID = globalCells_.toGlobal(cellID);
            allDist.append
            (
                Tuple2<label, scalar>(globalCellID, localDist[i].second())
            );
        }

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
                if
                (
                    stencilRegion != -1
                 && regionOfGlobalCell(procCells[i]) != stencilRegion
                )
                {
                    continue;
                }
                allDist.append
                (
                    Tuple2<label, scalar>
                    (
                        procCells[i],
                        magSqr(procCellCentres[i] - faceCentre)
                    )
                );
            }
        }

        Foam::sort
        (
            allDist,
            [](const Tuple2<label, scalar>& A, const Tuple2<label, scalar>& B)
            {
                return A.second() < B.second();
            }
        );
    }

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


void leastSquaresStencil::filterUnusedCandidates
(
    const List<labelList>& remoteCellsPerProc,
    const List<vectorField>& remoteCellsCentresPerProc
) const
{
    // Initialise storage
    remoteCellsPerProcPtr_.reset(new List<labelList>(Pstream::nProcs()));
    remoteCentresPerProcPtr_.reset(new List<vectorField>(Pstream::nProcs()));
    List<labelList>& storedRemoteCells = autoPtrRef(remoteCellsPerProcPtr_);
    List<vectorField>& storedRemoteCentres =
        autoPtrRef(remoteCentresPerProcPtr_);

    labelHashSet usedRemoteCells;

#ifdef OPENFOAM_COM
    label totalRemote = 0;
    forAll(remoteCellsPerProc, p)
    {
        totalRemote += remoteCellsPerProc[p].size();
    }
    usedRemoteCells.reserve(totalRemote);
#endif

    // Remote cells used by any face stencil, the neighbour-side interface
    // stencils included
    for (label sideI = 0; sideI < 2; ++sideI)
    {
        auto& facesStencil =
            compactListListCRef
            (
                sideI == 0
              ? this->facesStencil()
              : this->interfaceNeighbourFacesStencil()
            );

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

        storedRemoteCells[procI].transfer(filteredGlobalIDs.shrink());
        storedRemoteCentres[procI].transfer(filteredCellCentres.shrink());
    }
}


namespace
{

// Pack a list of lists into compact storage
autoPtr<CompactListList<label>> toCompactListList(const List<labelList>& lists)
{
    labelList sizes(lists.size());
    forAll(lists, i)
    {
        sizes[i] = lists[i].size();
    }

    autoPtr<CompactListList<label>> resultPtr
    (
        new CompactListList<label>(sizes)
    );
    CompactListList<label>& result = autoPtrRef(resultPtr);

    forAll(lists, i)
    {
        const labelList& src = lists[i];
        labelUList dst = result[i];

        forAll(src, j)
        {
            dst[j] = src[j];
        }
    }

    return resultPtr;
}

} // End anonymous namespace


void leastSquaresStencil::calcFacesStencil() const
{
    // Get average first halo depth per processor boundary
    const List<scalar> neighbourHaloDepth = calcFirstHaloDepth();

    // We will take max halo depth and scale it to get multi halo depth
    const scalar scaledHaloDepth = max(neighbourHaloDepth) * haloDepthScale_;

    // Absolute halo distance, independent of the partition box diagonal.
    // This replaces the former fractional inflate() call on OpenFOAM.org.
    orientedBoundBox ownedFacesBox = calcOwnedFacesBox();
    ownedFacesBox.grow(scaledHaloDepth);

    // Box for remote processors cells
    List<orientedBoundBox> allOwnedCellsBox(Pstream::nProcs());
    allOwnedCellsBox[Pstream::myProcNo()] = calcOwnedCellsBox();

    List<orientedBoundBox> allOwnedFacesBox(Pstream::nProcs());
    allOwnedFacesBox[Pstream::myProcNo()] = ownedFacesBox;

#ifdef OPENFOAM_COM
    Pstream::allGatherList(allOwnedCellsBox);
    Pstream::allGatherList(allOwnedFacesBox);
#else
    Pstream::gatherList(allOwnedCellsBox);
    Pstream::scatterList(allOwnedCellsBox);
    Pstream::gatherList(allOwnedFacesBox);
    Pstream::scatterList(allOwnedFacesBox);
#endif

    // Get processors that may have stencil cells, detected using overlap test
    labelList procToQuery =
       checkProcessorOverlap(allOwnedCellsBox, allOwnedFacesBox);

    // List of remote cells (per processor) written using global indexing
    List<labelList> remoteCells =
        remoteCandidates(allOwnedFacesBox, procToQuery);

    // Get cell centre for remote candidates
    List<vectorField> remoteCellsCentres =
        remoteCandidatesCellCentres(remoteCells, procToQuery);

    // With more than one material, the candidates' material regions
    calcRemoteRegions(remoteCells, procToQuery);

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

    // The neighbour side of the faces between two materials, built here so
    // that their remote cells are kept by filterUnusedCandidates below
    List<labelList> interfaceNeighbourStencil(mesh_.nFaces());

    if (cellRegion().size())
    {
        const labelList& region = cellRegion();
        const labelUList& own = mesh_.faceOwner();
        const labelUList& nei = mesh_.faceNeighbour();

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
        {
            const label neiRegion = region[nei[faceI]];

            if (neiRegion != region[own[faceI]])
            {
                interfaceNeighbourStencil[faceI] =
                    buildFacesStencil
                    (
                        faceI,
                        remoteCells,
                        remoteCellsCentres,
                        1e-8,
                        1.5,
                        1.3,
                        neiRegion
                    );
            }
        }
    }

    facesStencilPtr_.reset(toCompactListList(facesStencil).ptr());
    interfaceNeighbourFacesStencilPtr_.reset
    (
        toCompactListList(interfaceNeighbourStencil).ptr()
    );

    // Filter and store remote cells and cell centres
    filterUnusedCandidates(remoteCells, remoteCellsCentres);
}


void leastSquaresStencil::calcCellsStencil(const scalar relTol) const
{
    // Prerequisites
    auto& faceStencils =
        compactListListCRef(this->facesStencil());

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
                // With more than one material a face stencil on a material
                // interface belongs to one side only; keep this cell's own
                if
                (
                    regionOfGlobalCell(faceStencil[j])
                 == regionOfGlobalCell(globalCells_.toGlobal(cellI))
                )
                {
                    candidates.insert(faceStencil[j]);
                }
            }
        }

        // Remove the current cell itself from candidates. The reconstruction
        // scheme handles its contribution separately, so the stencil contains
        // neighbours only.
        const label globalCellID = globalCells_.toGlobal(cellI);

        if (candidates.found(globalCellID))
        {
            candidates.erase(globalCellID);
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

        // Required numnber of cell is for one smaller becouse of cell centre
        // which is also part of the stencil
        const label nReq = Nc_ - 1;

        if (dist.size() < nReq)
        {
            FatalErrorInFunction
                << "Cell " << cellI << " has only " << dist.size()
                << " candidates but requred minimum number is " << Nc_
                << exit(FatalError);
        }

        // Enlarge stencil to include cells with the same distance
        const scalar cut = dist[nReq-1].second();
        const scalar cutTol = cut*(1.0 + relTol);

        label nPick = nReq;
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

    CompactListList<label>& cellsStencilRef = autoPtrRef(cellsStencilPtr_);

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


leastSquaresStencil::leastSquaresStencil
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
    cellsStencilPtr_(),
    remoteCellsPerProcPtr_(),
    remoteCentresPerProcPtr_(),
    remoteCentresMapPtr_(),
    remoteCellLocationPtr_()
{
    if (!(haloDepthScale_ > 0) || !std::isfinite(haloDepthScale_))
    {
        FatalErrorInFunction
            << "haloDepthScale must be finite and greater than zero, found "
            << haloDepthScale_ << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<labelField> leastSquaresStencil::cellFacesStencilSize() const
{
    const fvMesh& mesh = mesh_;

    tmp<labelField> tpnf(new labelField(mesh.nCells()));
    labelField& pnf = tmpRef(tpnf);

    // Get each face stencil
    auto& facesStencil =
        compactListListCRef(this->facesStencil());

    // Loop over cell faces stencils and merge thems into one
    forAll(mesh.cells(), cellI)
    {
        labelHashSet cellStencil;

        const labelList& cellFaces = mesh.cells()[cellI];

        forAll(cellFaces, faceI)
        {
            label faceID = cellFaces[faceI];

            // Loop over face stencil and store face stencil
            forAll(facesStencil[faceID], i)
            {
                cellStencil.insert(facesStencil[faceID][i]);
            }
        }

        pnf[cellI] = cellStencil.size();
    }

    return tpnf;
}


const CompactListList<label>& leastSquaresStencil::facesStencil() const
{
    if (facesStencilPtr_.empty())
    {
        calcFacesStencil();
    }

    return autoPtrRef(facesStencilPtr_);
}


const CompactListList<label>&
leastSquaresStencil::interfaceNeighbourFacesStencil() const
{
    if (interfaceNeighbourFacesStencilPtr_.empty())
    {
        calcFacesStencil();
    }

    return autoPtrRef(interfaceNeighbourFacesStencilPtr_);
}


const CompactListList<label>& leastSquaresStencil::cellsStencil() const
{
    // Cells stencils are constructed from faces stencils
    if (facesStencilPtr_.empty())
    {
        calcFacesStencil();
    }

    if (cellsStencilPtr_.empty())
    {
        calcCellsStencil();
    }

    return autoPtrRef(cellsStencilPtr_);
}


const List<labelList>& leastSquaresStencil::remoteCellsPerProc() const
{
    if (remoteCellsPerProcPtr_.empty())
    {
        calcFacesStencil();
    }

    return autoPtrRef(remoteCellsPerProcPtr_);
}


const List<vectorField>& leastSquaresStencil::remoteCentresPerProc() const
{
    if (remoteCentresPerProcPtr_.empty())
    {
        calcFacesStencil();
    }

    return autoPtrRef(remoteCentresPerProcPtr_);
}


const boolList& leastSquaresStencil::procCells() const
{
    if (procCellsPtr_.empty())
    {
        calcProcessorCells();
    }

    return autoPtrRef(procCellsPtr_);
}


const Map<vector>& leastSquaresStencil::remoteCentresMap() const
{
    if (!remoteCentresMapPtr_.valid())
    {
        const labelListList& remoteCells = this->remoteCellsPerProc();
        const List<vectorField>& remoteCentres = this->remoteCentresPerProc();

        // Get number of entries to avoid dynamic allocation
        label nRemote = 0;
        forAll(remoteCells, procI)
        {
            nRemote += remoteCells[procI].size();
        }
        remoteCentresMapPtr_.set(new Map<vector>(2*nRemote));

        Map<vector>& remoteCentresMap = autoPtrRef(remoteCentresMapPtr_);

        // Fill map
        forAll(remoteCells, procI)
        {
            const labelList& curRemoteCells = remoteCells[procI];
            const vectorField& curRemoteCentres = remoteCentres[procI];

            if (curRemoteCells.size() != curRemoteCentres.size())
            {
                FatalErrorInFunction
                    << "remoteCellsPerProc and remoteCentresPerProc sizes do "
                    << "not match for processor " << procI << nl
                    << "Number of remote cells = " << curRemoteCells.size()
                    << nl << "Number of remote centres = "
                    << curRemoteCentres.size()
                    << abort(FatalError);
            }

            forAll(curRemoteCells, i)
            {
                remoteCentresMap.insert(curRemoteCells[i], curRemoteCentres[i]);
            }
        }
    }

    return autoPtrRef(remoteCentresMapPtr_);
}


const globalIndex& leastSquaresStencil::globalCells() const
{
    return globalCells_;
}


const Map<FixedList<label, 2>>&
leastSquaresStencil::remoteCellLocation() const
{
    if (!remoteCellLocationPtr_.valid())
    {
        const List<labelList>& remoteCells = remoteCellsPerProc();

        label nRemote = 0;
        forAll(remoteCells, procI)
        {
            nRemote += remoteCells[procI].size();
        }

        remoteCellLocationPtr_.set
        (
            new Map<FixedList<label, 2>>(2*nRemote)
        );
        Map<FixedList<label, 2>>& locMap =
            autoPtrRef(remoteCellLocationPtr_);

        forAll(remoteCells, procI)
        {
            const labelList& thisProcRemoteCells = remoteCells[procI];

            forAll(thisProcRemoteCells, i)
            {
                FixedList<label, 2> loc;
                loc[0] = procI;
                loc[1] = i;

                locMap.insert(thisProcRemoteCells[i], loc);
            }
        }
    }

    return autoPtrRef(remoteCellLocationPtr_);
}


void leastSquaresStencil::clear() const
{
    facesStencilPtr_.clear();
    cellsStencilPtr_.clear();
    remoteCellsPerProcPtr_.clear();
    remoteCentresPerProcPtr_.clear();
    procCellsPtr_.clear();
    remoteCentresMapPtr_.clear();
    remoteCellLocationPtr_.clear();
    cellRegionPtr_.clear();
    remoteRegionPtr_.clear();
    faceRegionPtr_.clear();
    interfaceNeighbourFacesStencilPtr_.clear();
}

} // End namespace Foam

// ************************************************************************* //
