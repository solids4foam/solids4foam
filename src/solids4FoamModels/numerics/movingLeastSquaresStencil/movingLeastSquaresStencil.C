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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(movingLeastSquaresStencil, 0);

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //




// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

List<scalar> movingLeastSquaresStencil::calcFirstHaloDepth() const
{
    // 1. Step: get local depth
    const scalarField& V = mesh_.V();

    labelHashSet haloCells;
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

                haloCells.insert(cellID);
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
        if (pBufs.recvDataCount(from))
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
)
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

    const pointField& C = mesh_.C();

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
    const vectorField& C = mesh_.C();

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

        for (label p = 0; p < Pstream::nProcs(); ++p)
        {
            if (!sBufs.recvDataCount(p) || p == Pstream::myProcNo())
            {
                continue;
            }
            UIPstream is(p, sBufs);
            treeBoundBox qb;
            is >> qb;

            incomingBoxesFromProc.insert(p, qb);
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
            if (!rBufs.recvDataCount(fromProc))
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
    const vectorField& C = mesh_.C();

    List<vectorField> remoteCellCentres(Pstream::nProcs());

    forAll(remoteCellCentres, procI)
    {
        remoteCellCentres[procI].clear();
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

    for (label fromProc = 0; fromProc < Pstream::nProcs(); ++fromProc)
    {
        if
        (
            fromProc == Pstream::myProcNo()
         || !reqBufs.recvDataCount(fromProc)
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

        if (!repBufs.recvDataCount(p))
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

labelList movingLeastSquaresStencil::buildFacesStencil
(
    const label faceI,
    const List<labelList>& remoteCells,
    const List<vectorField>& remoteCellsCentres
) const
{

    labelList  faceStencil;

    return faceStencil;
}


void movingLeastSquaresStencil::calcFacesStencil()
{
    // Get average first halo depth per processor boundary
    const List<scalar> neighbourHaloDepth = calcFirstHaloDepth();

    // We will take max halo depth and scale it to get multi halo depth
    const scalar scaledHaloDepth = max(neighbourHaloDepth) * haloDepthScale_;

    // Box for faces, augmented with expected multi-halo depth
    treeBoundBox ownedFacesBox = calcOwnedFacesBox();
    ownedFacesBox.grow(scaledHaloDepth);

    // Box for remote processors cells
    List<treeBoundBox> allOwnedCellsBox(Pstream::nProcs());
    allOwnedCellsBox[Pstream::myProcNo()] = calcOwnedCellsBox();
    Pstream::allGatherList(allOwnedCellsBox);

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

    // labelList sizes(facesStencil.size());
    // forAll(facesStencil, faceI)
    // {
    //     sizes[faceI] = facesStencil[faceI].size();
    // }


    // // Store face stencils using contigous compact storage
    // facesStencil_ = CompactListList<label>(facesStencil);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


movingLeastSquaresStencil::movingLeastSquaresStencil
(
    const fvMesh& mesh,
    const scalar& haloDepthScale,
    const label& N
)
:
    mesh_(mesh),
    globalCells_(mesh.nCells()),
    haloDepthScale_(haloDepthScale),
    N_(N),
    facesStencil_(),
    cellsStencil_()
{
    // Build face stencils at construction
    calcFacesStencil();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
