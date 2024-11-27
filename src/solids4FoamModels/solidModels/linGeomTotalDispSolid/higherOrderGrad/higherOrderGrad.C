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

#include "higherOrderGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "triangle.H"
#include "triFace.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(higherOrderGrad, 0);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void higherOrderGrad::makeGlobalCellStencils() const
{
    InfoInFunction
        << "start" << endl;

    if (globalCellStencilsPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set" << exit(FatalError);
    }

    // References from brevity and efficiency
    const fvMesh& mesh = mesh_;
    const label nCells = mesh.nCells();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const cellList& cells = mesh.cells();

    // Prepare and store processor neighbour face global cells, i.e. global cell
    // indices of the cells across the processor patches
    labelListList procPatchNeiGlobalCellIDs(mesh.boundaryMesh().size());
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

        // Send global cell IDs for each processor patch
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UOPstream toNeighbProc(procPatch.neighbProcNo(), pBufs);

                toNeighbProc
                    << globalCells_.toGlobal(pp.faceCells());
            }
        }

        pBufs.finishedSends(); // no-op for blocking

        // Receive data
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UIPstream fromNeighbProc(procPatch.neighbProcNo(), pBufs);

                procPatchNeiGlobalCellIDs[patchI].setSize(pp.size());

                fromNeighbProc
                    >> procPatchNeiGlobalCellIDs[patchI];
            }
        }
    }

    // Initialize visited sets and frontiers for all cells
    // For each local cell, we will find and store all cells in the stecil as
    // global cell indices
    List<labelHashSet> cellVisitedSets(nCells);
    List<labelList> cellFrontiers(nCells);

    // Initialize the visited sets and current frontiers
    for (label localCellI = 0; localCellI < nCells; ++localCellI)
    {
        // Get global cell ID
        const label globalCellI = globalCells_.toGlobal(localCellI);

        // Initialize the visited set with the cell itself
        cellVisitedSets[localCellI].insert(globalCellI);

        // Initialize current frontier with immediate neighbors
        DynamicList<label> neighborCellIDs;

        const labelList& cellFaces = cells[localCellI];

        forAll(cellFaces, faceI)
        {
            const label curFaceID = cellFaces[faceI];

            if (mesh.isInternalFace(curFaceID))
            {
                // Internal face
                const label own = owner[curFaceID];
                const label nei = neighbour[curFaceID];
                const label neiLocalCellI = (own == localCellI) ? nei : own;
                const label neiGlobalCellI =
                    globalCells_.toGlobal(neiLocalCellI);

                // Check if neiCellID is already in neighbors
                if (findIndex(neighborCellIDs, neiGlobalCellI) == -1)
                {
                    neighborCellIDs.append(neiGlobalCellI);
                }
            }
            else
            {
                // Boundary face
                const label patchID =
                    mesh.boundaryMesh().whichPatch(curFaceID);
                const polyPatch& pp = mesh.boundaryMesh()[patchID];

                if (isA<processorPolyPatch>(pp))
                {
                    // Lookup the global cell ID across the processor patch
                    const label localFaceI = curFaceID - pp.start();
                    const label neiGlobalCellI =
                        procPatchNeiGlobalCellIDs[patchID][localFaceI];

                    // Check if neiCellID is already in neighbors
                    if (findIndex(neighborCellIDs, neiGlobalCellI) == -1)
                    {
                        neighborCellIDs.append(neiGlobalCellI);
                    }
                }
                // else {} // Physical boundary, no action needed
            }
        }

        // Add immediate neighbors to the visited set and current frontier
        cellVisitedSets[localCellI].insert(neighborCellIDs);
        cellFrontiers[localCellI] = neighborCellIDs;
    }

    // Now perform N-1 iterations (we already have the first layer)
    for (label layer = 2; layer <= nLayers_; ++layer)
    {
        // Prepare the next layer's frontiers
        List<labelList> nextFrontiers(nCells);

        // Maps for inter-processor communication
        // Map from processor ID to Pair(originGlobalCellI, frontGlobalCellI)
        // What we want from other procs?
        //    cell-cells of a given cell on their proc
        // What to send?
        //     originGlobalCellID => cell whose stencil we are making
        //     frontGlobalCellID => whose cell-cells we want (as global IDs)
        //     myProcNo => so the other proc knows who to send the info back to
        //             may not be explicitly needed, as implicitly known
        Map<List<labelPair>> sendMap;

        for (label localCellI = 0; localCellI < nCells; ++localCellI)
        {
            const labelList& currentFrontier = cellFrontiers[localCellI];
            labelHashSet& visitedSet = cellVisitedSets[localCellI];
            labelList& nextFrontier = nextFrontiers[localCellI];

            forAll(currentFrontier, idx)
            {
                const label frontGlobalCellI = currentFrontier[idx];

                if (globalCells_.isLocal(frontGlobalCellI))
                {
                    const label frontLocalCellI =
                        globalCells_.toLocal(frontGlobalCellI);
                    const labelList& frontCellFaces = cells[frontLocalCellI];

                    forAll(frontCellFaces, fI)
                    {
                        const label curFaceID = frontCellFaces[fI];

                        if (mesh.isInternalFace(curFaceID))
                        {
                            // Internal face
                            const label own = owner[curFaceID];
                            const label nei = neighbour[curFaceID];
                            const label neiLocalCellID =
                                (own == frontLocalCellI) ? nei : own;
                            const label neiGlobalCellID =
                                globalCells_.toGlobal(neiLocalCellID);

                            if (!visitedSet.found(neiGlobalCellID))
                            {
                                visitedSet.insert(neiGlobalCellID);
                                nextFrontier.append(neiGlobalCellID);
                            }
                        }
                        else // Boundary face
                        {
                            const label patchID =
                                mesh.boundaryMesh().whichPatch(curFaceID);
                            if (patchID == -1)
                            {
                                FatalErrorInFunction
                                    << "patchID == -1 for face " << curFaceID
                                    << abort(FatalError);
                            }
                            const polyPatch& pp = mesh.boundaryMesh()[patchID];

                            if (isA<processorPolyPatch>(pp))
                            {
                                // Processor boundary face
                                //const processorPolyPatch& procPatch =
                                //    refCast<const processorPolyPatch>(pp);
                                //const label procNo = procPatch.neighbProcNo();

                               // Lookup the global cell ID across the processor patch
                               const label localFaceI = curFaceID - pp.start();
                               const label neiGlobalCellI =
                                   procPatchNeiGlobalCellIDs[patchID][localFaceI];

                               if (!visitedSet.found(neiGlobalCellI))
                               {
                                   visitedSet.insert(neiGlobalCellI);
                                   nextFrontier.append(neiGlobalCellI);
                               }
                            }
                            // else {} // Physical boundary, no action needed
                        }
                    }
                }
                else // frontGlobalCellI is on another proc
                {
                    // We need to request the cell-cells from the processor who
                    // owns this front cell

                    // Determine which processor owns this cell
                    const label procID =
                        globalCells_.whichProcID(frontGlobalCellI);

                    // Origin cell
                    const label globalCellI = globalCells_.toGlobal(localCellI);

                    // Record the origin cell and front cell for communication
                    sendMap(procID).append
                    (
                        labelPair(globalCellI, frontGlobalCellI)
                    );
                }
            }
        }

        // Handle inter-processor communication
        // Prepare data to send to neighboring processors
        Map<List<labelPair>> toSend(Pstream::nProcs());
        Map<List<labelPair>> toReceive(Pstream::nProcs());

        // Create toSend lists
        forAllIter(Map<List<labelPair>>, sendMap, iter)
        {
            const label procNo = iter.key();
            List<labelPair>& sendData = iter();

            toSend(procNo).transfer(sendData);
        }

        // Exchange data with neighboring processors
        Pstream::exchange<List<labelPair>, labelPair>
        (
            toSend, toReceive
        );

        // Clear send map for next communication
        sendMap.clear();

        // Process received data
        forAllConstIter(Map<List<labelPair>>, toReceive, iter)
        {
            const label procI = iter.key();
            const List<labelPair>& receivedData = iter();

            forAll(receivedData, idx)
            {
                const label globalCellI = receivedData[idx].first();
                const label frontGlobalCellI = receivedData[idx].second();

                if (!globalCells_.isLocal(frontGlobalCellI))
                {
                    FatalErrorInFunction
                        << "Global cell " << frontGlobalCellI
                        << " is not on this proc!" << abort(FatalError);
                }

                // Get local ID
                const label frontLocalCellI =
                    globalCells_.toLocal(frontGlobalCellI);

                // Prepare list of cell-cells for the frontGlobalCellI as
                // global IDs

                const labelList& cellFaces = cells[frontLocalCellI];

                forAll(cellFaces, fI)
                {
                    const label curFaceID = cellFaces[fI];

                    if (mesh.isInternalFace(curFaceID))
                    {
                        // Internal face
                        const label own = owner[curFaceID];
                        const label nei = neighbour[curFaceID];
                        const label neiLocalCellI =
                            (own == frontLocalCellI) ? nei : own;
                        const label neiGlobalCellI =
                            globalCells_.toGlobal(neiLocalCellI);

                        sendMap(procI).append
                        (
                            labelPair(globalCellI, neiGlobalCellI)
                        );
                    }
                    else
                    {
                        // Boundary face
                        const label patchID =
                            mesh.boundaryMesh().whichPatch(curFaceID);
                        const polyPatch& pp = mesh.boundaryMesh()[patchID];

                        if (isA<processorPolyPatch>(pp))
                        {
                            // Lookup the global cell ID across the processor patch
                            const label localFaceI = curFaceID - pp.start();
                            const label neiGlobalCellI =
                                procPatchNeiGlobalCellIDs[patchID][localFaceI];

                            sendMap(procI).append
                            (
                                labelPair(globalCellI, neiGlobalCellI)
                            );
                        }
                        // else {} // Physical boundary, no action needed
                    }
                }

                // const label globalCellI = receivedData[idx].first();
            }
        }

        // Handle inter-processor communication

        // Clear communication maps
        toSend.clear();
        toReceive.clear();

        // Populate toSend lists
        forAllIter(Map<List<labelPair>>, sendMap, iter)
        {
            const label procNo = iter.key();
            List<labelPair>& sendData = iter();

            toSend(procNo).transfer(sendData);
        }

        // Exchange data with neighboring processors
        Pstream::exchange<List<labelPair>, labelPair>
        (
            toSend, toReceive
        );

        // Finally, retreive data and add to local visited cells and next
        // frontier

        // Process received data
        // Here, we are receiving the cell-cells from other procs, which we had
        // requested. We will add these to our stencils.
        forAllConstIter(Map<List<labelPair>>, toReceive, iter)
        {
            //const label procI = iter.key();
            const List<labelPair>& receivedData = iter();

            forAll(receivedData, idx)
            {
                // Cell whose stencil we are creating
                const label globalCellI = receivedData[idx].first();

                if (!globalCells_.isLocal(globalCellI))
                {
                    FatalErrorInFunction
                        << "Global cell " << globalCellI
                        << " is not on this proc!" << abort(FatalError);
                }

                // Get local ID
                const label localCellI = globalCells_.toLocal(globalCellI);

                // Local visited cells and next front
                labelHashSet& visitedSet = cellVisitedSets[localCellI];
                labelList& nextFrontier = nextFrontiers[localCellI];

                // Cell to be added to the stencil of globallCellI
                const label neiGlobalCellI = receivedData[idx].second();

                // Add neiGlobalCellI to the visited cells and next front
                if (!visitedSet.found(neiGlobalCellI))
                {
                    visitedSet.insert(neiGlobalCellI);
                    nextFrontier.append(neiGlobalCellI);
                }
            }
        }

        // Update the frontiers for the next iteration
        cellFrontiers = nextFrontiers;
    }

    // At this point, cellVisitedSets contains the N-layer neighborhoods for
    // all cells
    globalCellStencilsPtr_.set(new labelListList(nCells));
    forAll(cellVisitedSets, localCellI)
    {
        globalCellStencilsPtr_()[localCellI] =
            cellVisitedSets[localCellI].toc();
    }

    InfoInFunction
        << "end" << endl;
}


void higherOrderGrad::makeGlobalFaceStencils() const
{
    if (globalFaceStencilsPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set" << exit(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    globalFaceStencilsPtr_.set(new labelListList(mesh.nFaces()));
    labelListList& faceStencils = globalFaceStencilsPtr_();

    // Create the cell stencils
    const labelListList& cellStencils = globalCellStencils();

    // The stencil for each face consists of the cell stencils of the cell(s)
    // containing that face
    // We will only keep unique entries
    forAll(faceStencils, faceI)
    {
        // We will use a set to check avoid adding duplicates
        labelHashSet curStencil;

        // Add owner stencil
        const label ownLocalCellI = faceOwner[faceI];
        curStencil.insert(cellStencils[ownLocalCellI]);

        // Add neighbour stencil
        // Note: we deal with processor patch faces afterwards
        if (mesh.isInternalFace(faceI))
        {
            const label neiLocalCellI = faceNeighbour[faceI];
            curStencil.merge(labelHashSet(cellStencils[neiLocalCellI]));
        }

        // Convert set to labelList
        faceStencils[faceI] = curStencil.toc();
    }

    // If the face is on a processor patch, we need to get the cell
    // stencils from the face cells on the neighbouring processor

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

        // Prepare and send stencils
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Prepare cell stencils for sending
                const labelUList& faceCells = pp.faceCells();
                labelListList patchCellStencils(faceCells.size());
                forAll(patchCellStencils, fI)
                {
                    const label localCellI = faceCells[fI];
                    patchCellStencils[fI] = cellStencils[localCellI];
                }

                // Send to neighbour processor
                UOPstream toNeighbProc(procPatch.neighbProcNo(), pBufs);
                toNeighbProc
                    << patchCellStencils;
            }
        }

        pBufs.finishedSends(); // no-op for blocking

        // Receive data
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                UIPstream fromNeighbProc(procPatch.neighbProcNo(), pBufs);

                labelListList receiveData(pp.size());

                fromNeighbProc
                    >> receiveData;

                // Merge stencils with local face cells
                const labelUList& faceCells = pp.faceCells();
                forAll(faceCells, fI)
                {
                    // Convert data to a set
                    labelHashSet receiveDataSet(receiveData[fI]);

                    // Convert local stencit to a set
                    const label faceID = pp.start() + fI;
                    labelHashSet curStencil(faceStencils[faceID]);

                    // Merge sets
                    curStencil.merge(receiveDataSet);

                    // Update local stencil
                    faceStencils[faceID] = curStencil.toc();
                }
            }
        }
    }
}


void higherOrderGrad::makeStencils() const
{
    //if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    if (Pstream::parRun())
    {
        notImplemented("not implemented for parallel run");
        // makeGlobalCellStencils works in parallel: next build coefficients
        // using these global stencils
    }

    if (stencilsPtr_ || stencilsBoundaryFacesPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const labelListList& cellCells = mesh.cellCells();

    stencilsPtr_.set(new List<DynamicList<label>>(mesh.nCells()));
    List<DynamicList<label>>& stencils = *stencilsPtr_;

    stencilsBoundaryFacesPtr_.set(new List<DynamicList<label>>(mesh.nCells()));
    List<DynamicList<label>>& stencilsBoundaryFaces = *stencilsBoundaryFacesPtr_;

    forAll(stencils, cellI)
    {
       DynamicList<label>& curStencil = stencils[cellI];
       curStencil.setCapacity(maxStencilSize_);
       const labelList& curCellCells = cellCells[cellI];

       labelHashSet stencilCells;
       labelHashSet prevLayer;

       // Add first layer of cells
       forAll(curCellCells, cI)
       {
           stencilCells.insert(curCellCells[cI]);
           prevLayer.insert(curCellCells[cI]);
       }

       // Remaining layers of cells
       for (int layerI = 1; layerI < nLayers_; layerI++)
       {
           labelList prevLayerCells(prevLayer.toc());
           labelHashSet curLayer;

           // Loop over previous layer and add one level of
           // layer neighbours
           for (const label cI : prevLayerCells)
           {
               const labelList& cellINei= mesh.cellCells()[cI];
               forAll (cellINei, nei)
               {
                   if (!stencilCells.found(cellINei[nei]))
                   {
                       curLayer.insert(cellINei[nei]);
                   }
               }
           }

           // Now we have curent layer which we need to add to the stencil
           // and the current layer will now be the previous layer for next
           // loop
           prevLayer.clear();
           prevLayer = curLayer;
           stencilCells.merge(curLayer);
       }

       // Neighbours of first layer will store cellI in stencil so we
       // need to remove it
       stencilCells.erase(cellI);

       curStencil.append(stencilCells.toc());
    }

    forAll(stencils, cellI)
    {
       // Note: stencils are sorted but I do not see problem in that
       stencils[cellI].shrink();
    }


    // Make cell boundary face stencils
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    forAll(boundaryMesh, patchI)
    {
        if
        (
            includePatchInStencils_[patchI]
         && boundaryMesh[patchI].type() != emptyPolyPatch::typeName
         && !boundaryMesh[patchI].coupled()
        )
        {
            const labelUList& faceCells = boundaryMesh[patchI].faceCells();

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                const label gI = faceI + boundaryMesh[patchI].start();
                stencilsBoundaryFaces[cellID].append(gI);
            }
        }
    }

    forAll(stencilsBoundaryFaces, cellI)
    {
        stencilsBoundaryFaces[cellI].shrink();
    }

    //if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }
}


const List<DynamicList<label>>& higherOrderGrad::stencils() const
{
    if (!stencilsPtr_)
    {
        makeStencils();
    }

    return *stencilsPtr_;
}


const labelListList& higherOrderGrad::globalCellStencils() const
{
    if (!globalCellStencilsPtr_)
    {
        makeGlobalCellStencils();
    }

    return globalCellStencilsPtr_();
}


const labelListList& higherOrderGrad::globalFaceStencils() const
{
    if (!globalFaceStencilsPtr_)
    {
        makeGlobalFaceStencils();
    }

    return globalFaceStencilsPtr_();
}


void higherOrderGrad::generateExponents
(
    const label N,
    DynamicList<FixedList<label, 3>>& exponents
) const
{
    if (N_ < 1)
    {
        FatalErrorInFunction
            << "N must be at least 1!" << exit(FatalError);
    }

    // Estimate the number of terms to set the capacity
    const label estimatedSize = (N + 1)*(N + 2)*(N + 3)/6;
    exponents.setCapacity(estimatedSize);

    // Add the constant term first
    exponents.append(FixedList<label, 3>{0, 0, 0});

    for (label n = 1; n <= N; ++n)
    {
        for (label i = n; i >= 0; --i)
        {
            for (label j = n - i; j >= 0; --j)
            {
                label k = n - i - j;
                if (i == 0 && j == 0 && k == 0)
                {
                    // Skip the constant term as it's already added
                    continue;
                }
                FixedList<label, 3> exponent = {i, j, k};
                exponents.append(exponent);
            }
        }
    }

    // Adjust capacity to actual size
    exponents.shrink();
}


void higherOrderGrad::calcQRCoeffs() const
{
    //if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    if (QRInterpCoeffsPtr_ || QRGradCoeffsPtr_)
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    QRInterpCoeffsPtr_.set(new List<DynamicList<scalar>>(mesh.nCells()));
    List<DynamicList<scalar>>& QRInterpCoeffs = *QRInterpCoeffsPtr_;

    QRGradCoeffsPtr_.set(new List<DynamicList<vector>>(mesh.nCells()));
    List<DynamicList<vector>>& QRGradCoeffs = *QRGradCoeffsPtr_;

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();
    const vectorField& CfI = mesh.Cf();

    // Calculate Taylor series exponents
    // 1 for zero order, 4 for 1 order, 10 for second order, etc.
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();
    if (debug)
    {
        Info<< "Np = " << Np << endl;
    }

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    List<DynamicList<scalar>> c(mesh.nCells());
    List<DynamicList<scalar>> cx(mesh.nCells());
    List<DynamicList<scalar>> cy(mesh.nCells());
    List<DynamicList<scalar>> cz(mesh.nCells());

    const List<DynamicList<label>>& stencilsBoundaryFaces =
        this->stencilsBoundaryFaces();
    const List<DynamicList<label>>& stencils = this->stencils();

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiCellID = curStencil[cI];
            const scalar d = mag(CI[neiCellID] - CI[cellI]);
            if (d > maxDist)
            {
                maxDist = d;
            }
        }

        // Loop over neighbours and construct matrix Q
        const label Nn =
            curStencil.size() + stencilsBoundaryFaces[cellI].size();

        // Use matrix format from Eigen/Dense library
        // Avoid initialisation to zero as we will set every entry below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < Nn; ++cI)
        {
            vector dx;
            if (cI < curStencil.size())
            {
                const label neiCellID = curStencil[cI];
                const vector& neiC = CI[neiCellID];
                dx = neiC - CI[cellI];
            }
            else
            {
                const label i = cI - curStencil.size();
                const label globalFaceID = stencilsBoundaryFaces[cellI][i];
                const vector& neiC = CfI[globalFaceID];
                dx = neiC - CI[cellI];
            }

            // Compute monomial values for each exponent
            for (label p = 0; p < Np; ++p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = exponent[2];

               // Compute factorial denominator
               const scalar factorialDenominator =
                   factorials[i]*factorials[j]*factorials[k];

               // Compute and assign monomial value with factorials
               // Note: the order of the quadratic and higher terms may not be
               // the same as the previous manual approach
               Q(p, cI) =
                   pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                  /factorialDenominator;
            }
        }

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);
        //W.setZero(); // no need to waste time initialising

        for (label cI = 0; cI < Nn; cI++)
        {
            scalar d;

            if (cI < curStencil.size())
            {
                const vector& C = CI[cellI];
                const label neiCellID = curStencil[cI];
                const vector& neiC = CI[neiCellID];
                d = mag(neiC - C);
            }
            else
            {
                // For boundary cells we need to add boundary face as
                // neigbour
                const vector& C = CI[cellI];
                const label i = cI - curStencil.size();
                const label globalFaceID = stencilsBoundaryFaces[cellI][i];
                const vector& neiC = CfI[globalFaceID];
                d = mag(neiC - C);
            }

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W.diagonal()[cI] = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> sqrtW =
            W.diagonal().cwiseSqrt().asDiagonal();
        const Eigen::MatrixXd Qhat =
            Q.array().rowwise()*sqrtW.diagonal().transpose().array();

        // B hat
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
            sqrtW.diagonal().asDiagonal();

        // QR decomposition
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

        // Q and R matrices
        const Eigen::MatrixXd O = qr.householderQ();
        const Eigen::MatrixXd& R = qr.matrixQR().triangularView<Eigen::Upper>();

        // Slice Rbar and Qbar, as we do not need full matrix
        // Note: auto is a reference type here (Rbar, Qbar are not copied)
        const auto Rbar = R.topLeftCorner(Np, Np);
        const auto Qbar = O.leftCols(Np);

        // Perform element-wise multiplication and convert to MatrixXd
        const Eigen::MatrixXd QbarBhat =
            (
                Qbar.transpose().array().rowwise()
               *Bhat.diagonal().transpose().array()
            ).matrix();

        // Solve to get A
        // const Eigen::MatrixXd A =
        //     Rbar.colPivHouseholderQr().solve(Qbar.transpose()*Bhat);
        // Solve using the modified QbarBhat
        const Eigen::MatrixXd A = Rbar.colPivHouseholderQr().solve(QbarBhat);

        // To be aware of interpolation accuracy we need to control the
        // condition number
        if (calcConditionNumber_)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd
            (
                Rbar, Eigen::ComputeFullU | Eigen::ComputeFullV
            );
            Eigen::VectorXd singularValues = svd.singularValues();

            conditionNumber()[cellI] =
                singularValues(0)
               /(singularValues(singularValues.size() - 1) + VSMALL);
        }

        c[cellI].setCapacity(A.cols());
        cx[cellI].setCapacity(A.cols());
        cy[cellI].setCapacity(A.cols());
        cz[cellI].setCapacity(A.cols());

        Eigen::RowVectorXd cRow = A.row(0);
        Eigen::RowVectorXd cxRow = A.row(1);
        Eigen::RowVectorXd cyRow = A.row(2);
        Eigen::RowVectorXd czRow = A.row(3);

        for (label i = 0; i < A.cols(); ++i)
        {
            c[cellI].append(cRow(i));
            cx[cellI].append(cxRow(i));
            cy[cellI].append(cyRow(i));
            cz[cellI].append(czRow(i));
        }

        c[cellI].shrink();
        cx[cellI].shrink();
        cy[cellI].shrink();
        cz[cellI].shrink();
    }

    forAll(QRInterpCoeffs, cellI)
    {
       const DynamicList<label>& curStencil = stencils[cellI];
       const label Nn = curStencil.size() + stencilsBoundaryFaces[cellI].size();

       QRInterpCoeffs[cellI].setCapacity(Nn);
       QRGradCoeffs[cellI].setCapacity(Nn);

       for (label I = 0; I < Nn; I++)
       {
           QRInterpCoeffs[cellI].append(c[cellI][I]);
           QRGradCoeffs[cellI].append
           (
               vector(cx[cellI][I], cy[cellI][I], cz[cellI][I])
           );
       }

       QRInterpCoeffs[cellI].shrink();
       QRGradCoeffs[cellI].shrink();
    }

    //if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }
}


void higherOrderGrad::calcGlobalQRCoeffs() const
{
    //if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    if (QRInterpCoeffsPtr_ || QRGradCoeffsPtr_)
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    QRInterpCoeffsPtr_.set(new List<DynamicList<scalar>>(mesh.nCells()));
    List<DynamicList<scalar>>& QRInterpCoeffs = *QRInterpCoeffsPtr_;

    QRGradCoeffsPtr_.set(new List<DynamicList<vector>>(mesh.nCells()));
    List<DynamicList<vector>>& QRGradCoeffs = *QRGradCoeffsPtr_;

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();

    // Collect CI for off-processor cells in the stencils
    Map<vector> globalCI;
    requestGlobalStencilData(CI, globalCI);

    // Calculate Taylor series exponents
    // 1 for zero order, 4 for 1 order, 10 for second order, etc.
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();
    if (debug)
    {
        Info<< "Np = " << Np << endl;
    }

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    List<DynamicList<scalar>> c(mesh.nCells());
    List<DynamicList<scalar>> cx(mesh.nCells());
    List<DynamicList<scalar>> cy(mesh.nCells());
    List<DynamicList<scalar>> cz(mesh.nCells());

    const List<labelList>& stencils = globalCellStencils();

    forAll(stencils, localCellI)
    {
        const labelList& curStencil = stencils[localCellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiGlobalCellID = curStencil[cI];
            scalar d;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);

                d = mag(CI[neiLocalCellID] - CI[localCellI]);
            }
            else
            {
                d = mag(globalCI[neiGlobalCellID] - CI[localCellI]);
            }

            maxDist = max(maxDist, d);
        }

        // Loop over neighbours and construct matrix Q
        const label Nn = curStencil.size();

        // Use matrix format from Eigen/Dense library
        // Avoid initialisation to zero as we will set every entry below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < Nn; ++cI)
        {
            const label neiGlobalCellID = curStencil[cI];
            vector dx;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);

                dx = CI[neiLocalCellID] - CI[localCellI];
            }
            else
            {
                dx = globalCI[neiGlobalCellID] - CI[localCellI];
            }

            // Compute monomial values for each exponent
            for (label p = 0; p < Np; ++p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = exponent[2];

               // Compute factorial denominator
               const scalar factorialDenominator =
                   factorials[i]*factorials[j]*factorials[k];

               // Compute and assign monomial value with factorials
               // Note: the order of the quadratic and higher terms may not be
               // the same as the previous manual approach
               Q(p, cI) =
                   pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                  /factorialDenominator;
            }
        }

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);
        //W.setZero(); // no need to waste time initialising

        for (label cI = 0; cI < Nn; cI++)
        {
            const label neiGlobalCellID = curStencil[cI];
            scalar d;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);

                d = mag(CI[neiLocalCellID] - CI[localCellI]);
            }
            else
            {
                d = mag(globalCI[neiGlobalCellID] - CI[localCellI]);
            }

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W.diagonal()[cI] = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> sqrtW =
            W.diagonal().cwiseSqrt().asDiagonal();
        const Eigen::MatrixXd Qhat =
            Q.array().rowwise()*sqrtW.diagonal().transpose().array();

        // B hat
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
            sqrtW.diagonal().asDiagonal();

        // QR decomposition
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

        // Q and R matrices
        const Eigen::MatrixXd O = qr.householderQ();
        const Eigen::MatrixXd& R = qr.matrixQR().triangularView<Eigen::Upper>();

        // Slice Rbar and Qbar, as we do not need full matrix
        // Note: auto is a reference type here (Rbar, Qbar are not copied)
        const auto Rbar = R.topLeftCorner(Np, Np);
        const auto Qbar = O.leftCols(Np);

        // Perform element-wise multiplication and convert to MatrixXd
        const Eigen::MatrixXd QbarBhat =
            (
                Qbar.transpose().array().rowwise()
               *Bhat.diagonal().transpose().array()
            ).matrix();

        // Solve to get A
        // const Eigen::MatrixXd A =
        //     Rbar.colPivHouseholderQr().solve(Qbar.transpose()*Bhat);
        // Solve using the modified QbarBhat
        const Eigen::MatrixXd A = Rbar.colPivHouseholderQr().solve(QbarBhat);

        // To be aware of interpolation accuracy we need to control the
        // condition number
        if (calcConditionNumber_)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd
            (
                Rbar, Eigen::ComputeFullU | Eigen::ComputeFullV
            );
            Eigen::VectorXd singularValues = svd.singularValues();

            conditionNumber()[localCellI] =
                singularValues(0)
               /(singularValues(singularValues.size() - 1) + VSMALL);
        }

        c[localCellI].setCapacity(A.cols());
        cx[localCellI].setCapacity(A.cols());
        cy[localCellI].setCapacity(A.cols());
        cz[localCellI].setCapacity(A.cols());

        Eigen::RowVectorXd cRow = A.row(0);
        Eigen::RowVectorXd cxRow = A.row(1);
        Eigen::RowVectorXd cyRow = A.row(2);
        Eigen::RowVectorXd czRow = A.row(3);

        for (label i = 0; i < A.cols(); ++i)
        {
            c[localCellI].append(cRow(i));
            cx[localCellI].append(cxRow(i));
            cy[localCellI].append(cyRow(i));
            cz[localCellI].append(czRow(i));
        }

        c[localCellI].shrink();
        cx[localCellI].shrink();
        cy[localCellI].shrink();
        cz[localCellI].shrink();
    }

    forAll(QRInterpCoeffs, localCellI)
    {
       const labelList& curStencil = stencils[localCellI];
       const label Nn = curStencil.size();

       QRInterpCoeffs[localCellI].setCapacity(Nn);
       QRGradCoeffs[localCellI].setCapacity(Nn);

       for (label I = 0; I < Nn; I++)
       {
           QRInterpCoeffs[localCellI].append(c[localCellI][I]);
           QRGradCoeffs[localCellI].append
           (
               vector(cx[localCellI][I], cy[localCellI][I], cz[localCellI][I])
           );
       }

       QRInterpCoeffs[localCellI].shrink();
       QRGradCoeffs[localCellI].shrink();
    }

    //if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }
}


void higherOrderGrad::calcGlobalQRFaceCoeffs() const
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    if (QRGradFaceCoeffsPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    QRGradFaceCoeffsPtr_.set(new List<DynamicList<vector>>(mesh.nFaces()));
    List<DynamicList<vector>>& QRGradCoeffs = *QRGradFaceCoeffsPtr_;

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();
    const vectorField& CfI = mesh.Cf();

    // Collect CI for off-processor cells in the stencils
    Map<vector> globalCI;
    requestGlobalStencilData(CI, globalCI);

    // Calculate Taylor series exponents
    // 1 for zero order, 4 for 1 order, 10 for second order, etc.
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();
    if (debug)
    {
        Info<< "Np = " << Np << endl;
    }

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    // Coefficients
    List<DynamicList<scalar>> c(mesh.nFaces());
    List<DynamicList<scalar>> cx(mesh.nFaces());
    List<DynamicList<scalar>> cy(mesh.nFaces());
    List<DynamicList<scalar>> cz(mesh.nFaces());

    const List<labelList>& stencils = globalFaceStencils();

    forAll(stencils, faceI)
    {
        const labelList& curStencil = stencils[faceI];

        // Centre of current face
        const vector& curCf = CfI[faceI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiGlobalCellID = curStencil[cI];

            scalar d;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);
                d = mag(CI[neiLocalCellID] - curCf);
            }
            else
            {
                d = mag(globalCI[neiGlobalCellID] - curCf);
            }

            maxDist = max(maxDist, d);
        }

        // Loop over neighbours and construct matrix Q
        const label Nn = curStencil.size();

        // Use matrix format from Eigen/Dense library
        // Avoid initialisation to zero as we will set every entry below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < Nn; ++cI)
        {
            const label neiGlobalCellID = curStencil[cI];
            vector dx;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);

                dx = CI[neiLocalCellID] - curCf;
            }
            else
            {
                dx = globalCI[neiGlobalCellID] - curCf;
            }

            // Compute monomial values for each exponent
            for (label p = 0; p < Np; ++p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = exponent[2];

               // Compute factorial denominator
               const scalar factorialDenominator =
                   factorials[i]*factorials[j]*factorials[k];

               // Compute and assign monomial value with factorials
               // Note: the order of the quadratic and higher terms may not be
               // the same as the previous manual approach
               Q(p, cI) =
                   pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                  /factorialDenominator;
            }
        }

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);
        //W.setZero(); // no need to waste time initialising

        for (label cI = 0; cI < Nn; cI++)
        {
            const label neiGlobalCellID = curStencil[cI];
            scalar d;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);

                d = mag(CI[neiLocalCellID] - curCf);
            }
            else
            {
                d = mag(globalCI[neiGlobalCellID] - curCf);
            }

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W.diagonal()[cI] = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic> sqrtW =
            W.diagonal().cwiseSqrt().asDiagonal();
        const Eigen::MatrixXd Qhat =
            Q.array().rowwise()*sqrtW.diagonal().transpose().array();

        // B hat
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
            sqrtW.diagonal().asDiagonal();

        // QR decomposition
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

        // Q and R matrices
        const Eigen::MatrixXd O = qr.householderQ();
        const Eigen::MatrixXd& R = qr.matrixQR().triangularView<Eigen::Upper>();

        // Slice Rbar and Qbar, as we do not need full matrix
        // Note: auto is a reference type here (Rbar, Qbar are not copied)
        const auto Rbar = R.topLeftCorner(Np, Np);
        const auto Qbar = O.leftCols(Np);

        // Perform element-wise multiplication and convert to MatrixXd
        const Eigen::MatrixXd QbarBhat =
            (
                Qbar.transpose().array().rowwise()
               *Bhat.diagonal().transpose().array()
            ).matrix();

        // Solve to get A
        // const Eigen::MatrixXd A =
        //     Rbar.colPivHouseholderQr().solve(Qbar.transpose()*Bhat);
        // Solve using the modified QbarBhat
        const Eigen::MatrixXd A = Rbar.colPivHouseholderQr().solve(QbarBhat);

        // TODO: how best to display a surface field? Maybe use fvc::average to
        // create a vol field
        // // To be aware of interpolation accuracy we need to control the
        // // condition number
        // if (calcConditionNumber_)
        // {
        //     Eigen::JacobiSVD<Eigen::MatrixXd> svd
        //     (
        //         Rbar, Eigen::ComputeFullU | Eigen::ComputeFullV
        //     );
        //     Eigen::VectorXd singularValues = svd.singularValues();

        //     conditionNumber()[localCellI] =
        //         singularValues(0)
        //        /(singularValues(singularValues.size() - 1) + VSMALL);
        // }

        c[faceI].setCapacity(A.cols());
        cx[faceI].setCapacity(A.cols());
        cy[faceI].setCapacity(A.cols());
        cz[faceI].setCapacity(A.cols());

        Eigen::RowVectorXd cRow = A.row(0);
        Eigen::RowVectorXd cxRow = A.row(1);
        Eigen::RowVectorXd cyRow = A.row(2);
        Eigen::RowVectorXd czRow = A.row(3);

        for (label i = 0; i < A.cols(); ++i)
        {
            c[faceI].append(cRow(i));
            cx[faceI].append(cxRow(i));
            cy[faceI].append(cyRow(i));
            cz[faceI].append(czRow(i));
        }

        c[faceI].shrink();
        cx[faceI].shrink();
        cy[faceI].shrink();
        cz[faceI].shrink();
    }

    forAll(QRGradCoeffs, faceI)
    {
       const labelList& curStencil = stencils[faceI];
       const label Nn = curStencil.size();

       QRGradCoeffs[faceI].setCapacity(Nn);

       for (label I = 0; I < Nn; I++)
       {
           QRGradCoeffs[faceI].append
           (
               vector(cx[faceI][I], cy[faceI][I], cz[faceI][I])
           );
       }

       QRGradCoeffs[faceI].shrink();
    }

    if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }
}


void higherOrderGrad::calcCholeskyCoeffs() const
{
    Info<< "higherOrderGrad::calcCholeskyCoeffs()" << endl;

    if (choleskyPtr_ || QhatPtr_ || sqrtWPtr_)
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    choleskyPtr_.set(new List<Eigen::LLT<Eigen::MatrixXd>>(mesh.nCells()));
    List<Eigen::LLT<Eigen::MatrixXd>>& cholesky = choleskyPtr_();

    QhatPtr_.set(new List<Eigen::MatrixXd>(mesh.nCells()));
    List<Eigen::MatrixXd>& Qhat = QhatPtr_();

    sqrtWPtr_.set
    (
        new List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>(mesh.nCells())
    );
    List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>& sqrtW = sqrtWPtr_();

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();
    const vectorField& CfI = mesh.Cf();

    // Calculate Taylor series exponents
    // 1 for zero order, 4 for 1 order, 10 for second order, etc.
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();
    if (debug)
    {
        Info<< "Np = " << Np << endl;
    }

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    List<DynamicList<scalar>> c(mesh.nCells());
    List<DynamicList<scalar>> cx(mesh.nCells());
    List<DynamicList<scalar>> cy(mesh.nCells());
    List<DynamicList<scalar>> cz(mesh.nCells());

    const List<DynamicList<label>>& stencilsBoundaryFaces =
        this->stencilsBoundaryFaces();
    const List<DynamicList<label>>& stencils = this->stencils();

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiCellID = curStencil[cI];
            const scalar d = mag(CI[neiCellID] - CI[cellI]);
            if (d > maxDist)
            {
                maxDist = d;
            }
        }

        // Loop over neighbours and construct matrix Q
        const label Nn = curStencil.size() + stencilsBoundaryFaces[cellI].size();

        // Use matrix format from Eigen/Dense library
        // Avoid initialisation to zero as we will set every entry below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < Nn; ++cI)
        {
            vector dx;
            if (cI < curStencil.size())
            {
                const label neiCellID = curStencil[cI];
                const vector& neiC = CI[neiCellID];
                dx = neiC - CI[cellI];
            }
            else
            {
                const label i = cI - curStencil.size();
                const label globalFaceID = stencilsBoundaryFaces[cellI][i];
                const vector& neiC = CfI[globalFaceID];
                dx = neiC - CI[cellI];
            }

            // Compute monomial values for each exponent
            for (label p = 0; p < Np; ++p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = exponent[2];

               // Compute factorial denominator
               const scalar factorialDenominator =
                   factorials[i]*factorials[j]*factorials[k];

               // Compute and assign monomial value with factorials
               // Note: the order of the quadratic and higher terms may not be
               // the same as the previous manual approach
               Q(p, cI) =
                   pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                  /factorialDenominator;
            }
        }

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);
        //W.setZero(); // no need to waste time initialising

        for (label cI = 0; cI < Nn; cI++)
        {
            scalar d;

            if (cI < curStencil.size())
            {
                const vector& C = CI[cellI];
                const label neiCellID = curStencil[cI];
                const vector& neiC = CI[neiCellID];
                d = mag(neiC - C);
            }
            else
            {
                // For boundary cells we need to add boundary face as
                // neigbour
                const vector& C = CI[cellI];
                const label i = cI - curStencil.size();
                const label globalFaceID = stencilsBoundaryFaces[cellI][i];
                const vector& neiC = CfI[globalFaceID];
                d = mag(neiC - C);
            }

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W.diagonal()[cI] = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        sqrtW[cellI] = W.diagonal().cwiseSqrt().asDiagonal();

        // B hat
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
            sqrtW[cellI].diagonal().asDiagonal();

        // Cholesky decomposition of the "normal equations"

        // Transpose Q to follow the standard convention
        // TODO: avoid this by assigning Q correctly from the start!
        // It may be clear to seperate QR and Cholesky into their own
        // functions
        Q = Q.transpose().eval();

        // Compute Q_hat = Q * W^{1/2}
        Qhat[cellI] = Q.array().colwise()*sqrtW[cellI].diagonal().array();

        // Compute N = Q_hat^T * Q_hat = Q^T W Q
        const Eigen::MatrixXd N = Qhat[cellI].transpose()*Qhat[cellI];

        if (debug)
        {
            Eigen::FullPivLU<Eigen::MatrixXd> lu(Q);
            int rank = lu.rank();
            Info<< "Rank of Q: " << rank << nl
                << "Q rows: " << Q.rows() << nl
                << "Q cols: " << Q.cols() << endl;

            if (rank < Q.cols())
            {
                WarningInFunction
                    << "Design matrix Q is rank-deficient!" << endl;
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(N);

            if (eigensolver.info() != Eigen::Success)
            {
                WarningInFunction
                    << "Eigenvalue computation failed!" << endl;

                Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();

                std::cout
                    << "Eigenvalues of N: " << eigenvalues.transpose()
                    << std::endl;
            }
        }

        // Perform Cholesky decomposition
        cholesky[cellI].compute(N);

        if (cholesky[cellI].info() != Eigen::Success)
        {
            FatalErrorInFunction
                << "Cholesky decomposition failed; "
                << "matrix is not positive definite."
                << exit(FatalError);
        }
    }

    Info<< "higherOrderGrad::calcCholeskyCoeffs(): end" << endl;
}


void higherOrderGrad::calcGlobalCholeskyCoeffs() const
{
    InfoInFunction
        << "start" << endl;

    if (choleskyPtr_ || QhatPtr_ || sqrtWPtr_)
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    choleskyPtr_.set(new List<Eigen::LLT<Eigen::MatrixXd>>(mesh.nCells()));
    List<Eigen::LLT<Eigen::MatrixXd>>& cholesky = choleskyPtr_();

    QhatPtr_.set(new List<Eigen::MatrixXd>(mesh.nCells()));
    List<Eigen::MatrixXd>& Qhat = QhatPtr_();

    sqrtWPtr_.set
    (
        new List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>(mesh.nCells())
    );
    List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>& sqrtW = sqrtWPtr_();

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();

    // Collect CI for off-processor cells in the stencils
    Map<vector> globalCI;
    requestGlobalStencilData(CI, globalCI);

    // Calculate Taylor series exponents
    // 1 for zero order, 4 for 1 order, 10 for second order, etc.
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();
    if (debug)
    {
        Info<< "Np = " << Np << endl;
    }

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    List<DynamicList<scalar>> c(mesh.nCells());
    List<DynamicList<scalar>> cx(mesh.nCells());
    List<DynamicList<scalar>> cy(mesh.nCells());
    List<DynamicList<scalar>> cz(mesh.nCells());

    const labelListList& stencils = this->globalCellStencils();

    forAll(stencils, localCellI)
    {
        const labelList& curStencil = stencils[localCellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiGlobalCellID = curStencil[cI];

            scalar d;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);
                d = mag(CI[neiLocalCellID] - CI[localCellI]);
            }
            else
            {
                d = mag(globalCI[neiGlobalCellID] - CI[localCellI]);
            }

            maxDist = max(maxDist, d);
        }

        // Loop over neighbours and construct matrix Q
        const label Nn = curStencil.size();

        // Use matrix format from Eigen/Dense library
        // Avoid initialisation to zero as we will set every entry below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < Nn; ++cI)
        {
            const label neiGlobalCellID = curStencil[cI];
            vector dx;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);
                dx = CI[neiLocalCellID] - CI[localCellI];
            }
            else
            {
                dx = globalCI[neiGlobalCellID] - CI[localCellI];
            }

            // Compute monomial values for each exponent
            for (label p = 0; p < Np; ++p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = exponent[2];

               // Compute factorial denominator
               const scalar factorialDenominator =
                   factorials[i]*factorials[j]*factorials[k];

               // Compute and assign monomial value with factorials
               // Note: the order of the quadratic and higher terms may not be
               // the same as the previous manual approach
               Q(p, cI) =
                   pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                  /factorialDenominator;
            }
        }

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);
        //W.setZero(); // no need to waste time initialising

        for (label cI = 0; cI < Nn; cI++)
        {
            const label neiGlobalCellID = curStencil[cI];
            scalar d;
            if (globalCells_.isLocal(neiGlobalCellID))
            {
                const label neiLocalCellID =
                    globalCells_.toLocal(neiGlobalCellID);

                d = mag(CI[neiLocalCellID] - CI[localCellI]);
            }
            else
            {
                d = mag(globalCI[neiGlobalCellID] - CI[localCellI]);
            }

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W.diagonal()[cI] = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        sqrtW[localCellI] = W.diagonal().cwiseSqrt().asDiagonal();

        // B hat
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
            sqrtW[localCellI].diagonal().asDiagonal();

        // Cholesky decomposition of the "normal equations"

        // Transpose Q to follow the standard convention
        // TODO: avoid this by assigning Q correctly from the start!
        Q = Q.transpose().eval();

        // Compute Q_hat = Q * W^{1/2}
        Qhat[localCellI] =
            Q.array().colwise()*sqrtW[localCellI].diagonal().array();

        // Compute N = Q_hat^T * Q_hat = Q^T W Q
        const Eigen::MatrixXd N = Qhat[localCellI].transpose()*Qhat[localCellI];

        if (debug)
        {
            Eigen::FullPivLU<Eigen::MatrixXd> lu(Q);
            int rank = lu.rank();
            Info<< "Rank of Q: " << rank << nl
                << "Q rows: " << Q.rows() << nl
                << "Q cols: " << Q.cols() << endl;

            if (rank < Q.cols())
            {
                WarningInFunction
                    << "Design matrix Q is rank-deficient!" << endl;
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(N);

            if (eigensolver.info() != Eigen::Success)
            {
                WarningInFunction
                    << "Eigenvalue computation failed!" << endl;

                Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();

                std::cout
                    << "Eigenvalues of N: " << eigenvalues.transpose()
                    << std::endl;
            }
        }

        // Perform Cholesky decomposition
        cholesky[localCellI].compute(N);

        if (cholesky[localCellI].info() != Eigen::Success)
        {
            FatalErrorInFunction
                << "Cholesky decomposition failed; "
                << "matrix is not positive definite."
                << exit(FatalError);
        }
    }

    InfoInFunction
        << "end" << endl;
}


const List<DynamicList<scalar>>& higherOrderGrad::QRInterpCoeffs() const
{
    if (!QRInterpCoeffsPtr_)
    {
        if (useGlobalStencils_)
        {
            calcGlobalQRCoeffs();
        }
        else
        {
            calcQRCoeffs();
        }
    }

    return QRInterpCoeffsPtr_();
}


const List<DynamicList<vector>>& higherOrderGrad::QRGradCoeffs() const
{
    if (!QRGradCoeffsPtr_)
    {
        if (useGlobalStencils_)
        {
            calcGlobalQRCoeffs();
        }
        else
        {
            calcQRCoeffs();
        }
    }

    return QRGradCoeffsPtr_();
}


const List<DynamicList<vector>>& higherOrderGrad::QRGradFaceCoeffs() const
{
    if (!QRGradFaceCoeffsPtr_)
    {
        calcGlobalQRFaceCoeffs();
    }

    return QRGradFaceCoeffsPtr_();
}


const List<DynamicList<label>>& higherOrderGrad::stencilsBoundaryFaces() const
{
    if (!stencilsBoundaryFacesPtr_)
    {
        makeStencils();
    }

    return stencilsBoundaryFacesPtr_();
}


const List<Eigen::LLT<Eigen::MatrixXd>>& higherOrderGrad::cholesky() const
{
    if (!choleskyPtr_)
    {
        if (useGlobalStencils_)
        {
            calcGlobalCholeskyCoeffs();
        }
        else
        {
            calcCholeskyCoeffs();
        }
    }

    return choleskyPtr_();
}


const List<Eigen::MatrixXd>& higherOrderGrad::Qhat() const
{
    if (!QhatPtr_)
    {
        if (useGlobalStencils_)
        {
            calcGlobalCholeskyCoeffs();
        }
        else
        {
            calcCholeskyCoeffs();
        }
    }

    return QhatPtr_();
}


const List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>&
higherOrderGrad::sqrtW() const
{
    if (!sqrtWPtr_)
    {
        if (useGlobalStencils_)
        {
            calcGlobalCholeskyCoeffs();
        }
        else
        {
            calcCholeskyCoeffs();
        }
    }

    return sqrtWPtr_();
}


volScalarField& higherOrderGrad::conditionNumber() const
{
    if (!conditionNumberPtr_)
    {
        makeConditionNumber();
    }

    return conditionNumberPtr_();
}


void higherOrderGrad::makeConditionNumber() const
{
    if (conditionNumberPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set" << exit(FatalError);
    }

    conditionNumberPtr_.set
    (
        new volScalarField
        (
           IOobject
           (
               "conditionNumber",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh_,
           dimensionedScalar("0", dimless, Zero)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

higherOrderGrad::higherOrderGrad
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
:
    mesh_(mesh),
    includePatchInStencils_(includePatchInStencils),
    N_(readInt(dict.lookup("N"))),
    nLayers_(readInt(dict.lookup("nLayers"))),
    k_(readScalar(dict.lookup("k"))),
    maxStencilSize_(readInt(dict.lookup("maxStencilSize"))),
    globalCells_(mesh.nCells()),
    useQRDecomposition_(dict.lookup("useQRDecomposition")),
    useGlobalStencils_(dict.lookup("useGlobalStencils")),
    calcConditionNumber_(dict.lookup("calcConditionNumber")),
    conditionNumberPtr_(),
    stencilsPtr_(),
    stencilsBoundaryFacesPtr_(),
    globalCellStencilsPtr_(),
    QRInterpCoeffsPtr_(),
    QRGradCoeffsPtr_(),
    QRGradFaceCoeffsPtr_(),
    choleskyPtr_(),
    QhatPtr_(),
    sqrtWPtr_()
{
    if (calcConditionNumber_)
    {
        if (!useQRDecomposition_)
        {
            FatalErrorInFunction
                << "useQRDecomposition must be 'on' when "
                << "`calcConditionNumber` is 'on'" << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

higherOrderGrad::~higherOrderGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volTensorField> higherOrderGrad::grad(const volVectorField& D) const
{
    auto test =  fGradGaussPoints(D);

    if (useQRDecomposition_)
    {
        if (useGlobalStencils_)
        {
            return gradGlobalQR(D);
        }
        else
        {
            return gradQR(D);
        }
    }
    else
    {
        if (useGlobalStencils_)
        {
            return gradGlobalCholesky(D);
        }
        else
        {
            return gradCholesky(D);
        }
    }
}


tmp<surfaceTensorField> higherOrderGrad::fGrad(const volVectorField& D) const
{
    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<surfaceTensorField> tgradD
    (
        new surfaceTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero)
        )
    );
    surfaceTensorField& gradD = tgradD.ref();
    tensorField& gradDI = gradD;

    const List<labelList>& stencils = globalFaceStencils();
    const List<DynamicList<vector>>& QRGradCoeffs = QRGradFaceCoeffs();
    const vectorField& DI = D;

    // Collect DI for off-processor cells in the stencils
    Map<vector> globalDI;
    requestGlobalStencilData(DI, globalDI);

    forAll(stencils, faceI)
    {
        const labelList& curStencil = stencils[faceI];
        const label Nn = curStencil.size();

        // Loop over stencil and multiply stencil cell values with
        // corresponding interpolation coefficient
        for (label cI = 0; cI < Nn; cI++)
        {
            const label neiGlobalCellI = curStencil[cI];

            if (globalCells_.isLocal(neiGlobalCellI))
            {
                const label neiLocalCellI =
                    globalCells_.toLocal(neiGlobalCellI);

                if (mesh.isInternalFace(faceI))
                {
                    gradDI[faceI] +=
                        QRGradCoeffs[faceI][cI]*DI[neiLocalCellI];
                }
                else
                {
                    // Boundary face
                    const label patchID =
                        mesh.boundaryMesh().whichPatch(faceI);
                    const polyPatch& pp = mesh.boundaryMesh()[patchID];

                    // Note: there is no special treatment for processor patches
                    // as all patches may use global data depending on their
                    // stencils
                    if (!isA<emptyPolyPatch>(pp))
                    {
                        const label localFaceI = faceI - pp.start();
                        gradD.boundaryFieldRef()[patchID][localFaceI] +=
                            QRGradCoeffs[faceI][cI]*DI[neiLocalCellI];
                    }
                }
            }
            else // global cell in the stencil
            {
                if (mesh.isInternalFace(faceI))
                {
                    gradD[faceI] +=
                        QRGradCoeffs[faceI][cI]*globalDI[neiGlobalCellI];
                }
                else
                {
                    // Boundary face
                    const label patchID =
                        mesh.boundaryMesh().whichPatch(faceI);
                    const polyPatch& pp = mesh.boundaryMesh()[patchID];

                    // Note: there is no special treatment for processor patches
                    // as all patches may use global data depending on their
                    // stencils
                    if (!isA<emptyPolyPatch>(pp))
                    {
                        const label localFaceI = faceI - pp.start();
                        gradD.boundaryFieldRef()[patchID][localFaceI] +=
                            QRGradCoeffs[faceI][cI]*globalDI[neiGlobalCellI];
                    }
                }
            }
        }
    }

    gradD.correctBoundaryConditions();

    return tgradD;
}


refPtr<List<List<tensor>>> higherOrderGrad::fGradGaussPoints(const volVectorField& D) const
{
    // Number of quadrature points per triangle.
    const label triQuadraturePtsNb = 3;

    const fvMesh& mesh = mesh_;

    const pointField& pts = mesh.points();

    // Prepare the return field
    // Philip: I'm using refPtr, tmp is more appropriate but I had a compilation
    // problems, when it is used for List<List<tensor>
    refPtr<List<List<tensor>>> tgradDGP(new List<List<tensor>>(mesh.nFaces()));
    List<List<tensor>>& gradDGP = tgradDGP.ref();

    forAll(gradDGP, i)
    {
        List<tensor>& faceGradGP = gradDGP[i];

        const label nbOfTriangles = mesh.faces()[i].size();
        const label nbOfGaussPoints = nbOfTriangles * triQuadraturePtsNb;

        faceGradGP.setSize(nbOfGaussPoints);
    }

    // 1. Stage - decompose faces into triangles. Store triangle points

    // Triangulate each face and store points of each triangle
    List<List<triPoints>> faceTri(mesh.nFaces());
    forAll(faceTri, i)
    {
        List<triPoints>& fT = faceTri[i];

        fT.setSize(mesh.faces()[i].size());
    }

    // Loop over faces and decompose each face, store triangles of each face
    for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
    {
        const face& f = mesh.faces()[faceI];

        const point fc = f.centre(pts);

        const label nPoints = f.size();

        label nextpI;
        for (label pI = 0; pI<nPoints; ++pI)
        {
            if (pI < f.size() - 1)
            {
                nextpI = pI + 1;
            }
            else
            {
                nextpI = 0;
            }

            const triPoints tri
            (
                pts[f[pI]],
                pts[f[nextpI]],
                fc
            );

            faceTri[faceI][pI] = tri;
        }
    }

    // 2. Stage - for each triangle calculate Gauss point locations and store
    //            corresponding weights

    // Gauss point locations on each face
    List<List<point>> faceGP(mesh.nFaces());
    forAll(faceGP, i)
    {
        List<point>& fGP = faceGP[i];
        fGP.setSize(mesh.faces()[i].size() * triQuadraturePtsNb);
    }

    // Gauss point weights
    List<List<scalar>> faceGPW(mesh.nFaces());
    forAll(faceGPW, i)
    {
        List<scalar>& fGPW = faceGPW[i];
        fGPW.setSize(mesh.faces()[i].size() * triQuadraturePtsNb);
    }

    // Loop over faces; loop over triangles. Store Gauss point locations and
    // weight
    forAll(faceTri, faceI)
    {
        List<triPoints>& fT = faceTri[faceI];

        forAll(fT, tI)
        {
            const triPoints& tp = fT[tI];

            // Sad bi tu isla funkcija koja mi vraca Gaussove tocke i tezine
        }
    }

    return tgradDGP;
}


tmp<volTensorField> higherOrderGrad::gradQR(const volVectorField& D) const
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<volTensorField> tgradD
    (
        new volTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero),
           "zeroGradient"
        )
    );
    volTensorField& gradD = tgradD.ref();

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const List<DynamicList<label>>& stencils = this->stencils();
    const List<DynamicList<label>>& stencilsBoundaryFaces =
        this->stencilsBoundaryFaces();
    const List<DynamicList<vector>>& QRGradCoeffs = this->QRGradCoeffs();

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];
        const label Nn = curStencil.size() + stencilsBoundaryFaces[cellI].size();

        // Loop over stencil and multiply stencil cell values with
        // corresponding interpolation coefficient
        for (label cI = 0; cI < Nn; cI++)
        {
            if (cI < curStencil.size())
            {
                gradD[cellI] += QRGradCoeffs[cellI][cI]*D[curStencil[cI]];
            }
            else
            {
                const label i = cI - curStencil.size();
                const label globalFaceID = stencilsBoundaryFaces[cellI][i];

                vector boundaryD = vector::zero;

                forAll(boundaryMesh, patchI)
                {
                    if
                    (
                        includePatchInStencils_[patchI]
                     && boundaryMesh[patchI].type() != emptyPolyPatch::typeName
                     && !boundaryMesh[patchI].coupled()
                    )
                    {
                        const label start = boundaryMesh[patchI].start();
                        const label nFaces = boundaryMesh[patchI].nFaces();

                        if (globalFaceID >= start && globalFaceID < start + nFaces)
                        {
                            const label k = globalFaceID - start;
                            boundaryD = D.boundaryField()[patchI][k];
                        }
                    }
                }

                gradD[cellI] += QRGradCoeffs[cellI][cI]*boundaryD;
            }
        }
    }

    gradD.correctBoundaryConditions();

    if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }

    return tgradD;
}


tmp<volTensorField> higherOrderGrad::gradGlobalQR(const volVectorField& D) const
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<volTensorField> tgradD
    (
        new volTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero),
           "zeroGradient"
        )
    );
    volTensorField& gradD = tgradD.ref();

    const List<labelList>& stencils = globalCellStencils();
    const List<DynamicList<vector>>& QRGradCoeffs = this->QRGradCoeffs();
    const vectorField& DI = D;

    // Collect DI for off-processor cells in the stencils
    Map<vector> globalDI;
    requestGlobalStencilData(DI, globalDI);

    forAll(stencils, localCellI)
    {
        const labelList& curStencil = stencils[localCellI];
        const label Nn = curStencil.size();

        // Loop over stencil and multiply stencil cell values with
        // corresponding interpolation coefficient
        for (label cI = 0; cI < Nn; cI++)
        {
            const label neiGlobalCellI = curStencil[cI];

            if (globalCells_.isLocal(neiGlobalCellI))
            {
                const label neiLocalCellI =
                    globalCells_.toLocal(neiGlobalCellI);

                gradD[localCellI] +=
                    QRGradCoeffs[localCellI][cI]*DI[neiLocalCellI];
            }
            else
            {
                gradD[localCellI] +=
                    QRGradCoeffs[localCellI][cI]*globalDI[neiGlobalCellI];
            }
        }
    }

    gradD.correctBoundaryConditions();

    if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }

    return tgradD;
}


tmp<volTensorField> higherOrderGrad::gradCholesky(const volVectorField& D) const
{
    if (debug)
    {
        Info<< "higherOrderGrad::gradCholesky(...)" << endl;
    }

    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<volTensorField> tgradD
    (
        new volTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero),
           "zeroGradient"
        )
    );
    volTensorField& gradD = tgradD.ref();

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const List<DynamicList<label>>& stencils = this->stencils();
    const List<DynamicList<label>>& stencilsBoundaryFaces =
        this->stencilsBoundaryFaces();
    const List<Eigen::LLT<Eigen::MatrixXd>>& cholesky = this->cholesky();
    const List<Eigen::MatrixXd>& Qhat = this->Qhat();
    const List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>& sqrtW =
        this->sqrtW();

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];
        const label Nn = curStencil.size() + stencilsBoundaryFaces[cellI].size();

        for (label cmptI = 0; cmptI < vector::nComponents; ++cmptI)
        {
            // Prepare right-hand side vector (y) for Cholesky decomposition
            Eigen::VectorXd y(Nn);

            for (label cI = 0; cI < Nn; cI++)
            {
                if (cI < curStencil.size())
                {
                    y[cI] = D[curStencil[cI]][cmptI];
                }
                else
                {
                    const label i = cI - curStencil.size();
                    const label globalFaceID = stencilsBoundaryFaces[cellI][i];

                    forAll(boundaryMesh, patchI)
                    {
                        if
                        (
                            includePatchInStencils_[patchI]
                         && boundaryMesh[patchI].type()
                         != emptyPolyPatch::typeName
                         && !boundaryMesh[patchI].coupled()
                        )
                        {
                            const label start = boundaryMesh[patchI].start();
                            const label nFaces = boundaryMesh[patchI].nFaces();

                            if
                            (
                                globalFaceID >= start
                             && globalFaceID < start + nFaces
                            )
                            {
                                const label k = globalFaceID - start;
                                y[cI] = D.boundaryField()[patchI][k][cmptI];
                            }
                        }
                    }
                }
            }

            // Compute y_hat = W^{1/2}*y
            const Eigen::VectorXd yhat =
                sqrtW[cellI].diagonal().array()*y.array();

            // Compute Q^T W y = Q_hat^T * y_hat
            const Eigen::VectorXd QTWy = Qhat[cellI].transpose()*yhat;

            // Solve for A
            const Eigen::VectorXd z = cholesky[cellI].matrixL().solve(QTWy);
            const Eigen::VectorXd A = cholesky[cellI].matrixU().solve(z);

            // Extract gradient components from A and assign the gradient field
            // Careful: they are column-wise
            gradD[cellI][3*0 + cmptI] = A[1];
            gradD[cellI][3*1 + cmptI] = A[2];
            gradD[cellI][3*2 + cmptI] = A[3];
        }
    }

    gradD.correctBoundaryConditions();

    if (debug)
    {
        Info<< "higherOrderGrad::gradCholesky(...): end" << endl;
    }

    return tgradD;
}


tmp<volTensorField> higherOrderGrad::gradGlobalCholesky
(
    const volVectorField& D
) const
{
    if (debug)
    {
        InfoInFunction
            << "start" << endl;
    }

    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<volTensorField> tgradD
    (
        new volTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero),
           "zeroGradient"
        )
    );
    volTensorField& gradD = tgradD.ref();

    const List<labelList>& stencils = globalCellStencils();
    const List<Eigen::LLT<Eigen::MatrixXd>>& cholesky = this->cholesky();
    const List<Eigen::MatrixXd>& Qhat = this->Qhat();
    const List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>& sqrtW =
        this->sqrtW();
    const vectorField& DI = D;

    // Collect DI for off-processor cells in the stencils
    Map<vector> globalDI;
    requestGlobalStencilData(DI, globalDI);

    forAll(stencils, localCellI)
    {
        const labelList& curStencil = stencils[localCellI];
        const label Nn = curStencil.size();

        for (label cmptI = 0; cmptI < vector::nComponents; ++cmptI)
        {
            // Prepare right-hand side vector (y) for Cholesky decomposition
            Eigen::VectorXd y(Nn);

            for (label cI = 0; cI < Nn; cI++)
            {
                const label globalCellID = curStencil[cI];

                if (globalCells_.isLocal(globalCellID))
                {
                    const label localCellID =
                        globalCells_.toLocal(globalCellID);

                    y[cI] = DI[localCellID][cmptI];
                }
                else
                {
                    y[cI] = globalDI[globalCellID][cmptI];
                }
            }

            // Compute y_hat = W^{1/2}*y
            const Eigen::VectorXd yhat =
                sqrtW[localCellI].diagonal().array()*y.array();

            // Compute Q^T W y = Q_hat^T * y_hat
            const Eigen::VectorXd QTWy = Qhat[localCellI].transpose()*yhat;

            // Solve for A
            const Eigen::VectorXd z = cholesky[localCellI].matrixL().solve(QTWy);
            const Eigen::VectorXd A = cholesky[localCellI].matrixU().solve(z);

            // Extract gradient components from A and assign the gradient field
            // Careful: they are column-wise
            gradD[localCellI][3*0 + cmptI] = A[1];
            gradD[localCellI][3*1 + cmptI] = A[2];
            gradD[localCellI][3*2 + cmptI] = A[3];
        }
    }

    gradD.correctBoundaryConditions();

    if (debug)
    {
        InfoInFunction
            << "end" << endl;
    }

    return tgradD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
