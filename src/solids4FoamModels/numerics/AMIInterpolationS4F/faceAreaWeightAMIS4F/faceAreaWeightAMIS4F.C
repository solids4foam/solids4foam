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

#ifdef OPENFOAM_COM

#include "faceAreaWeightAMIS4F.H"
#include "profiling.H"
//#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceAreaWeightAMIS4F, 0);
    addToRunTimeSelectionTable(AMIInterpolationS4F, faceAreaWeightAMIS4F, dict);
    addToRunTimeSelectionTable(AMIInterpolationS4F, faceAreaWeightAMIS4F, component);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

/*
    if (debug)
    {
        static label nAMI = 0;

        // Write out triangulated surfaces as OBJ files
        OBJstream srcTriObj("srcTris_" + Foam::name(nAMI) + ".obj");
        const pointField& srcPts = src.points();
        forAll(srcTris_, facei)
        {
            const DynamicList<face>& faces = srcTris_[facei];
            for (const face& f : faces)
            {
                srcTriObj.write
                (
                    triPointRef(srcPts[f[0]], srcPts[f[1]], srcPts[f[2]])
                );
            }
        }

        OBJstream tgtTriObj("tgtTris_" + Foam::name(nAMI) + ".obj");
        const pointField& tgtPts = tgt.points();
        forAll(tgtTris_, facei)
        {
            const DynamicList<face>& faces = tgtTris_[facei];
            for (const face& f : faces)
            {
                tgtTriObj.write
                (
                    triPointRef(tgtPts[f[0]], tgtPts[f[1]], tgtPts[f[2]])
                );
            }
        }

        ++nAMI;
    }
*/


void Foam::faceAreaWeightAMIS4F::calcAddressing
(
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<point>>& srcCtr,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght,
    label srcFacei,
    label tgtFacei
)
{
    addProfiling(ami, "faceAreaWeightAMIS4F::calcAddressing");

    // Construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nFacesRemaining = srcAddr.size();

    // List of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // List of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    // List to keep track of tgt faces used to seed src faces
    labelList seedFaces(nFacesRemaining, -1);
    seedFaces[srcFacei] = tgtFacei;

    // List to keep track of whether src face can be mapped
    bitSet mapFlag(nFacesRemaining, true);

    // Reset starting seed
    label startSeedi = 0;

    bool continueWalk = true;
    DynamicList<label> nonOverlapFaces;
    do
    {
        nbrFaces.clear();
        visitedFaces.clear();

        // Do advancing front starting from srcFacei,tgtFacei
        bool faceProcessed = processSourceFace
        (
            srcFacei,
            tgtFacei,

            nbrFaces,
            visitedFaces,

            srcAddr,
            srcWght,
            srcCtr,
            tgtAddr,
            tgtWght
        );

        mapFlag.unset(srcFacei);

        if (!faceProcessed)
        {
            nonOverlapFaces.append(srcFacei);
        }

        // Choose new src face from current src face neighbour
        continueWalk = setNextFaces
        (
            startSeedi,
            srcFacei,
            tgtFacei,
            mapFlag,
            seedFaces,
            visitedFaces,
            requireMatch_ && (lowWeightCorrection_ < 0)
            // pass in nonOverlapFaces for failed tree search?
        );
    } while (continueWalk);

    srcNonOverlap_.transfer(nonOverlapFaces);
}


bool Foam::faceAreaWeightAMIS4F::processSourceFace
(
    const label srcFacei,
    const label tgtStartFacei,

    // list of tgt face neighbour faces
    DynamicList<label>& nbrFaces,
    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label>& visitedFaces,

    // temporary storage for addressing, weights and centroid
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<point>>& srcCtr,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght
)
{
    addProfiling(ami, "faceAreaWeightAMIS4F::processSourceFace");

    if (tgtStartFacei == -1)
    {
        return false;
    }

    const auto& tgtPatch = this->tgtPatch();

    // append initial target face and neighbours
    nbrFaces.append(tgtStartFacei);
    appendNbrFaces(tgtStartFacei, tgtPatch, visitedFaces, nbrFaces);

    bool faceProcessed = false;

    label maxNeighbourFaces = nbrFaces.size();

    do
    {
        // process new target face
        label tgtFacei = nbrFaces.remove();
        visitedFaces.append(tgtFacei);

        scalar interArea = 0;
        vector interCentroid(Zero);
        calcInterArea(srcFacei, tgtFacei, interArea, interCentroid);

        // store when intersection fractional area > tolerance
        if (interArea/srcMagSf_[srcFacei] > faceAreaIntersect::tolerance())
        {
            srcAddr[srcFacei].append(tgtFacei);
            srcWght[srcFacei].append(interArea);
            srcCtr[srcFacei].append(interCentroid);

            tgtAddr[tgtFacei].append(srcFacei);
            tgtWght[tgtFacei].append(interArea);

            appendNbrFaces(tgtFacei, tgtPatch, visitedFaces, nbrFaces);

            faceProcessed = true;

            maxNeighbourFaces = max(maxNeighbourFaces, nbrFaces.size());
        }

    } while (nbrFaces.size() > 0);

    if (debug > 1)
    {
        DebugVar(maxNeighbourFaces);
    }

    return faceProcessed;
}


bool Foam::faceAreaWeightAMIS4F::setNextFaces
(
    label& startSeedi,
    label& srcFacei,
    label& tgtFacei,
    const bitSet& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    const bool errorOnNotFound
) const
{
    addProfiling(ami, "faceAreaWeightAMIS4F::setNextFaces");

    if (mapFlag.count() == 0)
    {
        // No more faces to map
        return false;
    }

    const labelList& srcNbrFaces = this->srcPatch().faceFaces()[srcFacei];

    // Initialise tgtFacei
    tgtFacei = -1;

    // Set possible seeds for later use
    bool valuesSet = false;
    for (label faceS: srcNbrFaces)
    {
        if (mapFlag.test(faceS) && seedFaces[faceS] == -1)
        {
            for (label faceT : visitedFaces)
            {
                const scalar threshold =
                    srcMagSf_[faceS]*faceAreaIntersect::tolerance();

                // Store when intersection fractional area > threshold
                if (overlaps(faceS, faceT, threshold))
                {
                    seedFaces[faceS] = faceT;

                    if (!valuesSet)
                    {
                        srcFacei = faceS;
                        tgtFacei = faceT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    if (valuesSet)
    {
        return true;
    }

    // Set next src and tgt faces if not set above
    // - try to use existing seed
    label facei = startSeedi;
    if (!mapFlag.test(startSeedi))
    {
        facei = mapFlag.find_next(facei);
    }
    const label startSeedi0 = facei;

    bool foundNextSeed = false;
    while (facei != -1)
    {
        if (!foundNextSeed)
        {
            startSeedi = facei;
            foundNextSeed = true;
        }

        if (seedFaces[facei] != -1)
        {
            srcFacei = facei;
            tgtFacei = seedFaces[facei];

            return true;
        }

        facei = mapFlag.find_next(facei);
    }

    // Perform new search to find match
    if (debug)
    {
        Pout<< "Advancing front stalled: searching for new "
            << "target face" << endl;
    }

    facei = startSeedi0;
    while (facei != -1)
    {
        srcFacei = facei;
        tgtFacei = findTargetFace(srcFacei, visitedFaces);

        if (tgtFacei >= 0)
        {
            return true;
        }

        facei = mapFlag.find_next(facei);
    }

    if (errorOnNotFound)
    {
        FatalErrorInFunction
            << "Unable to set target face for source face " << srcFacei
            << abort(FatalError);
    }

    return false;
}


void Foam::faceAreaWeightAMIS4F::calcInterArea
(
    const label srcFacei,
    const label tgtFacei,
    scalar& area,
    vector& centroid
) const
{
    addProfiling(ami, "faceAreaWeightAMIS4F::interArea");

    // Quick reject if either face has zero area
    if
    (
        (srcMagSf_[srcFacei] < ROOTVSMALL)
     || (tgtMagSf_[tgtFacei] < ROOTVSMALL)
    )
    {
        return;
    }

    const auto& srcPatch = this->srcPatch();
    const auto& tgtPatch = this->tgtPatch();

    const pointField& srcPoints = srcPatch.points();
    const pointField& tgtPoints = tgtPatch.points();

    // Create intersection object
    faceAreaIntersect inter
    (
        srcPoints,
        tgtPoints,
        srcTris_[srcFacei],
        tgtTris_[tgtFacei],
        reverseTarget_,
        AMIInterpolationS4F::cacheIntersections_
    );

    // Crude resultant norm
    vector n(-srcPatch.faceNormals()[srcFacei]);
    if (reverseTarget_)
    {
        n -= tgtPatch.faceNormals()[tgtFacei];
    }
    else
    {
        n += tgtPatch.faceNormals()[tgtFacei];
    }
    scalar magN = mag(n);

    const face& src = srcPatch[srcFacei];
    const face& tgt = tgtPatch[tgtFacei];

    if (magN > ROOTVSMALL)
    {
        inter.calc(src, tgt, n/magN, area, centroid);
    }
    else
    {
        WarningInFunction
            << "Invalid normal for source face " << srcFacei
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFacei
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << endl;
    }

    // if (AMIInterpolationS4F::cacheIntersections_ && debug)
    // {
    //     static OBJstream tris("intersectionTris.obj");
    //     const auto& triPts = inter.triangles();
    //     for (const auto& tp : triPts)
    //     {
    //         tris.write(triPointRef(tp[0], tp[1], tp[2]), false);
    //     }
    // }

    if ((debug > 1) && (area > 0))
    {
        writeIntersectionOBJ(area, src, tgt, srcPoints, tgtPoints);
    }
}


bool Foam::faceAreaWeightAMIS4F::overlaps
(
    const label srcFacei,
    const label tgtFacei,
    const scalar threshold
) const
{
    // Quick reject if either face has zero area
    if
    (
        (srcMagSf_[srcFacei] < ROOTVSMALL)
     || (tgtMagSf_[tgtFacei] < ROOTVSMALL)
    )
    {
        return false;
    }

    const auto& srcPatch = this->srcPatch();
    const auto& tgtPatch = this->tgtPatch();

    const pointField& srcPoints = srcPatch.points();
    const pointField& tgtPoints = tgtPatch.points();

    faceAreaIntersect inter
    (
        srcPoints,
        tgtPoints,
        srcTris_[srcFacei],
        tgtTris_[tgtFacei],
        reverseTarget_,
        AMIInterpolationS4F::cacheIntersections_
    );

    // Crude resultant norm
    vector n(-srcPatch.faceNormals()[srcFacei]);
    if (reverseTarget_)
    {
        n -= tgtPatch.faceNormals()[tgtFacei];
    }
    else
    {
        n += tgtPatch.faceNormals()[tgtFacei];
    }
    scalar magN = mag(n);

    const face& src = srcPatch[srcFacei];
    const face& tgt = tgtPatch[tgtFacei];

    if (magN > ROOTVSMALL)
    {
        return inter.overlaps(src, tgt, n/magN, threshold);
    }
    else
    {
        WarningInFunction
            << "Invalid normal for source face " << srcFacei
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFacei
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << ", srcN " << srcPatch.faceNormals()[srcFacei]
            << ", tgtN " << tgtPatch.faceNormals()[srcFacei]
            << ", n " << n
            << ", magN " << magN
            << endl;
    }

    return false;
}


void Foam::faceAreaWeightAMIS4F::restartUncoveredSourceFace
(
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<point>>& srcCtr,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght
)
{
    addProfiling(ami, "faceAreaWeightAMIS4F::restartUncoveredSourceFace");

    // Note: exclude faces in srcNonOverlap_ for ACMI?

    label nBelowMinWeight = 0;
    const scalar minWeight = 0.95;

    // List of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // List of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    const auto& srcPatch = this->srcPatch();

    forAll(srcWght, srcFacei)
    {
        const scalar s = sum(srcWght[srcFacei]);
        const scalar t = s/srcMagSf_[srcFacei];

        if (t < minWeight)
        {
            ++nBelowMinWeight;

            const face& f = srcPatch[srcFacei];

            forAll(f, fpi)
            {
                const label tgtFacei =
                    findTargetFace(srcFacei, srcAddr[srcFacei], fpi);

                if (tgtFacei != -1)
                {
                    nbrFaces.clear();
                    visitedFaces = srcAddr[srcFacei];

                    (void)processSourceFace
                    (
                        srcFacei,
                        tgtFacei,

                        nbrFaces,
                        visitedFaces,

                        srcAddr,
                        srcWght,
                        srcCtr,
                        tgtAddr,
                        tgtWght
                    );
                }
            }
        }
    }

    if (debug && nBelowMinWeight)
    {
        WarningInFunction
            << "Restarted search on " << nBelowMinWeight
            << " faces since sum of weights < " << minWeight
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceAreaWeightAMIS4F::faceAreaWeightAMIS4F
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    advancingFrontAMIS4F(dict, reverseTarget),
    restartUncoveredSourceFace_
    (
        dict.getOrDefault("restartUncoveredSourceFace", true)
    )
{}


Foam::faceAreaWeightAMIS4F::faceAreaWeightAMIS4F
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection,
    const faceAreaIntersect::triangulationMode triMode,
    const bool restartUncoveredSourceFace
)
:
    advancingFrontAMIS4F
    (
        requireMatch,
        reverseTarget,
        lowWeightCorrection,
        triMode
    ),
    restartUncoveredSourceFace_(restartUncoveredSourceFace)
{}


Foam::faceAreaWeightAMIS4F::faceAreaWeightAMIS4F(const faceAreaWeightAMIS4F& ami)
:
    advancingFrontAMIS4F(ami),
    restartUncoveredSourceFace_(ami.restartUncoveredSourceFace_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceAreaWeightAMIS4F::calculate
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (upToDate_)
    {
        return false;
    }

    addProfiling(ami, "faceAreaWeightAMIS4F::calculate");

    advancingFrontAMIS4F::calculate(srcPatch, tgtPatch, surfPtr);

    label srcFacei = 0;
    label tgtFacei = 0;

    bool ok = initialiseWalk(srcFacei, tgtFacei);

    srcCentroids_.setSize(srcAddress_.size());

    const auto& src = this->srcPatch();
    const auto& tgt = this->tgtPatch(); // might be the extended patch!

    // Temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(src.size());
    List<DynamicList<scalar>> srcWght(srcAddr.size());
    List<DynamicList<point>> srcCtr(srcAddr.size());
    List<DynamicList<label>> tgtAddr(tgt.size());
    List<DynamicList<scalar>> tgtWght(tgtAddr.size());

    if (ok)
    {
        calcAddressing
        (
            srcAddr,
            srcWght,
            srcCtr,
            tgtAddr,
            tgtWght,
            srcFacei,
            tgtFacei
        );

        if (debug && !srcNonOverlap_.empty())
        {
            Pout<< "    AMI: " << srcNonOverlap_.size()
                << " non-overlap faces identified"
                << endl;
        }

        // Check for badly covered faces
        if (restartUncoveredSourceFace_) // && requireMatch_???
        {
            restartUncoveredSourceFace
            (
                srcAddr,
                srcWght,
                srcCtr,
                tgtAddr,
                tgtWght
            );
        }
    }

    // Transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress_[i].transfer(srcAddr[i]);
        srcWeights_[i].transfer(srcWght[i]);
        srcCentroids_[i].transfer(srcCtr[i]);
    }

    forAll(tgtAddr, i)
    {
        tgtAddress_[i].transfer(tgtAddr[i]);
        tgtWeights_[i].transfer(tgtWght[i]);
    }

    if (distributed())
    {
        const standAlonePatch& srcPatch0 = this->srcPatch0();
        const standAlonePatch& tgtPatch0 = this->tgtPatch0();

        // Create global indexing for each original patch
        globalIndex globalSrcFaces(srcPatch0.size());
        globalIndex globalTgtFaces(tgtPatch0.size());

        for (labelList& addressing : srcAddress_)
        {
            for (label& addr : addressing)
            {
                addr = extendedTgtFaceIDs_[addr];
            }
        }

        for (labelList& addressing : tgtAddress_)
        {
            globalSrcFaces.inplaceToGlobal(addressing);
        }

        // Send data back to originating procs. Note that contributions
        // from different processors get added (ListOps::appendEqOp)

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtPatch0.size(),
            extendedTgtMapPtr_->constructMap(),
            false,                      // has flip
            extendedTgtMapPtr_->subMap(),
            false,                      // has flip
            tgtAddress_,
            labelList(),
            ListOps::appendEqOp<label>(),
            flipOp()                    // flip operation
        );

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtPatch0.size(),
            extendedTgtMapPtr_->constructMap(),
            false,
            extendedTgtMapPtr_->subMap(),
            false,
            tgtWeights_,
            scalarList(),
            ListOps::appendEqOp<scalar>(),
            flipOp()
        );

        // Note: using patch face areas calculated by the AMI method
        extendedTgtMapPtr_->reverseDistribute(tgtPatch0.size(), tgtMagSf_);

        // Cache maps and reset addresses
        List<Map<label>> cMapSrc;
        srcMapPtr_.reset
        (
            new mapDistribute(globalSrcFaces, tgtAddress_, cMapSrc)
        );

        List<Map<label>> cMapTgt;
        tgtMapPtr_.reset
        (
            new mapDistribute(globalTgtFaces, srcAddress_, cMapTgt)
        );
    }

    // Convert the weights from areas to normalised values
    normaliseWeights(conformal(), true);

    upToDate_ = true;

    return upToDate_;
}


void Foam::faceAreaWeightAMIS4F::write(Ostream& os) const
{
    advancingFrontAMIS4F::write(os);

    if (restartUncoveredSourceFace_)
    {
        os.writeEntry
        (
            "restartUncoveredSourceFace",
            restartUncoveredSourceFace_
        );
    }
}

#endif // OPENFOAM_COM

// ************************************************************************* //
