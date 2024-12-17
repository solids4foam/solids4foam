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

#include "AMIInterpolationS4F.H"
#include "meshTools.H"
#include "mapDistribute.H"
#include "flipOp.H"
#include "profiling.H"
#include "triPointRef.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AMIInterpolationS4F, 0);
    defineRunTimeSelectionTable(AMIInterpolationS4F, dict);
    defineRunTimeSelectionTable(AMIInterpolationS4F, component);
}

bool Foam::AMIInterpolationS4F::cacheIntersections_ = false;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indexedOctree<Foam::AMIInterpolationS4F::treeType>>
Foam::AMIInterpolationS4F::createTree
(
    const standAlonePatch& patch
) const
{
    treeBoundBox bb(patch.points(), patch.meshPoints());
    bb.inflate(0.01);

    return autoPtr<indexedOctree<treeType>>::New
    (
        treeType
        (
            false,
            patch,
            indexedOctree<treeType>::perturbTol()
        ),
        bb,                         // overall search domain
        8,                          // maxLevel
        10,                         // leaf size
        3.0                         // duplicity
    );
}


Foam::label Foam::AMIInterpolationS4F::calcDistribution
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch
) const
{
    label proci = 0;

    if (Pstream::parRun())
    {
        labelList facesPresentOnProc(Pstream::nProcs(), Zero);
        if ((srcPatch.size() > 0) || (tgtPatch.size() > 0))
        {
            facesPresentOnProc[Pstream::myProcNo()] = 1;
        }
        else
        {
            facesPresentOnProc[Pstream::myProcNo()] = 0;
        }

        Pstream::gatherList(facesPresentOnProc);
        Pstream::scatterList(facesPresentOnProc);

        label nHaveFaces = sum(facesPresentOnProc);

        if (nHaveFaces > 1)
        {
            proci = -1;
            if (debug)
            {
                InfoInFunction
                    << "AMI split across multiple processors" << endl;
            }
        }
        else if (nHaveFaces == 1)
        {
            proci = facesPresentOnProc.find(1);
            if (debug)
            {
                InfoInFunction
                    << "AMI local to processor" << proci << endl;
            }
        }
    }


    // Either not parallel or no faces on any processor
    return proci;
}


void Foam::AMIInterpolationS4F::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    addProfiling(ami, "AMIInterpolationS4F::projectPointsToSurface");

    DebugInfo<< "AMI: projecting points to surface" << endl;

    List<pointIndexHit> nearInfo;

    surf.findNearest(pts, scalarField(pts.size(), GREAT), nearInfo);

    label nMiss = 0;
    forAll(nearInfo, i)
    {
        const pointIndexHit& pi = nearInfo[i];

        if (pi.hit())
        {
            pts[i] = pi.hitPoint();
        }
        else
        {
            // Point remains unchanged
            ++nMiss;
        }
    }

    if (nMiss > 0)
    {
        FatalErrorInFunction
            << "Error projecting points to surface: "
            << nMiss << " faces could not be determined"
            << abort(FatalError);
    }
}


void Foam::AMIInterpolationS4F::normaliseWeights
(
    const scalarList& patchAreas,
    const word& patchName,
    const labelListList& addr,
    scalarListList& wght,
    scalarField& wghtSum,
    const bool conformal,
    const bool output,
    const scalar lowWeightTol
)
{
    addProfiling(ami, "AMIInterpolationS4F::normaliseWeights");

    // Normalise the weights
    wghtSum.setSize(wght.size(), 0.0);
    //label nLowWeight = 0;

    forAll(wght, facei)
    {
        scalarList& w = wght[facei];

        if (w.size())
        {
            scalar denom = patchAreas[facei];

            scalar s = sum(w);
            scalar t = s/denom;
            if (conformal)
            {
                denom = s;
            }

            forAll(w, i)
            {
                w[i] /= denom;
            }

            wghtSum[facei] = t;
            // if (t < lowWeightTol)
            // {
            //     ++nLowWeight;
            // }
        }
        else
        {
            wghtSum[facei] = 0;
        }
    }

    // if (output)
    // {
    //     const label nFace = returnReduce(wght.size(), sumOp<label>());

    //     if (nFace)
    //     {
    //         Info<< indent
    //             << "AMI: Patch " << patchName
    //             << " sum(weights)"
    //             << " min:" << gMin(wghtSum)
    //             << " max:" << gMax(wghtSum)
    //             << " average:" << gAverage(wghtSum) << nl;

    //         const label nLow = returnReduce(nLowWeight, sumOp<label>());

    //         if (nLow)
    //         {
    //             Info<< indent
    //                 << "AMI: Patch " << patchName
    //                 << " identified " << nLow
    //                 << " faces with weights less than " << lowWeightTol
    //                 << endl;
    //         }
    //     }
    // }
}


void Foam::AMIInterpolationS4F::agglomerate
(
    const autoPtr<mapDistribute>& targetMapPtr,
    const scalarList& fineSrcMagSf,
    const labelListList& fineSrcAddress,
    const scalarListList& fineSrcWeights,

    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing,

    scalarList& srcMagSf,
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarField& srcWeightsSum,
    autoPtr<mapDistribute>& tgtMap
)
{
    addProfiling(ami, "AMIInterpolationS4F::agglomerate");

    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label targetCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    // Agglomerate face areas
    {
        srcMagSf.setSize(sourceRestrictAddressing.size(), 0.0);
        forAll(sourceRestrictAddressing, facei)
        {
            label coarseFacei = sourceRestrictAddressing[facei];
            srcMagSf[coarseFacei] += fineSrcMagSf[facei];
        }
    }

    // Agglomerate weights and indices
    if (targetMapPtr)
    {
        const mapDistribute& map = *targetMapPtr;

        // Get all restriction addressing.
        labelList allRestrict(targetRestrictAddressing);
        map.distribute(allRestrict);

        // So now we have agglomeration of the target side in
        // allRestrict:
        //  0..size-1 : local agglomeration (= targetRestrictAddressing
        //              (but potentially permutated))
        //  size..    : agglomeration data from other processors


        // The trickiness in this algorithm is finding out the compaction
        // of the remote data (i.e. allocation of the coarse 'slots'). We could
        // either send across the slot compaction maps or just make sure
        // that we allocate the slots in exactly the same order on both sending
        // and receiving side (e.g. if the submap is set up to send 4 items,
        // the constructMap is also set up to receive 4 items.


        // Short note about the various types of indices:
        // - face indices : indices into the geometry.
        // - coarse face indices : how the faces get agglomerated
        // - transferred data : how mapDistribute sends/receives data,
        // - slots : indices into data after distribution (e.g. stencil,
        //           srcAddress/tgtAddress). Note: for fully local addressing
        //           the slots are equal to face indices.
        // A mapDistribute has:
        // - a subMap : these are face indices
        // - a constructMap : these are from 'transferred-data' to slots

        labelListList tgtSubMap(Pstream::nProcs());

        // Local subMap is just identity
        {
            tgtSubMap[Pstream::myProcNo()] = identity(targetCoarseSize);
        }

        forAll(map.subMap(), proci)
        {
            if (proci != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element.
                // The important bit is to loop over the data (and hand out
                // compact indices ) in 'transferred data' order. This
                // guarantees that we're doing exactly the
                // same on sending and receiving side - e.g. the fourth element
                // in the subMap is the fourth element received in the
                // constructMap

                const labelList& elems = map.subMap()[proci];
                const labelList& elemsMap =
                    map.constructMap()[Pstream::myProcNo()];
                labelList& newSubMap = tgtSubMap[proci];
                newSubMap.setSize(elems.size());

                labelList oldToNew(targetCoarseSize, -1);
                label newi = 0;

                for (const label elemi : elems)
                {
                    label fineElem = elemsMap[elemi];
                    label coarseElem = allRestrict[fineElem];
                    if (oldToNew[coarseElem] == -1)
                    {
                        oldToNew[coarseElem] = newi;
                        newSubMap[newi] = coarseElem;
                        ++newi;
                    }
                }
                newSubMap.setSize(newi);
            }
        }

        // Reconstruct constructMap by combining entries. Note that order
        // of handing out indices should be the same as loop above to compact
        // the sending map

        labelListList tgtConstructMap(Pstream::nProcs());

        // Local constructMap is just identity
        {
            tgtConstructMap[Pstream::myProcNo()] = identity(targetCoarseSize);
        }

        labelList tgtCompactMap(map.constructSize());

        {
            // Note that in special cases (e.g. 'appending' two AMIs) the
            // local size after distributing can be longer than the number
            // of faces. I.e. it duplicates elements.
            // Since we don't know this size instead we loop over all
            // reachable elements (using the local constructMap)

            const labelList& elemsMap = map.constructMap()[Pstream::myProcNo()];
            for (const label fineElem : elemsMap)
            {
                label coarseElem = allRestrict[fineElem];
                tgtCompactMap[fineElem] = coarseElem;
            }
        }

        label compacti = targetCoarseSize;

        // Compact data from other processors
        forAll(map.constructMap(), proci)
        {
            if (proci != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element. All
                // elements now are remote data so we cannot use any local
                // data here - use allRestrict instead.
                const labelList& elems = map.constructMap()[proci];

                labelList& newConstructMap = tgtConstructMap[proci];
                newConstructMap.setSize(elems.size());

                if (elems.size())
                {
                    // Get the maximum target coarse size for this set of
                    // received data.
                    label remoteTargetCoarseSize = labelMin;
                    for (const label elemi : elems)
                    {
                        remoteTargetCoarseSize = max
                        (
                            remoteTargetCoarseSize,
                            allRestrict[elemi]
                        );
                    }
                    remoteTargetCoarseSize += 1;

                    // Combine locally data coming from proci
                    labelList oldToNew(remoteTargetCoarseSize, -1);
                    label newi = 0;

                    for (const label fineElem : elems)
                    {
                        // fineElem now points to section from proci
                        label coarseElem = allRestrict[fineElem];
                        if (oldToNew[coarseElem] == -1)
                        {
                            oldToNew[coarseElem] = newi;
                            tgtCompactMap[fineElem] = compacti;
                            newConstructMap[newi] = compacti++;
                            ++newi;
                        }
                        else
                        {
                            // Get compact index
                            label compacti = oldToNew[coarseElem];
                            tgtCompactMap[fineElem] = newConstructMap[compacti];
                        }
                    }
                    newConstructMap.setSize(newi);
                }
            }
        }

        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, facei)
        {
            // All the elements contributing to facei. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[facei];
            const scalarList& weights = fineSrcWeights[facei];
            const scalar fineArea = fineSrcMagSf[facei];

            label coarseFacei = sourceRestrictAddressing[facei];

            labelList& newElems = srcAddress[coarseFacei];
            scalarList& newWeights = srcWeights[coarseFacei];

            forAll(elems, i)
            {
                label elemi = elems[i];
                label coarseElemi = tgtCompactMap[elemi];

                label index = newElems.find(coarseElemi);
                if (index == -1)
                {
                    newElems.append(coarseElemi);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }

        tgtMap.reset
        (
            new mapDistribute
            (
                compacti,
                std::move(tgtSubMap),
                std::move(tgtConstructMap)
            )
        );
    }
    else
    {
        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);

        forAll(fineSrcAddress, facei)
        {
            // All the elements contributing to facei. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[facei];
            const scalarList& weights = fineSrcWeights[facei];
            const scalar fineArea = fineSrcMagSf[facei];

            label coarseFacei = sourceRestrictAddressing[facei];

            labelList& newElems = srcAddress[coarseFacei];
            scalarList& newWeights = srcWeights[coarseFacei];

            forAll(elems, i)
            {
                const label elemi = elems[i];
                const label coarseElemi = targetRestrictAddressing[elemi];

                const label index = newElems.find(coarseElemi);
                if (index == -1)
                {
                    newElems.append(coarseElemi);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
        }
    }

    // Weights normalisation
    normaliseWeights
    (
        srcMagSf,
        "source",
        srcAddress,
        srcWeights,
        srcWeightsSum,
        true,
        false,
        -1
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AMIInterpolationS4F::AMIInterpolationS4F
(
    const dictionary& dict,
    const bool reverseTarget,
     const bool useGlobalPolyPatch
)
:
    requireMatch_(dict.getOrDefault("requireMatch", true)),
    reverseTarget_(dict.getOrDefault("reverseTarget", reverseTarget)),
    lowWeightCorrection_(dict.getOrDefault<scalar>("lowWeightCorrection", -1)),
    singlePatchProc_(-999),
//	slavePointAddressingPtr_(NULL),//JN: from extendedggi
//	slavePointDistancePtr_(NULL),//JN: from extendedggi
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{}

Foam::AMIInterpolationS4F::AMIInterpolationS4F
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch,
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection,
     const bool useGlobalPolyPatch
)
:
    requireMatch_(requireMatch),
    reverseTarget_(reverseTarget),
    lowWeightCorrection_(lowWeightCorrection),
    singlePatchProc_(-999),
//	slavePointAddressingPtr_(NULL),//JN: from extendedggi
//	slavePointDistancePtr_(NULL),//JN: from extendedggi
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcPatchPts_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtPatchPts_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{
	calculate(srcPatch, tgtPatch);
}

Foam::AMIInterpolationS4F::AMIInterpolationS4F
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection,
     const bool useGlobalPolyPatch
)
:
    requireMatch_(requireMatch),
    reverseTarget_(reverseTarget),
    lowWeightCorrection_(lowWeightCorrection),
    singlePatchProc_(-999),
//	slavePointAddressingPtr_(NULL),//JN: from extendedggi
//	slavePointDistancePtr_(NULL),//JN: from extendedggi
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcPatchPts_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtPatchPts_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{}


Foam::AMIInterpolationS4F::AMIInterpolationS4F
(
    const AMIInterpolationS4F& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing,
     const bool useGlobalPolyPatch
)
:
    requireMatch_(fineAMI.requireMatch_),
    reverseTarget_(fineAMI.reverseTarget_),
    lowWeightCorrection_(-1.0),
    singlePatchProc_(fineAMI.singlePatchProc_),
//	slavePointAddressingPtr_(NULL),//JN: from extendedggi
//	slavePointDistancePtr_(NULL),//JN: from extendedggi
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    srcPatchPts_(),
    srcMapPtr_(nullptr),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    tgtPatchPts_(),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label neighbourCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    if (debug & 2)
    {
        Pout<< "AMI: Creating addressing and weights as agglomeration of AMI :"
            << " source:" << fineAMI.srcAddress().size()
            << " target:" << fineAMI.tgtAddress().size()
            << " coarse source size:" << sourceCoarseSize
            << " neighbour source size:" << neighbourCoarseSize
            << endl;
    }

    if
    (
        fineAMI.srcAddress().size() != sourceRestrictAddressing.size()
     || fineAMI.tgtAddress().size() != targetRestrictAddressing.size()
    )
    {
        FatalErrorInFunction
            << "Size mismatch." << nl
            << "Source patch size:" << fineAMI.srcAddress().size() << nl
            << "Source agglomeration size:"
            << sourceRestrictAddressing.size() << nl
            << "Target patch size:" << fineAMI.tgtAddress().size() << nl
            << "Target agglomeration size:"
            << targetRestrictAddressing.size()
            << exit(FatalError);
    }


    // Agglomerate addresses and weights

    agglomerate
    (
        fineAMI.tgtMapPtr_,
        fineAMI.srcMagSf(),
        fineAMI.srcAddress(),
        fineAMI.srcWeights(),

        sourceRestrictAddressing,
        targetRestrictAddressing,

        srcMagSf_,
        srcAddress_,
        srcWeights_,
        srcWeightsSum_,
        tgtMapPtr_
    );

    agglomerate
    (
        fineAMI.srcMapPtr_,
        fineAMI.tgtMagSf(),
        fineAMI.tgtAddress(),
        fineAMI.tgtWeights(),

        targetRestrictAddressing,
        sourceRestrictAddressing,

        tgtMagSf_,
        tgtAddress_,
        tgtWeights_,
        tgtWeightsSum_,
        srcMapPtr_
    );
}


Foam::AMIInterpolationS4F::AMIInterpolationS4F(
	const AMIInterpolationS4F& ami,
     const bool useGlobalPolyPatch
)
:
    requireMatch_(ami.requireMatch_),
    reverseTarget_(ami.reverseTarget_),
    lowWeightCorrection_(ami.lowWeightCorrection_),
    singlePatchProc_(ami.singlePatchProc_),
//	slavePointAddressingPtr_(NULL),//JN: from extendedggi
//	slavePointDistancePtr_(NULL),//JN: from extendedggi
    srcMagSf_(ami.srcMagSf_),
    srcAddress_(ami.srcAddress_),
    srcWeights_(ami.srcWeights_),
    srcWeightsSum_(ami.srcWeightsSum_),
    srcCentroids_(ami.srcCentroids_),
    srcMapPtr_(nullptr),
    tgtMagSf_(ami.tgtMagSf_),
    tgtAddress_(ami.tgtAddress_),
    tgtWeights_(ami.tgtWeights_),
    tgtWeightsSum_(ami.tgtWeightsSum_),
    tgtCentroids_(ami.tgtCentroids_),
    tgtMapPtr_(nullptr),
    upToDate_(false),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::AMIInterpolationS4F::calculate
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

    addProfiling(ami, "AMIInterpolationS4F::calculate");

    if (surfPtr)
    {
        srcPatchPts_ = srcPatch.points();
        projectPointsToSurface(surfPtr(), srcPatchPts_);
        tsrcPatch0_ = refPtr<standAlonePatch>::New
        (
            SubList<face>(srcPatch),
            srcPatchPts_
        );

        tgtPatchPts_ = tgtPatch.points();
        projectPointsToSurface(surfPtr(), tgtPatchPts_);
        ttgtPatch0_ = refPtr<standAlonePatch>::New
        (
            SubList<face>(tgtPatch),
            tgtPatchPts_
        );
    }
    else
    {
        tsrcPatch0_.cref(srcPatch);
        ttgtPatch0_.cref(tgtPatch);
    }

    label srcTotalSize = returnReduce(srcPatch.size(), sumOp<label>());
    //label tgtTotalSize = returnReduce(tgtPatch.size(), sumOp<label>());

    if (srcTotalSize == 0)
    {
        DebugInfo<< "AMI: no source faces present - no addressing constructed"
            << endl;

        return false;
    }

    // Philip disable
    // Info<< indent
    //     << "AMI: Creating addressing and weights between "
    //     << srcTotalSize << " source faces and "
    //     << tgtTotalSize << " target faces"
    //     << endl;

    singlePatchProc_ = calcDistribution(srcPatch, tgtPatch);

    if (useGlobalPolyPatch_)
    {
        if (debug)
        {
            Info<< indent
                << "AMI: using globalPolyPatch" << endl;
        }

        singlePatchProc_ = Pstream::myProcNo();
    }

    if (debug)
    {
        Info<< "AMIInterpolationS4F:" << nl
            << "    singlePatchProc:" << singlePatchProc_ << nl
            << endl;
    }

    return true;
}


void Foam::AMIInterpolationS4F::reset
(
    autoPtr<mapDistribute>&& srcToTgtMap,
    autoPtr<mapDistribute>&& tgtToSrcMap,
    labelListList&& srcAddress,
    scalarListList&& srcWeights,
    labelListList&& tgtAddress,
    scalarListList&& tgtWeights
)
{
    DebugInFunction<< endl;

    srcAddress_.transfer(srcAddress);
    srcWeights_.transfer(srcWeights);
    tgtAddress_.transfer(tgtAddress);
    tgtWeights_.transfer(tgtWeights);

    // Reset the sums of the weights
    srcWeightsSum_.setSize(srcWeights_.size());
    forAll(srcWeights_, facei)
    {
        srcWeightsSum_[facei] = sum(srcWeights_[facei]);
    }

    tgtWeightsSum_.setSize(tgtWeights_.size());
    forAll(tgtWeights_, facei)
    {
        tgtWeightsSum_[facei] = sum(tgtWeights_[facei]);
    }

    srcMapPtr_ = std::move(srcToTgtMap);
    tgtMapPtr_ = std::move(tgtToSrcMap);

    upToDate_ = true;
}


void Foam::AMIInterpolationS4F::append
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch
)
{
    addProfiling(ami, "AMIInterpolationS4F::append");

    // Create a new interpolation
    auto newPtr = clone();
    newPtr->calculate(srcPatch, tgtPatch);

    // If parallel then combine the mapDistribution and re-index
    if (distributed())
    {
        labelListList& srcSubMap = srcMapPtr_->subMap();
        labelListList& srcConstructMap = srcMapPtr_->constructMap();

        labelListList& tgtSubMap = tgtMapPtr_->subMap();
        labelListList& tgtConstructMap = tgtMapPtr_->constructMap();

        labelListList& newSrcSubMap = newPtr->srcMapPtr_->subMap();
        labelListList& newSrcConstructMap = newPtr->srcMapPtr_->constructMap();

        labelListList& newTgtSubMap = newPtr->tgtMapPtr_->subMap();
        labelListList& newTgtConstructMap = newPtr->tgtMapPtr_->constructMap();

        // Re-calculate the source indices
        {
            labelList mapMap(0), newMapMap(0);
            forAll(srcSubMap, proci)
            {
                mapMap.append
                (
                    identity
                    (
                        srcConstructMap[proci].size(),
                        (mapMap.size() + newMapMap.size())
                    )
                );
                newMapMap.append
                (
                    identity
                    (
                        newSrcConstructMap[proci].size(),
                        (mapMap.size() + newMapMap.size())
                    )
                );
            }

            forAll(srcSubMap, proci)
            {
                forAll(srcConstructMap[proci], srci)
                {
                    srcConstructMap[proci][srci] =
                        mapMap[srcConstructMap[proci][srci]];
                }
            }

            forAll(srcSubMap, proci)
            {
                forAll(newSrcConstructMap[proci], srci)
                {
                    newSrcConstructMap[proci][srci] =
                        newMapMap[newSrcConstructMap[proci][srci]];
                }
            }

            forAll(tgtAddress_, tgti)
            {
                forAll(tgtAddress_[tgti], tgtj)
                {
                    tgtAddress_[tgti][tgtj] = mapMap[tgtAddress_[tgti][tgtj]];
                }
            }

            forAll(newPtr->tgtAddress_, tgti)
            {
                forAll(newPtr->tgtAddress_[tgti], tgtj)
                {
                    newPtr->tgtAddress_[tgti][tgtj] =
                        newMapMap[newPtr->tgtAddress_[tgti][tgtj]];
                }
            }
        }

        // Re-calculate the target indices
        {
            labelList mapMap(0), newMapMap(0);
            forAll(srcSubMap, proci)
            {
                mapMap.append
                (
                    identity
                    (
                        tgtConstructMap[proci].size(),
                        (mapMap.size() + newMapMap.size())
                    )
                );
                newMapMap.append
                (
                    identity
                    (
                        newTgtConstructMap[proci].size(),
                        (mapMap.size() + newMapMap.size())
                    )
                );
            }

            forAll(srcSubMap, proci)
            {
                forAll(tgtConstructMap[proci], tgti)
                {
                    tgtConstructMap[proci][tgti] =
                        mapMap[tgtConstructMap[proci][tgti]];
                }
            }

            forAll(srcSubMap, proci)
            {
                forAll(newTgtConstructMap[proci], tgti)
                {
                    newTgtConstructMap[proci][tgti] =
                        newMapMap[newTgtConstructMap[proci][tgti]];
                }
            }

            forAll(srcAddress_, srci)
            {
                forAll(srcAddress_[srci], srcj)
                {
                    srcAddress_[srci][srcj] =
                        mapMap[srcAddress_[srci][srcj]];
                }
            }

            forAll(newPtr->srcAddress_, srci)
            {
                forAll(newPtr->srcAddress_[srci], srcj)
                {
                    newPtr->srcAddress_[srci][srcj] =
                        newMapMap[newPtr->srcAddress_[srci][srcj]];
                }
            }
        }

        // Sum the construction sizes
        srcMapPtr_->constructSize() += newPtr->srcMapPtr_->constructSize();
        tgtMapPtr_->constructSize() += newPtr->tgtMapPtr_->constructSize();

        // Combine the maps
        forAll(srcSubMap, proci)
        {
            srcSubMap[proci].append(newSrcSubMap[proci]);
            srcConstructMap[proci].append(newSrcConstructMap[proci]);

            tgtSubMap[proci].append(newTgtSubMap[proci]);
            tgtConstructMap[proci].append(newTgtConstructMap[proci]);
        }
    }

    // Combine new and current source data
    forAll(srcMagSf_, srcFacei)
    {
        srcAddress_[srcFacei].append(newPtr->srcAddress()[srcFacei]);
        srcWeights_[srcFacei].append(newPtr->srcWeights()[srcFacei]);
        srcWeightsSum_[srcFacei] += newPtr->srcWeightsSum()[srcFacei];
    }

    // Combine new and current target data
    forAll(tgtMagSf_, tgtFacei)
    {
        tgtAddress_[tgtFacei].append(newPtr->tgtAddress()[tgtFacei]);
        tgtWeights_[tgtFacei].append(newPtr->tgtWeights()[tgtFacei]);
        tgtWeightsSum_[tgtFacei] += newPtr->tgtWeightsSum()[tgtFacei];
    }
}


void Foam::AMIInterpolationS4F::normaliseWeights
(
    const bool conformal,
    const bool output
)
{
    normaliseWeights
    (
        srcMagSf_,
        "source",
        srcAddress_,
        srcWeights_,
        srcWeightsSum_,
        conformal,
        output,
        lowWeightCorrection_
    );

    normaliseWeights
    (
        tgtMagSf_,
        "target",
        tgtAddress_,
        tgtWeights_,
        tgtWeightsSum_,
        conformal,
        output,
        lowWeightCorrection_
    );
}


Foam::label Foam::AMIInterpolationS4F::srcPointFace
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch,
    const vector& n,
    const label tgtFacei,
    point& tgtPoint
)
const
{
    const pointField& srcPoints = srcPatch.points();

    // Source face addresses that intersect target face tgtFacei
    const labelList& addr = tgtAddress_[tgtFacei];

    pointHit nearest;
    nearest.setDistance(GREAT);
    label nearestFacei = -1;

    for (const label srcFacei : addr)
    {
        const face& f = srcPatch[srcFacei];

        pointHit ray =
            f.ray(tgtPoint, n, srcPoints, intersection::algorithm::VISIBLE);

        if (ray.hit())
        {
            tgtPoint = ray.rawPoint();
            return srcFacei;
        }
        else if (ray.distance() < nearest.distance())
        {

            nearest = ray;
            nearestFacei = srcFacei;
        }
    }

    if (nearest.hit() || nearest.eligibleMiss())
    {
        tgtPoint = nearest.rawPoint();
        return nearestFacei;
    }

    return -1;
}


Foam::label Foam::AMIInterpolationS4F::tgtPointFace
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch,
    const vector& n,
    const label srcFacei,
    point& srcPoint
)
const
{
    const pointField& tgtPoints = tgtPatch.points();

    pointHit nearest;
    nearest.setDistance(GREAT);
    label nearestFacei = -1;

    // Target face addresses that intersect source face srcFacei
    const labelList& addr = srcAddress_[srcFacei];

    for (const label tgtFacei : addr)
    {
        const face& f = tgtPatch[tgtFacei];

        pointHit ray =
            f.ray(srcPoint, n, tgtPoints, intersection::algorithm::VISIBLE);

        if (ray.hit())
        {
            srcPoint = ray.rawPoint();
            return tgtFacei;
        }
        const pointHit near = f.nearestPoint(srcPoint, tgtPoints);

        if (near.distance() < nearest.distance())
        {
            nearest = near;
            nearestFacei = tgtFacei;
        }
    }
    if (nearest.hit() || nearest.eligibleMiss())
    {

        srcPoint = nearest.rawPoint();
        return nearestFacei;
    }

    return -1;
}


bool Foam::AMIInterpolationS4F::checkSymmetricWeights(const bool log) const
{
    if (Pstream::parRun() && (singlePatchProc_ == -1))
    {
        Log << "Checks only valid for serial running (currently)" << endl;

        return true;
    }

    bool symmetricSrc = true;

    Log << "    Checking for missing src face in tgt lists" << nl;

    forAll(srcAddress_, srcFacei)
    {
        const labelList& tgtIds = srcAddress_[srcFacei];
        for (const label tgtFacei : tgtIds)
        {
            if (!tgtAddress_[tgtFacei].found(srcFacei))
            {
                symmetricSrc = false;

                Log << "       srcFacei:" << srcFacei
                    << " not found in tgtToSrc list for tgtFacei:"
                    << tgtFacei << nl;
            }
        }
    }

    if (symmetricSrc)
    {
        Log << "    - symmetric" << endl;
    }

    bool symmetricTgt = true;

    Log << "    Checking for missing tgt face in src lists" << nl;

    forAll(tgtAddress_, tgtFacei)
    {
        const labelList& srcIds = tgtAddress_[tgtFacei];
        for (const label srcFacei : srcIds)
        {
            if (!srcAddress_[srcFacei].found(tgtFacei))
            {
                symmetricTgt = false;

                Log << "       tgtFacei:" << tgtFacei
                    << " not found in srcToTgt list for srcFacei:"
                    << srcFacei << nl;
            }
        }
    }

    if (symmetricTgt)
    {
        Log << "    - symmetric" << endl;
    }

    return symmetricSrc && symmetricTgt;
}


void Foam::AMIInterpolationS4F::writeFaceConnectivity
(
    const standAlonePatch& srcPatch,
    const standAlonePatch& tgtPatch,
    const labelListList& srcAddress
)
const
{
    OFstream os("faceConnectivity" + Foam::name(Pstream::myProcNo()) + ".obj");

    label pti = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];

        for (const label tgtPti : addr)
        {
            const point& tgtPt = tgtPatch.faceCentres()[tgtPti];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << pti << " " << pti + 1 << endl;

            pti += 2;
        }
    }
}


void Foam::AMIInterpolationS4F::write(Ostream& os) const
{
    os.writeEntry("AMIMethod", type());

    if (reverseTarget_)
    {
        os.writeEntry("flipNormals", reverseTarget_);
    }

    if (lowWeightCorrection_ > 0)
    {
        os.writeEntry("lowWeightCorrection", lowWeightCorrection_);
    }
}

#endif // OPENFOAM_COM

// ************************************************************************* //
