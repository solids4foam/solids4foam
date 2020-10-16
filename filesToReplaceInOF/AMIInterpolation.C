/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation.H"
#include "AMIMethod.H"
#include "meshTools.H"
#include "mapDistribute.H"
#include "flipOp.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
const Foam::Enum
<
    typename Foam::AMIInterpolation<SourcePatch, TargetPatch>::
    interpolationMethod
>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolationMethodNames_
({
    { interpolationMethod::imDirect, "directAMI" },
    { interpolationMethod::imMapNearest, "mapNearestAMI" },
    { interpolationMethod::imFaceAreaWeight, "faceAreaWeightAMI" },
    { interpolationMethod::imPartialFaceAreaWeight, "partialFaceAreaWeightAMI" }
});

template<class SourcePatch, class TargetPatch>
bool Foam::AMIInterpolation<SourcePatch, TargetPatch>::cacheIntersections_ =
    false;

template<class SourcePatch, class TargetPatch>
template<class Patch>
Foam::tmp<Foam::scalarField>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::patchMagSf
(
    const Patch& patch,
    const faceAreaIntersect::triangulationMode triMode
)
{
    tmp<scalarField> tResult(new scalarField(patch.size(), Zero));
    scalarField& result = tResult.ref();

    const pointField& patchPoints = patch.localPoints();

    faceList patchFaceTris;

    forAll(result, patchFacei)
    {
        faceAreaIntersect::triangulate
        (
            patch.localFaces()[patchFacei],
            patchPoints,
            triMode,
            patchFaceTris
        );

        forAll(patchFaceTris, i)
        {
            result[patchFacei] +=
                triPointRef
                (
                    patchPoints[patchFaceTris[i][0]],
                    patchPoints[patchFaceTris[i][1]],
                    patchPoints[patchFaceTris[i][2]]
                ).mag();
        }
    }

    return tResult;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    addProfiling(ami, "AMIInterpolation::projectPointsToSurface");

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


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
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
    addProfiling(ami, "AMIInterpolation::normaliseWeights");

    // Normalise the weights
    wghtSum.setSize(wght.size(), 0.0);
    label nLowWeight = 0;

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

            if (t < lowWeightTol)
            {
                nLowWeight++;
            }
        }
        else
        {
            wghtSum[facei] = 0;
        }
    }


    if (output)
    {
        const label nFace = returnReduce(wght.size(), sumOp<label>());

        if (nFace)
        {
            Info<< indent
                << "AMI: Patch " << patchName
                << " sum(weights)"
                << " min:" << gMin(wghtSum)
                << " max:" << gMax(wghtSum)
                << " average:" << gAverage(wghtSum) << nl;

            const label nLow = returnReduce(nLowWeight, sumOp<label>());

            if (nLow)
            {
                Info<< indent
                    << "AMI: Patch " << patchName
                    << " identified " << nLow
                    << " faces with weights less than " << lowWeightTol
                    << endl;
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::agglomerate
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
    addProfiling(ami, "AMIInterpolation::agglomerate");

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
    if (targetMapPtr.valid())
    {
        const mapDistribute& map = targetMapPtr();

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
        // - a constructMap : these are from 'transferred-date' to slots

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
                label elemi = elems[i];
                label coarseElemi = targetRestrictAddressing[elemi];

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


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::constructFromSurface
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (surfPtr.valid())
    {
        // Create new patches for source and target
        pointField srcPoints = srcPatch.points();
        SourcePatch srcPatch0
        (
            SubList<face>
            (
                srcPatch,
                srcPatch.size(),
                0
            ),
            srcPoints
        );

        if (debug)
        {
            OFstream os("amiSrcPoints.obj");
            for (const point& pt : srcPoints)
            {
                meshTools::writeOBJ(os, pt);
            }
        }

        pointField tgtPoints = tgtPatch.points();
        TargetPatch tgtPatch0
        (
            SubList<face>
            (
                tgtPatch,
                tgtPatch.size(),
                0
            ),
            tgtPoints
        );

        if (debug)
        {
            OFstream os("amiTgtPoints.obj");
            for (const point& pt : tgtPoints)
            {
                meshTools::writeOBJ(os, pt);
            }
        }


        // Map source and target patches onto projection surface
        projectPointsToSurface(surfPtr(), srcPoints);
        projectPointsToSurface(surfPtr(), tgtPoints);


        // Calculate AMI interpolation
        update(srcPatch0, tgtPatch0);
    }
    else
    {
        update(srcPatch, tgtPatch);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const bool useGlobalPolyPatch
)
:
    methodName_(interpolationMethodNames_[method]),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{
    update(srcPatch, tgtPatch);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const word& methodName,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const bool useGlobalPolyPatch
)
:
    methodName_(methodName),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{
    update(srcPatch, tgtPatch);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const bool useGlobalPolyPatch
)
:
    methodName_(interpolationMethodNames_[method]),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{
    constructFromSurface(srcPatch, tgtPatch, surfPtr);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const word& methodName,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const bool useGlobalPolyPatch
)
:
    methodName_(methodName),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    useGlobalPolyPatch_(useGlobalPolyPatch)
{
    constructFromSurface(srcPatch, tgtPatch, surfPtr);
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const AMIInterpolation<SourcePatch, TargetPatch>& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing
)
:
    methodName_(fineAMI.methodName_),
    reverseTarget_(fineAMI.reverseTarget_),
    requireMatch_(fineAMI.requireMatch_),
    singlePatchProc_(fineAMI.singlePatchProc_),
    lowWeightCorrection_(-1.0),
    srcMagSf_(),
    srcAddress_(),
    srcWeights_(),
    srcWeightsSum_(),
    tgtMagSf_(),
    tgtAddress_(),
    tgtWeights_(),
    tgtWeightsSum_(),
    triMode_(fineAMI.triMode_),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    useGlobalPolyPatch_(false)
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


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::~AMIInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::update
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
)
{
    addProfiling(ami, "AMIInterpolation::update");

    label srcTotalSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtTotalSize = returnReduce(tgtPatch.size(), sumOp<label>());

    if (srcTotalSize == 0)
    {
        DebugInfo<< "AMI: no source faces present - no addressing constructed"
            << endl;

        return;
    }

    Info<< indent
        << "AMI: Creating addressing and weights between "
        << srcTotalSize << " source faces and "
        << tgtTotalSize << " target faces"
        << endl;

    // Calculate if patches present on multiple processors
    singlePatchProc_ = calcDistribution(srcPatch, tgtPatch);

    if (useGlobalPolyPatch_)
    {
        Info<< indent
            << "AMI: using globalPolyPatch" << endl;
        singlePatchProc_ = Pstream::myProcNo();
    }

    if (singlePatchProc_ == -1)
    {
        // Convert local addressing to global addressing
        globalIndex globalSrcFaces(srcPatch.size());
        globalIndex globalTgtFaces(tgtPatch.size());

        // Create processor map of overlapping faces. This map gets
        // (possibly remote) faces from the tgtPatch such that they (together)
        // cover all of the srcPatch
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcPatch, tgtPatch);
        const mapDistribute& map = mapPtr();

        // Create new target patch that fully encompasses source patch

        // Faces and points
        faceList newTgtFaces;
        pointField newTgtPoints;

        // Original faces from tgtPatch (in globalIndexing since might be
        // remote)
        labelList tgtFaceIDs;
        distributeAndMergePatches
        (
            map,
            tgtPatch,
            globalTgtFaces,
            newTgtFaces,
            newTgtPoints,
            tgtFaceIDs
        );

        const TargetPatch
            newTgtPatch
            (
                SubList<face>
                (
                    newTgtFaces,
                    newTgtFaces.size()
                ),
                newTgtPoints
            );

        // Calculate AMI interpolation
        autoPtr<AMIMethod<SourcePatch, TargetPatch>> AMIPtr
        (
            AMIMethod<SourcePatch, TargetPatch>::New
            (
                methodName_,
                srcPatch,
                newTgtPatch,
                triMode_,
                reverseTarget_,
                requireMatch_ && (lowWeightCorrection_ < 0)
            )
        );

        AMIPtr->calculate
        (
            srcAddress_,
            srcWeights_,
            tgtAddress_,
            tgtWeights_
        );


        // Note: using patch face areas calculated by the AMI method
        // - TODO: move into the calculate or normalise method?
        AMIPtr->setMagSf(tgtPatch, map, srcMagSf_, tgtMagSf_);


        // Now
        // ~~~
        //  srcAddress_ :   per srcPatch face a list of the newTgtPatch (not
        //                  tgtPatch) faces it overlaps
        //  tgtAddress_ :   per newTgtPatch (not tgtPatch) face a list of the
        //                  srcPatch faces it overlaps

        if (debug)
        {
            writeFaceConnectivity(srcPatch, newTgtPatch, srcAddress_);
        }


        // Rework newTgtPatch indices into globalIndices of tgtPatch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        for (labelList& addressing : srcAddress_)
        {
            for (label& addr : addressing)
            {
                addr = tgtFaceIDs[addr];
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
            tgtPatch.size(),
            map.constructMap(),
            false,                      // has flip
            map.subMap(),
            false,                      // has flip
            tgtAddress_,
            ListOps::appendEqOp<label>(),
            flipOp(),                   // flip operation
            labelList()
        );

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtWeights_,
            ListOps::appendEqOp<scalar>(),
            flipOp(),
            scalarList()
        );

        // weights normalisation
        AMIPtr->normaliseWeights(true, *this);

        // Cache maps and reset addresses
        List<Map<label>> cMap;
        srcMapPtr_.reset(new mapDistribute(globalSrcFaces, tgtAddress_, cMap));
        tgtMapPtr_.reset(new mapDistribute(globalTgtFaces, srcAddress_, cMap));
    }
    else
    {
        // Calculate AMI interpolation
        autoPtr<AMIMethod<SourcePatch, TargetPatch>> AMIPtr
        (
            AMIMethod<SourcePatch, TargetPatch>::New
            (
                methodName_,
                srcPatch,
                tgtPatch,
                triMode_,
                reverseTarget_,
                requireMatch_ && (lowWeightCorrection_ < 0)
            )
        );

        AMIPtr->calculate
        (
            srcAddress_,
            srcWeights_,
            tgtAddress_,
            tgtWeights_
        );

        srcMagSf_.transfer(AMIPtr->srcMagSf());
        tgtMagSf_.transfer(AMIPtr->tgtMagSf());

        AMIPtr->normaliseWeights(true, *this);
    }

    if (debug)
    {
        Info<< "AMIInterpolation : Constructed addressing and weights" << nl
            << "    triMode        :"
            << faceAreaIntersect::triangulationModeNames_[triMode_] << nl
            << "    singlePatchProc:" << singlePatchProc_ << nl
            << "    srcMagSf       :" << gSum(srcMagSf_) << nl
            << "    tgtMagSf       :" << gSum(tgtMagSf_) << nl
            << endl;
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::append
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
)
{
    addProfiling(ami, "AMIInterpolation::append");

    // Create a new interpolation
    autoPtr<AMIInterpolation<SourcePatch, TargetPatch>> newPtr
    (
        new AMIInterpolation<SourcePatch, TargetPatch>
        (
            srcPatch,
            tgtPatch,
            triMode_,
            requireMatch_,
            methodName_,
            lowWeightCorrection_,
            reverseTarget_
        )
    );

    // If parallel then combine the mapDistribution and re-index
    if (singlePatchProc_ == -1)
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
                    tgtAddress_[tgti][tgtj] =
                        mapMap[tgtAddress_[tgti][tgtj]];
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


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
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


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToTarget");

    if (fld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != tgtAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    target patch   = " << tgtAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(tgtAddress_.size());

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, facei)
        {
            if (tgtWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = tgtAddress_[facei];
                const scalarList& weights = tgtWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            if (tgtWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = tgtAddress_[facei];
                const scalarList& weights = tgtWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToSource");

    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != srcAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    source patch   = " << srcAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(srcAddress_.size());

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, facei)
        {
            if (srcWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = srcAddress_[facei];
                const scalarList& weights = srcWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            if (srcWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = srcAddress_[facei];
                const scalarList& weights = srcWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            srcAddress_.size(),
            Zero
        )
    );

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            tgtAddress_.size(),
            Zero
        )
    );

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
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

        pointHit ray = f.ray(tgtPoint, n, srcPoints);

        if (ray.hit())
        {
            // tgtPoint = ray.rawPoint();
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
        // tgtPoint = nearest.rawPoint();
        return nearestFacei;
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
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

        pointHit ray = f.ray(srcPoint, n, tgtPoints);

        if (ray.hit() || ray.eligibleMiss())
        {
            // srcPoint = ray.rawPoint();
            return tgtFacei;
        }
        else if (ray.distance() < nearest.distance())
        {
            nearest = ray;
            nearestFacei = tgtFacei;
        }
    }

    if (nearest.hit() || nearest.eligibleMiss())
    {
        // srcPoint = nearest.rawPoint();
        return nearestFacei;
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeFaceConnectivity
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
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


// ************************************************************************* //
