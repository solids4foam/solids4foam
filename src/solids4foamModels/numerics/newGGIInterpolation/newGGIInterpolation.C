/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.

Contributor:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "newGGIInterpolationTemplate.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::clearOut()
{
    deleteDemandDrivenData(masterAddrPtr_);
    deleteDemandDrivenData(masterWeightsPtr_);
    deleteDemandDrivenData(slaveAddrPtr_);
    deleteDemandDrivenData(slaveWeightsPtr_);

    deleteDemandDrivenData(uncoveredMasterAddrPtr_);
    deleteDemandDrivenData(uncoveredSlaveAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class MasterPatch, class SlavePatch>
newGGIInterpolation<MasterPatch, SlavePatch>::newGGIInterpolation
(
    const MasterPatch& masterPatch,
    const SlavePatch&  slavePatch,
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep,
    const scalar masterNonOverlapFaceTol,
    const scalar slaveNonOverlapFaceTol,
    const bool rescaleGGIWeightingFactors,
    const quickReject reject,
    const boundBox& regionOfInterest
)
:
    masterPatch_(masterPatch),
    slavePatch_(slavePatch),
    forwardT_(forwardT),
    reverseT_(reverseT),
    forwardSep_(forwardSep),
    masterNonOverlapFaceTol_(masterNonOverlapFaceTol),
    slaveNonOverlapFaceTol_(slaveNonOverlapFaceTol),
    rescaleGGIWeightingFactors_(rescaleGGIWeightingFactors),
    reject_(reject),
    regionOfInterest_(regionOfInterest),
    masterAddrPtr_(NULL),
    masterWeightsPtr_(NULL),
    slaveAddrPtr_(NULL),
    slaveWeightsPtr_(NULL),
    uncoveredMasterAddrPtr_(NULL),
    uncoveredSlaveAddrPtr_(NULL)
{
    // Check size of transform.  They should be equal to slave patch size
    // if the transform is not constant
    if (forwardT_.size() > 1 || reverseT_.size() > 1)
    {
        if
        (
            forwardT_.size() != slavePatch_.size()
         || reverseT_.size() != masterPatch_.size()
        )
        {
            FatalErrorIn
            (
                "newGGIInterpolation<MasterPatch, SlavePatch>::"
                "newGGIInterpolation"
            )   << "Incorrectly defined transform: forwardT: "
                << forwardT_.size() << " patch: " << slavePatch_.size()
                << " reverseT: " << reverseT_.size()
                << " patch: " << masterPatch_.size()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
newGGIInterpolation<MasterPatch, SlavePatch>::~newGGIInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const labelListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterAddr() const
{
    if (!masterAddrPtr_)
    {
        calcAddressing();
    }

    return *masterAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::masterWeights() const
{
    if (!masterWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveAddr() const
{
    if (!slaveAddrPtr_)
    {
        calcAddressing();
    }

    return *slaveAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
newGGIInterpolation<MasterPatch, SlavePatch>::slaveWeights() const
{
    if (!slaveWeightsPtr_)
    {
        calcAddressing();
    }

    return *slaveWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelList&
newGGIInterpolation<MasterPatch, SlavePatch>::uncoveredMasterFaces() const
{
    if (!uncoveredMasterAddrPtr_)
    {
        calcAddressing();
    }

    return *uncoveredMasterAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelList&
newGGIInterpolation<MasterPatch, SlavePatch>::uncoveredSlaveFaces() const
{
    if (!uncoveredSlaveAddrPtr_)
    {
        calcAddressing();
    }

    return *uncoveredSlaveAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
bool newGGIInterpolation<MasterPatch, SlavePatch>::movePoints
(
    const tensorField& forwardT,
    const tensorField& reverseT,
    const vectorField& forwardSep
)
{
    this->forwardT_ = forwardT;
    this->reverseT_ = reverseT;
    this->forwardSep_ = forwardSep;

    clearOut();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newGGIInterpolationPolygonIntersection.C"
#   include "newGGIInterpolationQuickRejectTests.C"
#   include "newGGIInterpolationWeights.C"
#   include "newGGIInterpolate.C"

// ************************************************************************* //
