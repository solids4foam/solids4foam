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

inline const Foam::standAlonePatch& Foam::AMIInterpolationS4F::srcPatch0() const
{
    if (!tsrcPatch0_.valid())
    {
        FatalErrorInFunction
            << "tsrcPatch0Ptr_ not set"
            << abort(FatalError);
    }

    return tsrcPatch0_();
}


inline const Foam::standAlonePatch& Foam::AMIInterpolationS4F::tgtPatch0() const
{

    if (!ttgtPatch0_.valid())
    {
        FatalErrorInFunction
            << "ttgtPatch0Ptr_ not set"
            << abort(FatalError);
    }

    return ttgtPatch0_();
}


inline bool Foam::AMIInterpolationS4F::upToDate() const
{
    return upToDate_;
}


inline bool& Foam::AMIInterpolationS4F::upToDate()
{
    return upToDate_;
}


inline bool Foam::AMIInterpolationS4F::distributed() const
{
    return singlePatchProc_ == -1;
}


inline bool Foam::AMIInterpolationS4F::requireMatch() const
{
    return requireMatch_ && lowWeightCorrection_ < 0;
}


inline bool Foam::AMIInterpolationS4F::setRequireMatch(const bool flag)
{
    requireMatch_ = flag;
    return requireMatch_;
}


inline bool Foam::AMIInterpolationS4F::reverseTarget() const
{
    return reverseTarget_;
}


inline Foam::scalar Foam::AMIInterpolationS4F::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


inline bool Foam::AMIInterpolationS4F::applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


inline Foam::label Foam::AMIInterpolationS4F::singlePatchProc() const
{
    return singlePatchProc_;
}


inline const Foam::List<Foam::scalar>& Foam::AMIInterpolationS4F::srcMagSf() const
{
    return srcMagSf_;
}


inline Foam::List<Foam::scalar>& Foam::AMIInterpolationS4F::srcMagSf()
{
    return srcMagSf_;
}


inline const Foam::labelListList& Foam::AMIInterpolationS4F::srcAddress() const
{
    return srcAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolationS4F::srcAddress()
{
    return srcAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolationS4F::srcWeights() const
{
    return srcWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolationS4F::srcWeights()
{
    return srcWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolationS4F::srcWeightsSum() const
{
    return srcWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolationS4F::srcWeightsSum()
{
    return srcWeightsSum_;
}


inline const Foam::pointListList& Foam::AMIInterpolationS4F::srcCentroids() const
{
    return srcCentroids_;
}


inline Foam::pointListList& Foam::AMIInterpolationS4F::srcCentroids()
{
    return srcCentroids_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolationS4F::srcMap() const
{
    return *srcMapPtr_;
}


inline const Foam::List<Foam::scalar>& Foam::AMIInterpolationS4F::tgtMagSf() const
{
    return tgtMagSf_;
}


inline Foam::List<Foam::scalar>& Foam::AMIInterpolationS4F::tgtMagSf()
{
    return tgtMagSf_;
}


inline const Foam::labelListList& Foam::AMIInterpolationS4F::tgtAddress() const
{
    return tgtAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolationS4F::tgtAddress()
{
    return tgtAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolationS4F::tgtWeights() const
{
    return tgtWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolationS4F::tgtWeights()
{
    return tgtWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolationS4F::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolationS4F::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolationS4F::tgtMap() const
{
    return *tgtMapPtr_;
}


// ************************************************************************* //
