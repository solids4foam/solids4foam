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

#include "profiling.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::AMIInterpolationS4F::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolationS4F::interpolateToTarget");

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

    if (distributed())
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


template<class Type, class CombineOp>
void Foam::AMIInterpolationS4F::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolationS4F::interpolateToSource");

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
                << " but supplied default field size is not equal to source "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    source patch   = " << srcAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(srcAddress_.size());

    if (distributed())
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


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(srcAddress_.size(), Zero);

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(tgtAddress_.size(), Zero);

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolationS4F::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


// ************************************************************************* //
