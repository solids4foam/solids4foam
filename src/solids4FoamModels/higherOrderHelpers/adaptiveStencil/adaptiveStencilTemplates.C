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

#include "DynamicList.H"
#include "cpuTime.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CalcEstimate, class StencilSize, class Enlarge>
Foam::scalarField Foam::adaptiveStencil::adapt
(
    const word& stencilType,
    const label nStencils,
    const CalcEstimate& calcEstimate,
    const StencilSize& stencilSize,
    const Enlarge& enlarge
) const
{
    cpuTime timer;

    // Stencil sizes requested by the reconstruction scheme
    labelList baseSize(nStencils);
    forAll(baseSize, i)
    {
        baseSize[i] = stencilSize(i);
    }

    // The estimate is calculated for all stencils in the first pass, and
    // only for the enlarged stencils afterwards
    boolList active(nStencils, true);
    scalarField estimate(nStencils, 0.0);
    conditioning request(active, estimate);

    // Every stencil grows by extraCellsPerPass in each pass, or reaches its
    // limit, so this number of passes is not exceeded
    const label maxPasses =
        (maxExtraCells_ + extraCellsPerPass_ - 1)/extraCellsPerPass_;

    label nPasses = 0;

    while (true)
    {
        // Parallel communication may be required here, so this is called on
        // all processors, also if there are no local active stencils
        calcEstimate(request);

        // Stencils to be enlarged
        DynamicList<label> poor;

        forAll(active, i)
        {
            if (active[i])
            {
                if
                (
                    estimate[i] > maxConditionNumber_
                 && stencilSize(i) - baseSize[i] < maxExtraCells_
                )
                {
                    poor.append(i);
                }
                else
                {
                    active[i] = false;
                }
            }
        }

        if
        (
            returnReduce(label(poor.size()), sumOp<label>()) == 0
         || nPasses >= maxPasses
        )
        {
            break;
        }

        // Growth requested for each stencil is limited to its remaining
        // allowance, so that maxExtraCells is not exceeded
        labelList oldSize(poor.size());
        labelList nExtra(poor.size());
        forAll(poor, i)
        {
            oldSize[i] = stencilSize(poor[i]);

            const label remaining =
                maxExtraCells_ - (oldSize[i] - baseSize[poor[i]]);

            nExtra[i] = min(extraCellsPerPass_, remaining);
        }

        enlarge(poor, nExtra);

        // Stencils that did not grow have no more candidate cells
        forAll(poor, i)
        {
            if (stencilSize(poor[i]) == oldSize[i])
            {
                active[poor[i]] = false;
            }
        }

        ++nPasses;
    }

    // Report
    label nUsed = 0;
    label nEnlarged = 0;
    label nPoor = 0;
    label maxSize = 0;
    scalar maxEstimate = 0.0;

    forAll(baseSize, i)
    {
        const label size = stencilSize(i);

        if (size > 0)
        {
            ++nUsed;
        }
        if (size > baseSize[i])
        {
            ++nEnlarged;
        }
        if (estimate[i] > maxConditionNumber_)
        {
            ++nPoor;
        }

        maxSize = max(maxSize, size);
        maxEstimate = max(maxEstimate, estimate[i]);
    }

    reduce(nUsed, sumOp<label>());
    reduce(nEnlarged, sumOp<label>());
    reduce(nPoor, sumOp<label>());
    reduce(maxSize, maxOp<label>());
    reduce(maxEstimate, maxOp<scalar>());

    Info<< "Adaptive " << stencilType << " stencils: " << nEnlarged
        << " of " << nUsed << " enlarged in " << nPasses << " passes, "
        << timer.elapsedCpuTime() << " s" << nl
        << "    largest stencil: " << maxSize << " cells" << nl
        << "    largest condition number estimate: " << maxEstimate << nl
        << "    stencils above maxConditionNumber ("
        << maxConditionNumber_ << "): " << nPoor << endl;

    return estimate;
}


// ************************************************************************* //
