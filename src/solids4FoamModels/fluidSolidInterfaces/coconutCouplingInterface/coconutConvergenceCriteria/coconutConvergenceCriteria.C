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

#include "coconutConvergenceCriteria.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coconutConvergenceCriteria, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coconutConvergenceCriteria::coconutConvergenceCriteria
(
    const scalar tolerance,
    const int maxIter
)
:
    tolerance_(tolerance),
    iter_(0),
    maxIter_(maxIter),
    initialNorm_(0.0),
    lastNorm_(0.0),
    isInitialNormSet_(false)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coconutConvergenceCriteria::initializeSolutionStep()
{
    iter_ = 0;
    initialNorm_ = 0.0;
    lastNorm_ = 0.0;
    isInitialNormSet_ = false;
}


void coconutConvergenceCriteria::update(const vectorField& r)
{
    // Increment the iteration
    iter_++;

    // Calculate the norm of the residual vector (assume order is 2)
    // lastNorm_ = sum(mag(r)); // order 1
    lastNorm_ = sum(magSqr(r)); // order 2

    if (!isInitialNormSet_)
    {
        initialNorm_ = lastNorm_;
        isInitialNormSet_ = true;

        // Check if the initial norm is too small
        if (initialNorm_ < SMALL)
        {
            FatalErrorInFunction
                << "Initial norm is too small" << exit(FatalError);
        }
    }
}


bool coconutConvergenceCriteria::isSatisfied() const
{
    if (!isInitialNormSet_)
    {
        return false;
    }
    else
    {
        return (lastNorm_/initialNorm_) < tolerance_ || iter_ > maxIter_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
