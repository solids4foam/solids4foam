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

#include "fsiDriven.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(fsiDriven, 0);
    addToRunTimeSelectionTable(immersedBodyMotion, fsiDriven, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::fsiDriven::fsiDriven(const dictionary& dict)
:
    immersedBodyMotion(dict),
    points_(),
    velocities_(),
    configuration_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedBodyMotions::fsiDriven::setPointsAndVelocities
(
    const pointField& points,
    const vectorField& velocities
)
{
    if (points.size() != velocities.size())
    {
        FatalErrorInFunction
            << "The numbers of points (" << points.size()
            << ") and velocities (" << velocities.size() << ") differ"
            << abort(FatalError);
    }

    points_ = points;
    velocities_ = velocities;
    ++configuration_;
}


Foam::tmp<Foam::pointField> Foam::immersedBodyMotions::fsiDriven::points
(
    const pointField& points0,
    const scalar t
) const
{
    if (points_.empty())
    {
        return tmp<pointField>::New(points0);
    }

    if (points_.size() != points0.size())
    {
        FatalErrorInFunction
            << "The number of points set by the fluid-solid interface ("
            << points_.size() << ") differs from the number of surface "
            << "points (" << points0.size() << ")" << abort(FatalError);
    }

    return tmp<pointField>::New(points_);
}


Foam::tmp<Foam::vectorField> Foam::immersedBodyMotions::fsiDriven::velocity
(
    const pointField& x,
    const scalar t
) const
{
    FatalErrorInFunction
        << "The velocity of the deforming motion " << typeName
        << " is only defined at the surface points" << abort(FatalError);

    return tmp<vectorField>::New(x.size(), Zero);
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::fsiDriven::pointVelocities
(
    const pointField& points0,
    const scalar t
) const
{
    if (velocities_.empty())
    {
        return tmp<vectorField>::New(points0.size(), Zero);
    }

    return tmp<vectorField>::New(velocities_);
}


// ************************************************************************* //
