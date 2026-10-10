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

#include "solidBodyRotation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(solidBodyRotation, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        solidBodyRotation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::solidBodyRotation::solidBodyRotation
(
    const dictionary& dict
)
:
    immersedBodyMotion(dict),
    origin_(dict.get<point>("origin")),
    axis_(normalised(dict.get<vector>("axis"))),
    omega_(dict.get<scalar>("omega")),
    startTime_(dict.getOrDefault<scalar>("startTime", 0))
{
    if (mag(axis_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "The axis must not be zero" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::solidBodyRotation::points
(
    const pointField& points0,
    const scalar t
) const
{
    // Rodrigues' rotation formula
    const scalar theta = omega_*(t - startTime_);
    const scalar c = Foam::cos(theta);
    const scalar s = Foam::sin(theta);
    const tensor R
    (
        c*I + s*tensor
        (
            0, -axis_.z(), axis_.y(),
            axis_.z(), 0, -axis_.x(),
            -axis_.y(), axis_.x(), 0
        )
      + (1 - c)*sqr(axis_)
    );

    return origin_ + (R & (points0 - origin_));
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::solidBodyRotation::velocity
(
    const pointField& x,
    const scalar t
) const
{
    return (omega_*axis_) ^ (x - origin_);
}


// ************************************************************************* //
