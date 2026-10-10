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

#include "sinusoidalTranslation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(sinusoidalTranslation, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        sinusoidalTranslation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::immersedBodyMotions::sinusoidalTranslation::omega() const
{
    return constant::mathematical::twoPi/period_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::sinusoidalTranslation::sinusoidalTranslation
(
    const dictionary& dict
)
:
    immersedBodyMotion(dict),
    amplitude_(dict.get<scalar>("amplitude")),
    period_(dict.getCheck<scalar>("period", scalarMinMax::ge(SMALL))),
    phase_(dict.getOrDefault<scalar>("phase", 0)),
    direction_
    (
        normalised(dict.getOrDefault<vector>("direction", vector(1, 0, 0)))
    ),
    offset_(dict.getOrDefault<vector>("offset", Zero))
{
    if (mag(direction_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "The direction must not be zero" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::sinusoidalTranslation::points
(
    const pointField& points0,
    const scalar t
) const
{
    const vector displacement =
        amplitude_*Foam::sin(omega()*t + phase_)*direction_ + offset_;

    return points0 + displacement;
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::sinusoidalTranslation::velocity
(
    const pointField& x,
    const scalar t
) const
{
    return tmp<vectorField>::New
    (
        x.size(),
        amplitude_*omega()*Foam::cos(omega()*t + phase_)*direction_
    );
}


// ************************************************************************* //
