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

#include "customProfileBend.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(customProfileBend, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        customProfileBend,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::customProfileBend::displacedPoints
(
    const pointField& points0,
    const scalar s,
    const bool addPoints
) const
{
    const scalar H = yMax_ - yMin_;

    tmp<pointField> tpts
    (
        addPoints
      ? new pointField(points0)
      : new pointField(points0.size(), Zero)
    );
    pointField& pts = tpts.ref();

    forAll(points0, i)
    {
        const scalar xi = min(max((points0[i].y() - yMin_)/H, 0), 1);
        const scalar g = xi*xi*xi*(10 - 15*xi + 6*xi*xi);
        const scalar dg = 30*xi*xi*(1 - 2*xi + xi*xi);

        pts[i].x() += amplitudeX_*g*s;
        pts[i].y() -= (points0[i].x() - xCenter_)*amplitudeY_/H*dg*s;
    }

    return tpts;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::customProfileBend::customProfileBend
(
    const dictionary& dict
)
:
    immersedBodyMotion(dict),
    amplitudeX_(dict.get<scalar>("amplitudeX")),
    amplitudeY_(dict.get<scalar>("amplitudeY")),
    period_(dict.getCheck<scalar>("period", scalarMinMax::ge(SMALL))),
    phase_(dict.getOrDefault<scalar>("phase", 0)),
    yMin_(dict.get<scalar>("yMin")),
    yMax_(dict.get<scalar>("yMax")),
    xCenter_(dict.get<scalar>("xCenter"))
{
    if (yMax_ <= yMin_)
    {
        FatalIOErrorInFunction(dict)
            << "yMax must be greater than yMin" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::customProfileBend::points
(
    const pointField& points0,
    const scalar t
) const
{
    const scalar omega = constant::mathematical::twoPi/period_;
    return displacedPoints(points0, Foam::sin(omega*t + phase_), true);
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::customProfileBend::velocity
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
Foam::immersedBodyMotions::customProfileBend::pointVelocities
(
    const pointField& points0,
    const scalar t
) const
{
    const scalar omega = constant::mathematical::twoPi/period_;
    return
        displacedPoints
        (
            points0,
            omega*Foam::cos(omega*t + phase_),
            false
        );
}


// ************************************************************************* //
