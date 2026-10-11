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

#include "quadraticBend.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(quadraticBend, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        quadraticBend,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::quadraticBend::quadraticBend
(
    const dictionary& dict
)
:
    immersedBodyMotion(dict),
    amplitude_(dict.get<scalar>("amplitude")),
    period_(dict.getCheck<scalar>("period", scalarMinMax::ge(SMALL))),
    phase_(dict.getOrDefault<scalar>("phase", 0)),
    yMin_(dict.get<scalar>("yMin")),
    yMax_(dict.get<scalar>("yMax")),
    direction_
    (
        normalised(dict.getOrDefault<vector>("direction", vector(1, 0, 0)))
    )
{
    if (yMax_ <= yMin_ || mag(direction_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "yMax must be greater than yMin and the direction must not be "
            << "zero" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::quadraticBend::points
(
    const pointField& points0,
    const scalar t
) const
{
    const scalar omega = constant::mathematical::twoPi/period_;
    const scalar s = Foam::sin(omega*t + phase_);

    tmp<pointField> tpts(new pointField(points0));
    pointField& pts = tpts.ref();
    forAll(pts, i)
    {
        const scalar xi =
            min(max((points0[i].y() - yMin_)/(yMax_ - yMin_), 0), 1);
        pts[i] += amplitude_*sqr(xi)*s*direction_;
    }
    return tpts;
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::quadraticBend::velocity
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
Foam::immersedBodyMotions::quadraticBend::pointVelocities
(
    const pointField& points0,
    const scalar t
) const
{
    const scalar omega = constant::mathematical::twoPi/period_;
    const scalar c = omega*Foam::cos(omega*t + phase_);

    tmp<vectorField> tU(new vectorField(points0.size()));
    vectorField& U = tU.ref();
    forAll(U, i)
    {
        const scalar xi =
            min(max((points0[i].y() - yMin_)/(yMax_ - yMin_), 0), 1);
        U[i] = amplitude_*sqr(xi)*c*direction_;
    }
    return tU;
}


// ************************************************************************* //
