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

#include "valveMorph.H"
#include "valveSliceAxis.H"
#include "valveTimeLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(valveMorph, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        valveMorph,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBodyMotions::valveMorph::checkPoints
(
    const pointField& points0
) const
{
    if (points0.size() != openPoints_.size())
    {
        FatalErrorInFunction
            << "The open surface " << openSurface_ << " has "
            << openPoints_.size() << " vertices, but the valve surface has "
            << points0.size() << ": the two surfaces must have the same "
            << "vertices in the same order" << exit(FatalError);
    }
}


Foam::scalar Foam::immersedBodyMotions::valveMorph::alpha
(
    const scalar t
) const
{
    const scalar phi = valveTimeLaw::phase(t, period_);

    return
        valveTimeLaw::edge(phi, openWindow_)
       *(1 - valveTimeLaw::edge(phi, closeWindow_));
}


Foam::scalar Foam::immersedBodyMotions::valveMorph::alphaDerivative
(
    const scalar t
) const
{
    const scalar phi = valveTimeLaw::phase(t, period_);

    return
    (
        valveTimeLaw::edgeDerivative(phi, openWindow_)
       *(1 - valveTimeLaw::edge(phi, closeWindow_))
      - valveTimeLaw::edge(phi, openWindow_)
       *valveTimeLaw::edgeDerivative(phi, closeWindow_)
    )/period_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::valveMorph::valveMorph(const dictionary& dict)
:
    immersedBodyMotion(dict),
    openSurface_(dict.get<fileName>("openSurface")),
    openPoints_
    (
        triSurface(valveSliceAxis::surfaceFile(openSurface_)).points()
    ),
    period_(dict.getCheck<scalar>("period", scalarMinMax::ge(SMALL))),
    openWindow_(valveTimeLaw::readWindow(dict, "openWindow")),
    closeWindow_(valveTimeLaw::readWindow(dict, "closeWindow"))
{
    if (openWindow_.second() > closeWindow_.first())
    {
        FatalIOErrorInFunction(dict)
            << "The openWindow " << openWindow_ << " must end before the "
            << "closeWindow " << closeWindow_ << " begins"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::valveMorph::points
(
    const pointField& points0,
    const scalar t
) const
{
    checkPoints(points0);

    return points0 + alpha(t)*(openPoints_ - points0);
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::valveMorph::velocity
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
Foam::immersedBodyMotions::valveMorph::pointVelocities
(
    const pointField& points0,
    const scalar t
) const
{
    checkPoints(points0);

    return alphaDerivative(t)*(openPoints_ - points0);
}


// ************************************************************************* //
