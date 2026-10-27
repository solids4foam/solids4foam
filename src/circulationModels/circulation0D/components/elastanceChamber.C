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

Description
    Reads the parameters of the time-varying elastance heart chamber and
    evaluates its half-cosine elastance.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "elastanceChamber.H"
#include "mathematicalConstants.H"
#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(elastanceChamber, 0);
    addToRunTimeSelectionTable
    (
        circulationComponent,
        elastanceChamber,
        dictionary
    );
}

Foam::elastanceChamber::elastanceChamber
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    V0_(ctx.read(dict, "V0", circulationContext::quantity::volume)),
    P0_(ctx.read(dict, "P0", circulationContext::quantity::pressure, 0)),
    Emin_(ctx.read(dict, "Emin", circulationContext::quantity::elastance)),
    Emax_(ctx.read(dict, "Emax", circulationContext::quantity::elastance)),
    Tes_(dict.get<scalar>("Tes")),
    Ted_(dict.get<scalar>("Ted")),
    Eshift_(dict.getOrDefault<scalar>("Eshift", 0)),
    T_(dict.getOrDefault<scalar>("T", ctx.period())),
    Vinit_(V0_)
{
    if (Emin_ < 0 || Emax_ < Emin_)
    {
        fail("need 0 <= Emin <= Emax");
    }

    if (Tes_ <= 0 || Ted_ <= Tes_)
    {
        fail("need 0 < Tes < Ted");
    }

    if (Ted_ > T_)
    {
        fail("Ted exceeds the cycle length T");
    }

    if (dict.found("initialVolume"))
    {
        Vinit_ =
            ctx.read
            (
                dict,
                "initialVolume",
                circulationContext::quantity::volume
            );
    }
}


Foam::scalar Foam::elastanceChamber::elastance(const scalar t) const
{
    using constant::mathematical::pi;

    const scalar shifted = t + (1 - Eshift_)*T_;
    const scalar ti = shifted - T_*std::floor(shifted/T_);   // mod(), >= 0

    scalar f = 0;
    if (ti <= Tes_)
    {
        f = 0.5*(1 - Foam::cos(pi*ti/Tes_));
    }
    else if (ti <= Ted_)
    {
        f = 0.5*(1 + Foam::cos(pi*(ti - Tes_)/(Ted_ - Tes_)));
    }

    return Emin_ + (Emax_ - Emin_)*f;
}


void Foam::elastanceChamber::parameters(DynamicList<parameter>& p) const
{
    p.append({"V0", circulationContext::quantity::volume, V0_});
    p.append({"P0", circulationContext::quantity::pressure, P0_});
    p.append({"Emin", circulationContext::quantity::elastance, Emin_});
    p.append({"Emax", circulationContext::quantity::elastance, Emax_});
}

// ************************************************************************* //
