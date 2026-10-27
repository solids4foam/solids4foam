/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "elastanceChamber.H"
#include "mathematicalConstants.H"
#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(elastanceChamber, 0);
    addToRunTimeSelectionTable(circulationComponent, elastanceChamber, dictionary);
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
    if (Emin_ < 0 || Emax_ < Emin_) fail("need 0 <= Emin <= Emax");
    if (Tes_ <= 0 || Ted_ <= Tes_) fail("need 0 < Tes < Ted");
    if (Ted_ > T_) fail("Ted exceeds the cycle length T");

    if (dict.found("initialVolume"))
    {
        Vinit_ = ctx.read(dict, "initialVolume", circulationContext::quantity::volume);
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
