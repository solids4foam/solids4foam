/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "mixedChamber.H"
#include "mathematicalConstants.H"

#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(mixedChamber, 0);
    addToRunTimeSelectionTable(circulationComponent, mixedChamber, dictionary);
}

Foam::mixedChamber::mixedChamber
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    V0_(ctx.read(dict, "V0", circulationContext::quantity::volume)),
    E1_(0),
    E2_(0),
    P0_(0),
    K_(0),
    Vinit_(V0_),
    sineSquared_(true),
    offset_(0),
    duration_(0),
    T_(ctx.period()),
    Tes_(0),
    Ted_(0),
    Eshift_(0)
{
    if (dict.isDict("active"))
    {
        const dictionary& a = dict.subDict("active");
        E1_ = ctx.read(a, "E1", circulationContext::quantity::elastance, 0);
        E2_ = ctx.read(a, "E2", circulationContext::quantity::elastanceQuadratic, 0);
    }

    if (dict.isDict("passive"))
    {
        const dictionary& p = dict.subDict("passive");
        P0_ = ctx.read(p, "P0", circulationContext::quantity::pressure, 0);
        K_ = ctx.read(p, "K", circulationContext::quantity::inverseVolume, 0);
    }

    const dictionary& act = dict.subDict("activation");
    const word type = act.get<word>("type");
    T_ = act.getOrDefault<scalar>("T", ctx.period());

    if (type == "sineSquared")
    {
        sineSquared_ = true;
        offset_ = act.getOrDefault<scalar>("offset", 0);
        duration_ = act.get<scalar>("duration");
        if (duration_ <= 0 || duration_ > T_) fail("activation duration must be in (0, T]");
    }
    else if (type == "halfCosine")
    {
        sineSquared_ = false;
        Tes_ = act.get<scalar>("Tes");
        Ted_ = act.get<scalar>("Ted");
        Eshift_ = act.getOrDefault<scalar>("Eshift", 0);
        if (Tes_ <= 0 || Ted_ <= Tes_ || Ted_ > T_) fail("need 0 < Tes < Ted <= T");
    }
    else
    {
        fail("activation type must be sineSquared or halfCosine, not " + type);
    }

    if (dict.found("initialVolume"))
    {
        Vinit_ = ctx.read(dict, "initialVolume", circulationContext::quantity::volume);
    }
}


Foam::scalar Foam::mixedChamber::activation(const scalar t) const
{
    using constant::mathematical::pi;

    if (sineSquared_)
    {
        // fmod keeps the sign of its first argument: shift into [0, T)
        scalar phase = std::fmod(t - offset_, T_);
        if (phase < 0) phase += T_;

        if (phase >= 0 && phase <= duration_)
        {
            const scalar s = Foam::sin(pi*phase/duration_);
            return s*s;
        }
        return 0;
    }

    const scalar shifted = t + (1 - Eshift_)*T_;
    const scalar ti = shifted - T_*std::floor(shifted/T_);

    if (ti <= Tes_) return 0.5*(1 - Foam::cos(pi*ti/Tes_));
    if (ti <= Ted_) return 0.5*(1 + Foam::cos(pi*(ti - Tes_)/(Ted_ - Tes_)));
    return 0;
}


Foam::scalar Foam::mixedChamber::pressure(const scalar t, const scalar V) const
{
    const scalar dV = V - V0_;

    return
        activation(t)*(E1_*dV + E2_*dV*dV)
      + P0_*(Foam::exp(K_*dV) - 1);
}


void Foam::mixedChamber::parameters(DynamicList<parameter>& p) const
{
    p.append({"V0", circulationContext::quantity::volume, V0_});
    p.append({"E1", circulationContext::quantity::elastance, E1_});
    p.append({"E2", circulationContext::quantity::elastanceQuadratic, E2_});
    p.append({"P0", circulationContext::quantity::pressure, P0_});
    p.append({"K", circulationContext::quantity::inverseVolume, K_});
}
