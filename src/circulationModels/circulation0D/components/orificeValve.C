/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "orificeValve.H"

namespace Foam
{
    defineTypeNameAndDebug(orificeValve, 0);
    addToRunTimeSelectionTable(circulationComponent, orificeValve, dictionary);
}

Foam::orificeValve::orificeValve
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    CQ_
    (
        dict.found("Rhat")
      ? 1/Foam::sqrt(ctx.read(dict, "Rhat", circulationContext::quantity::quadraticResistance))
      : ctx.read(dict, "CQ", circulationContext::quantity::valveCoefficient)
    ),
    Rs_(ctx.read(dict, "R", circulationContext::quantity::resistance, 0)),
    eps_(ctx.read(dict, "epsilon", circulationContext::quantity::pressure, 0)),
    reverseFraction_(dict.getOrDefault<scalar>("reverseFraction", 0))
{
    if (CQ_ <= 0) fail("CQ (or Rhat) must be positive");
    if (Rs_ < 0) fail("R must not be negative");
    if (Rs_ > 0 && dict.getOrDefault<scalar>("epsilon", 0) > 0) fail("R and epsilon cannot be combined");
    if (eps_ < 0) fail("epsilon must not be negative");
    if (reverseFraction_ < 0 || reverseFraction_ > 1)
    {
        fail("reverseFraction must be in [0,1]");
    }
}


Foam::scalar Foam::orificeValve::flow(const scalar, const scalar dP) const
{
    if (eps_ > 0)
    {
        // Smooth diode: sqrt of a strictly positive max(dP,0) and a sigmoid gate
        const scalar dPpos = 0.5*(dP + Foam::sqrt(dP*dP + eps_*eps_));
        const scalar a = Foam::min(Foam::max(dP/eps_, scalar(-50)), scalar(50));
        const scalar gate = 1.0/(1.0 + Foam::exp(-a));
        return CQ_*Foam::sqrt(Foam::max(dPpos, VSMALL))*gate;
    }

    if (dP >= 0)
    {
        if (Rs_ > 0)
        {
            // Positive root of Rhat Q^2 + R Q = dP, Rhat = 1/CQ^2
            const scalar Rhat = 1/(CQ_*CQ_);
            return (-Rs_ + Foam::sqrt(Rs_*Rs_ + 4*Rhat*dP))/(2*Rhat);
        }
        return CQ_*Foam::sqrt(dP);
    }

    return -reverseFraction_*CQ_*Foam::sqrt(-dP);
}


Foam::scalar Foam::orificeValve::dFlowdDP(const scalar t, const scalar dP) const
{
    // Used by the nodal solve (and the finite-difference Jacobian). The slope
    // of the square root is infinite where the valve opens; the floor keeps
    // it finite without hiding the stiffness of the opening.
    if (eps_ > 0)
    {
        const scalar h = Foam::max(1e-3, 1e-6*Foam::mag(dP));
        return (flow(t, dP + h) - flow(t, dP - h))/(2*h);
    }

    const scalar floorP = 1e-9;     // Pa

    if (dP > 0)
    {
        if (Rs_ > 0)
        {
            // Q = (-R + sqrt(R^2 + 4 Rhat dP))/(2 Rhat), Rhat = 1/CQ^2
            const scalar Rhat = 1/(CQ_*CQ_);
            return 1/Foam::sqrt(Rs_*Rs_ + 4*Rhat*dP);
        }
        return CQ_/(2*Foam::sqrt(Foam::max(dP, floorP)));
    }

    return
        reverseFraction_ > 0
      ? reverseFraction_*CQ_/(2*Foam::sqrt(Foam::max(-dP, floorP)))
      : scalar(0);
}


Foam::scalar Foam::orificeValve::dPfromQ
(
    const scalar CQ,
    const scalar Q,
    const scalar fr
)
{
    if (Q >= 0)
    {
        const scalar s = Q/CQ;
        return s*s;
    }

    const scalar frc = Foam::min(Foam::max(fr, scalar(1e-3)), scalar(1));
    const scalar s = -Q/(frc*CQ);
    return -s*s;
}


void Foam::orificeValve::parameters(DynamicList<parameter>& p) const
{
    p.append({"CQ", circulationContext::quantity::valveCoefficient, CQ_});
    if (Rs_ > 0) p.append({"R", circulationContext::quantity::resistance, Rs_});
}
