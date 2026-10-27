/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "compliance.H"

namespace Foam
{
    defineTypeNameAndDebug(compliance, 0);
    addToRunTimeSelectionTable(circulationComponent, compliance, dictionary);
}

Foam::compliance::compliance
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    C_(ctx.read(dict, "C", circulationContext::quantity::compliance)),
    V0_(ctx.read(dict, "V0", circulationContext::quantity::volume, 0)),
    P0_(ctx.read(dict, "P0", circulationContext::quantity::pressure, 0)),
    Vinit_(V0_)
{
    if (C_ <= 0) fail("C must be positive");

    if (dict.found("initialVolume"))
    {
        Vinit_ = ctx.read(dict, "initialVolume", circulationContext::quantity::volume);
    }
    else if (dict.found("initialPressure"))
    {
        Vinit_ =
            V0_
          + C_*(ctx.read(dict, "initialPressure", circulationContext::quantity::pressure) - P0_);
    }
}


void Foam::compliance::parameters(DynamicList<parameter>& p) const
{
    p.append({"C", circulationContext::quantity::compliance, C_});
    p.append({"V0", circulationContext::quantity::volume, V0_});
    p.append({"P0", circulationContext::quantity::pressure, P0_});
}
