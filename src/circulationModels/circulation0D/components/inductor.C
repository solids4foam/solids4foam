/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "inductor.H"

namespace Foam
{
    defineTypeNameAndDebug(inductor, 0);
    addToRunTimeSelectionTable(circulationComponent, inductor, dictionary);
}

Foam::inductor::inductor
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    L_(ctx.read(dict, "L", circulationContext::quantity::inertance)),
    q0_(ctx.read(dict, "initialFlow", circulationContext::quantity::flow, 0))
{
    if (L_ <= 0) fail("L must be positive");
}


void Foam::inductor::parameters(DynamicList<parameter>& p) const
{
    p.append({"L", circulationContext::quantity::inertance, L_});
}
