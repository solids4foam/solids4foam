/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "resistor.H"

namespace Foam
{
    defineTypeNameAndDebug(resistor, 0);
    addToRunTimeSelectionTable(circulationComponent, resistor, dictionary);
}

Foam::resistor::resistor
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    R_(ctx.read(dict, "R", circulationContext::quantity::resistance))
{
    if (R_ <= 0) fail("R must be positive");
}


void Foam::resistor::parameters(DynamicList<parameter>& p) const
{
    p.append({"R", circulationContext::quantity::resistance, R_});
}
