/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "ground.H"

namespace Foam
{
    defineTypeNameAndDebug(ground, 0);
    addToRunTimeSelectionTable(circulationComponent, ground, dictionary);
}

Foam::ground::ground
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    P_(ctx.read(dict, "P", circulationContext::quantity::pressure, 0))
{}
