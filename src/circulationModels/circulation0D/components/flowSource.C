/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "flowSource.H"

namespace Foam
{
    defineTypeNameAndDebug(flowSource, 0);
    addToRunTimeSelectionTable(circulationComponent, flowSource, dictionary);
}

Foam::flowSource::flowSource
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    wave_(dict, ctx, circulationContext::quantity::flow)
{}
