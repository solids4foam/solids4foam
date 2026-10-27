/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "pressureSource.H"

namespace Foam
{
    defineTypeNameAndDebug(pressureSource, 0);
    addToRunTimeSelectionTable(circulationComponent, pressureSource, dictionary);
}

Foam::pressureSource::pressureSource
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    wave_(dict, ctx, circulationContext::quantity::pressure)
{}
