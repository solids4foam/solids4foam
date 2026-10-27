/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "valvedChamber.H"
#include "circulationBuilder.H"
#include "elastanceChamber.H"
#include "orificeValve.H"

namespace Foam
{
    defineTypeNameAndDebug(valvedChamber, 0);
    addToRunTimeSelectionTable(circulationComponent, valvedChamber, dictionary);
}

Foam::valvedChamber::valvedChamber
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx)
{}

void Foam::valvedChamber::expand(circulationBuilder& b) const
{
    const word mv = name_ + ".inletValve";
    const word ch = name_ + ".chamber";
    const word av = name_ + ".outletValve";

    b.add(autoPtr<circulationComponent>(new orificeValve(mv, dict_.subDict("inletValve"), ctx_)));
    b.add(autoPtr<circulationComponent>(new elastanceChamber(ch, dict_.subDict("chamber"), ctx_)));
    b.add(autoPtr<circulationComponent>(new orificeValve(av, dict_.subDict("outletValve"), ctx_)));

    b.connect(mv + ".out", ch + ".in");
    b.connect(ch + ".out", av + ".in");

    b.alias(name_ + ".in", mv + ".in");
    b.alias(name_ + ".out", av + ".out");
}
