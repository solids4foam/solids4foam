/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "CR.H"
#include "circulationBuilder.H"
#include "compliance.H"
#include "resistor.H"

namespace Foam
{
    defineTypeNameAndDebug(CR, 0);
    addToRunTimeSelectionTable(circulationComponent, CR, dictionary);
}

Foam::CR::CR(const word& name, const dictionary& dict, const circulationContext& ctx)
:
    circulationComponent(name, dict, ctx)
{}

void Foam::CR::expand(circulationBuilder& b) const
{
    b.add(autoPtr<circulationComponent>(new compliance(name_ + ".C", dict_, ctx_)));
    b.add(autoPtr<circulationComponent>(new resistor(name_ + ".R", dict_, ctx_)));

    b.connect(name_ + ".C.out", name_ + ".R.in");

    b.alias(name_ + ".in", name_ + ".C.in");
    b.alias(name_ + ".out", name_ + ".R.out");
}
