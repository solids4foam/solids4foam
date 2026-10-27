/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "CRL.H"
#include "circulationBuilder.H"
#include "compliance.H"
#include "resistor.H"
#include "inductor.H"

namespace Foam
{
    defineTypeNameAndDebug(CRL, 0);
    addToRunTimeSelectionTable(circulationComponent, CRL, dictionary);
}

Foam::CRL::CRL(const word& name, const dictionary& dict, const circulationContext& ctx)
:
    circulationComponent(name, dict, ctx)
{}

void Foam::CRL::expand(circulationBuilder& b) const
{
    b.add(autoPtr<circulationComponent>(new compliance(name_ + ".C", dict_, ctx_)));
    b.add(autoPtr<circulationComponent>(new resistor(name_ + ".R", dict_, ctx_)));
    b.add(autoPtr<circulationComponent>(new inductor(name_ + ".L", dict_, ctx_)));

    b.connect(name_ + ".C.out", name_ + ".R.in");
    b.connect(name_ + ".R.out", name_ + ".L.in");

    b.alias(name_ + ".in", name_ + ".C.in");
    b.alias(name_ + ".out", name_ + ".L.out");
}
