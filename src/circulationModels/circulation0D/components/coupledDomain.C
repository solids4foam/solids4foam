/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "coupledDomain.H"
#include "circulationBuilder.H"
#include "externalPort.H"

namespace Foam
{
    defineTypeNameAndDebug(coupledDomain, 0);
    addToRunTimeSelectionTable(circulationComponent, coupledDomain, dictionary);
}

Foam::coupledDomain::coupledDomain
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx)
{}

void Foam::coupledDomain::expand(circulationBuilder& b) const
{
    const word provider = dict_.getOrDefault<word>("provider", "external");

    b.add
    (
        autoPtr<circulationComponent>
        (
            new externalPort(name_ + ".inlet", dict_.subDict("inlet"), ctx_, "inlet", provider)
        )
    );
    b.add
    (
        autoPtr<circulationComponent>
        (
            new externalPort(name_ + ".outlet", dict_.subDict("outlet"), ctx_, "outlet", provider)
        )
    );

    b.alias(name_ + ".in", name_ + ".inlet.in");
    b.alias(name_ + ".out", name_ + ".outlet.out");
}
