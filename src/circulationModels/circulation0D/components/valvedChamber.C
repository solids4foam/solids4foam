/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Expands a valvedChamber into its inlet valve, elastance chamber and
    outlet valve in series.

Author
    David Kelly, UCD.

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

    b.add
    (
        autoPtr<circulationComponent>
        (
            new orificeValve(mv, dict_.subDict("inletValve"), ctx_)
        )
    );
    b.add
    (
        autoPtr<circulationComponent>
        (
            new elastanceChamber(ch, dict_.subDict("chamber"), ctx_)
        )
    );
    b.add
    (
        autoPtr<circulationComponent>
        (
            new orificeValve(av, dict_.subDict("outletValve"), ctx_)
        )
    );

    b.connect(mv + ".out", ch + ".in");
    b.connect(ch + ".out", av + ".in");

    b.alias(name_ + ".in", mv + ".in");
    b.alias(name_ + ".out", av + ".out");
}

// ************************************************************************* //
