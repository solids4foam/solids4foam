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
    Expands the CR vessel into a compliance followed by a resistor.

Author
    David Kelly, UCD.

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

Foam::CR::CR
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx)
{}

void Foam::CR::expand(circulationBuilder& b) const
{
    b.add
    (
        autoPtr<circulationComponent>
        (
            new compliance(name_ + ".C", dict_, ctx_)
        )
    );
    b.add
    (
        autoPtr<circulationComponent>
        (
            new resistor(name_ + ".R", dict_, ctx_)
        )
    );

    b.connect(name_ + ".C.out", name_ + ".R.in");

    b.alias(name_ + ".in", name_ + ".C.in");
    b.alias(name_ + ".out", name_ + ".R.out");
}

// ************************************************************************* //
