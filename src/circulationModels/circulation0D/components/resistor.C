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
    Reads the resistance of the resistor and reports it for the set-up
    check.

Author
    David Kelly, UCD.

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
    if (R_ <= 0)
    {
        fail("R must be positive");
    }
}


void Foam::resistor::parameters(DynamicList<parameter>& p) const
{
    p.append({"R", circulationContext::quantity::resistance, R_});
}

// ************************************************************************* //
