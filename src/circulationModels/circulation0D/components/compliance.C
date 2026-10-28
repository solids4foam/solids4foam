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
    Reads the parameters of the compliance, a node that stores volume, and
    reports them for the set-up check.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "compliance.H"

namespace Foam
{
    defineTypeNameAndDebug(compliance, 0);
    addToRunTimeSelectionTable(circulationComponent, compliance, dictionary);
}

Foam::compliance::compliance
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    C_(ctx.read(dict, "C", circulationContext::quantity::compliance)),
    V0_(ctx.read(dict, "V0", circulationContext::quantity::volume, 0)),
    P0_(ctx.read(dict, "P0", circulationContext::quantity::pressure, 0)),
    Vinit_(V0_)
{
    if (C_ <= 0)
    {
        fail("C must be positive");
    }

    if (dict.found("initialVolume"))
    {
        Vinit_ =
            ctx.read
            (
                dict,
                "initialVolume",
                circulationContext::quantity::volume
            );
    }
    else if (dict.found("initialPressure"))
    {
        Vinit_ =
            V0_
          + C_
           *(
                ctx.read
                (
                    dict,
                    "initialPressure",
                    circulationContext::quantity::pressure
                )
              - P0_
            );
    }
}


void Foam::compliance::parameters(DynamicList<parameter>& p) const
{
    p.append({"C", circulationContext::quantity::compliance, C_});
    p.append({"V0", circulationContext::quantity::volume, V0_});
    p.append({"P0", circulationContext::quantity::pressure, P0_});
}

// ************************************************************************* //
