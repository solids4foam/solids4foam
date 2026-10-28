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
    Reads the parameters of the inductor, which represents blood inertia,
    and reports them for the set-up check.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "inductor.H"

namespace Foam
{
    defineTypeNameAndDebug(inductor, 0);
    addToRunTimeSelectionTable(circulationComponent, inductor, dictionary);
}

Foam::inductor::inductor
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    L_(ctx.read(dict, "L", circulationContext::quantity::inertance)),
    q0_(ctx.read(dict, "initialFlow", circulationContext::quantity::flow, 0))
{
    if (L_ <= 0)
    {
        fail("L must be positive");
    }
}


void Foam::inductor::parameters(DynamicList<parameter>& p) const
{
    p.append({"L", circulationContext::quantity::inertance, L_});
}

// ************************************************************************* //
