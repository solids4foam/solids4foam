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
    Registers the ground component, a node held at a fixed pressure.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "ground.H"

namespace Foam
{
    defineTypeNameAndDebug(ground, 0);
    addToRunTimeSelectionTable(circulationComponent, ground, dictionary);
}

Foam::ground::ground
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx),
    P_(ctx.read(dict, "P", circulationContext::quantity::pressure, 0))
{}

// ************************************************************************* //
