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
    Collects the primitive components, connections and terminal aliases
    while the composite components expand, and resolves terminal names.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationBuilder.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationBuilder::circulationBuilder(const circulationContext& ctx)
:
    ctx_(ctx)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circulationBuilder::add(autoPtr<circulationComponent> component)
{
    const word name = component->name();

    if (!names_.insert(name))
    {
        FatalErrorInFunction
            << "Component name '" << name << "' is used twice"
            << exit(FatalError);
    }

    if (component->kind() == circulationComponent::role::composite)
    {
        component->expand(*this);
    }
    else
    {
        primitives_.append(component.release());
    }
}


void Foam::circulationBuilder::connect(const word& from, const word& to)
{
    links_.append(Pair<word>(from, to));
}


void Foam::circulationBuilder::alias
(
    const word& moduleTerminal,
    const word& innerTerminal
)
{
    if (!aliases_.insert(moduleTerminal, innerTerminal))
    {
        FatalErrorInFunction
            << "Terminal '" << moduleTerminal << "' aliased twice"
            << exit(FatalError);
    }
}


Foam::word Foam::circulationBuilder::resolve(const word& terminal) const
{
    word t = terminal;

    for (label hops = 0; hops < 32; ++hops)
    {
        const auto iter = aliases_.cfind(t);
        if (!iter.good())
        {
            return t;
        }

        t = iter.val();
    }

    FatalErrorInFunction
        << "Terminal alias loop starting at '" << terminal << "'"
        << exit(FatalError);

    return t;
}

// ************************************************************************* //
