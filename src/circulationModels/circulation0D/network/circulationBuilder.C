/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
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
        if (!iter.good()) return t;
        t = iter.val();
    }

    FatalErrorInFunction
        << "Terminal alias loop starting at '" << terminal << "'"
        << exit(FatalError);

    return t;
}

// ************************************************************************* //
