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
    Run-time selection of circulation components by their "type", and the
    default numerical members of the base class, which stop with an error
    for a component whose role does not provide them.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationComponent.H"
#include "circulationModule.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(circulationComponent, 0);
    defineRunTimeSelectionTable(circulationComponent, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationComponent::circulationComponent
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    name_(name),
    dict_(dict),
    ctx_(ctx)
{}


Foam::autoPtr<Foam::circulationComponent> Foam::circulationComponent::New
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
{
    const word type(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        // Not a C++ component: a module (built-in preset or the case's own)?
        if (const dictionary* module = ctx.findModule(type))
        {
            return autoPtr<circulationComponent>
            (
                new circulationModule(name, dict, ctx, *module)
            );
        }

        FatalIOErrorInFunction(dict)
            << "Unknown circulation component type '" << type
            << "' for component '" << name << "'." << nl
            << "Components:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << "Modules (presets and this case's \"modules\"):" << nl
            << ctx.moduleNames() << nl
            << exit(FatalIOError);
    }

    return autoPtr<circulationComponent>(ctorPtr(name, dict, ctx));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circulationComponent::fail(const string& message) const
{
    FatalIOErrorInFunction(dict_)
        << "Component '" << name_ << "': " << message << nl
        << exit(FatalIOError);

    std::abort(); // not reached; satisfies [[noreturn]]
}


#define CIRC_NOT_IMPLEMENTED(what)                                            \
    fail(string("'") + type() + "' does not provide " what);                  \
    return 0

Foam::scalar Foam::circulationComponent::flow(const scalar, const scalar) const
{
    CIRC_NOT_IMPLEMENTED("flow()");
}

Foam::scalar Foam::circulationComponent::dFlowdDP
(
    const scalar,
    const scalar
) const
{
    CIRC_NOT_IMPLEMENTED("dFlowdDP()");
}

Foam::scalar Foam::circulationComponent::inertance() const
{
    CIRC_NOT_IMPLEMENTED("inertance()");
}

Foam::scalar Foam::circulationComponent::fixedFlow(const scalar) const
{
    CIRC_NOT_IMPLEMENTED("fixedFlow()");
}

Foam::scalar Foam::circulationComponent::pressure
(
    const scalar,
    const scalar
) const
{
    CIRC_NOT_IMPLEMENTED("pressure()");
}

Foam::scalar Foam::circulationComponent::initialVolume() const
{
    CIRC_NOT_IMPLEMENTED("initialVolume()");
}

Foam::scalar Foam::circulationComponent::fixedPressure(const scalar) const
{
    CIRC_NOT_IMPLEMENTED("fixedPressure()");
}

Foam::scalar Foam::circulationComponent::injectedFlow(const scalar) const
{
    CIRC_NOT_IMPLEMENTED("injectedFlow()");
}

Foam::scalar Foam::circulationComponent::portFlow(const scalar) const
{
    CIRC_NOT_IMPLEMENTED("portFlow()");
}

Foam::scalar Foam::circulationComponent::patchPressure
(
    const scalar,
    const scalar
) const
{
    CIRC_NOT_IMPLEMENTED("patchPressure()");
}

#undef CIRC_NOT_IMPLEMENTED

// ************************************************************************* //
