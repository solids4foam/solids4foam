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
    Instantiates a sub-circuit module: substitutes its parameters into the
    module dictionary and expands it into components, connections and
    terminal aliases.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationModule.H"
#include "circulationBuilder.H"
#include "primitiveEntry.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- Copy "in" into "out" replacing every word token "@name" by the tokens of
//  the parameter "name"
void substitute
(
    const dictionary& in,
    const dictionary& params,
    const word& moduleName,
    dictionary& out
)
{
    for (const entry& e : in)
    {
        if (e.isDict())
        {
            dictionary sub;
            substitute(e.dict(), params, moduleName, sub);
            out.add(e.keyword(), sub);
            continue;
        }

        DynamicList<token> toks;
        const ITstream& is = e.stream();

        for (const token& t : is)
        {
            if
            (
                t.isWord()
             && t.wordToken().size() > 1
             && t.wordToken()[0] == '@'
            )
            {
                const word pname(t.wordToken().substr(1), false);
                const entry* pe = params.findEntry(pname, keyType::LITERAL);

                if (!pe)
                {
                    FatalIOErrorInFunction(in)
                        << "Module '" << moduleName << "' refers to the"
                        << " parameter @" << pname << ", which it does not"
                        << " declare in required or defaults"
                        << exit(FatalIOError);
                }

                for (const token& pt : pe->stream())
                {
                    toks.append(pt);
                }
            }
            else
            {
                toks.append(t);
            }
        }

        out.add(new primitiveEntry(e.keyword(), List<token>(toks)));
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationModule::circulationModule
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx,
    const dictionary& definition
)
:
    circulationComponent(name, dict, ctx),
    moduleType_(definition.dictName())
{
    const word modName = definition.dictName();

    // Units the module's own numbers are written in
    if
    (
        definition.found("units")
     && definition.get<word>("units") != ctx.name()
    )
    {
        fail
        (
            "module '" + modName + "' is written in units "
          + definition.get<word>("units") + " but this case uses "
          + ctx.name() + ". Set  units " + definition.get<word>("units")
          + ";  or give every parameter explicitly"
        );
    }

    const wordList required
    (
        definition.getOrDefault<wordList>("required", wordList())
    );
    const dictionary defaults(definition.subOrEmptyDict("defaults"));

    // Parameters: the defaults, overridden by the instance's entries
    dictionary params;
    for (const entry& e : defaults)
    {
        params.add(e.clone(params).release());
    }

    wordHashSet known;
    for (const word& r : required)
    {
        known.insert(r);
    }
    for (const entry& e : defaults)
    {
        known.insert(e.keyword());
    }

    for (const entry& e : dict)
    {
        if (e.keyword() == "type")
        {
            continue;
        }

        if (!known.found(e.keyword()))
        {
            string valid;
            for (const word& k : known.sortedToc())
            {
                valid += " " + k;
            }
            fail
            (
                "module '" + modName + "' has no parameter '" + e.keyword()
              + "'. Parameters:" + valid
            );
        }

        params.remove(e.keyword());
        params.add(e.clone(params).release());
    }

    for (const word& r : required)
    {
        if (!params.found(r))
        {
            fail("module '" + modName + "' needs the parameter '" + r + "'");
        }
    }

    // Everything the module defines, with the parameters in place
    for
    (
        const word& part :
        {word("components"), word("connections"), word("terminals")}
    )
    {
        if (!definition.found(part))
        {
            continue;
        }

        dictionary one, oneOut;
        one.add
        (
            definition.findEntry(part, keyType::LITERAL)->clone(one).release()
        );
        substitute(one, params, modName, oneOut);
        for (const entry& e : oneOut)
        {
            resolved_.add(e.clone(resolved_).release());
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::circulationModule::type() const
{
    return moduleType_;
}


void Foam::circulationModule::expand(circulationBuilder& b) const
{
    // Parts
    if (resolved_.isDict("components"))
    {
        for (const entry& e : resolved_.subDict("components"))
        {
            if (!e.isDict())
            {
                fail
                (
                    "component '" + e.keyword()
                  + "' of the module must be a dictionary"
                );
            }

            b.add
            (
                circulationComponent::New
                (
                    name_ + "." + e.keyword(),
                    e.dict(),
                    ctx_
                )
            );
        }
    }

    // Connections, with names relative to the module: a chain (A B C) joins
    // A.out to B.in and B.out to C.in; "A.in"/"A.out" name a terminal
    if (resolved_.found("connections"))
    {
        const List<wordList> chains(resolved_.lookup("connections"));

        for (const wordList& chain : chains)
        {
            if (chain.size() < 2)
            {
                fail("a module connection needs at least two names");
            }

            for (label i = 0; i + 1 < chain.size(); ++i)
            {
                const word& a = chain[i];
                const word& c = chain[i + 1];
                const bool aExplicit = a.find('.') != std::string::npos;
                const bool cExplicit = c.find('.') != std::string::npos;

                b.connect
                (
                    name_ + "." + (aExplicit ? a : word(a + ".out")),
                    name_ + "." + (cExplicit ? c : word(c + ".in"))
                );
            }
        }
    }

    // Terminals of the module as a whole
    if (resolved_.isDict("terminals"))
    {
        for (const entry& e : resolved_.subDict("terminals"))
        {
            b.alias
            (
                name_ + "." + e.keyword(),
                name_ + "." + e.stream()[0].wordToken()
            );
        }
    }
}

// ************************************************************************* //
