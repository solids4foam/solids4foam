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
    Expands a coupledDomain into an externalPort for each of its coupled
    patches, with the valve law of the case between the lumped node and the
    patch.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "coupledDomain.H"
#include "circulationBuilder.H"
#include "externalPort.H"

namespace Foam
{
    defineTypeNameAndDebug(coupledDomain, 0);
    addToRunTimeSelectionTable(circulationComponent, coupledDomain, dictionary);
}

Foam::coupledDomain::coupledDomain
(
    const word& name,
    const dictionary& dict,
    const circulationContext& ctx
)
:
    circulationComponent(name, dict, ctx)
{}

void Foam::coupledDomain::expand(circulationBuilder& b) const
{
    const word provider = dict_.getOrDefault<word>("provider", "external");

    b.add
    (
        autoPtr<circulationComponent>
        (
            new externalPort
            (
                name_ + ".inlet",
                dict_.subDict("inlet"),
                ctx_,
                "inlet",
                provider
            )
        )
    );
    b.add
    (
        autoPtr<circulationComponent>
        (
            new externalPort
            (
                name_ + ".outlet",
                dict_.subDict("outlet"),
                ctx_,
                "outlet",
                provider
            )
        )
    );

    b.alias(name_ + ".in", name_ + ".inlet.in");
    b.alias(name_ + ".out", name_ + ".outlet.out");
}

// ************************************************************************* //
