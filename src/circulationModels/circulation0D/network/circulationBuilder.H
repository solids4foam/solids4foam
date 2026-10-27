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

Class
    Foam::circulationBuilder

Description
    Collects the primitive components, connections and terminal aliases while
    a circulation network is being read. Composite components (CR, CRL, a
    valvedChamber, a coupledDomain) call back into it from expand() to add
    their parts; a terminal of a composite is an alias for a terminal of one
    of its parts, so that "SAS.in" can be connected exactly like the
    terminal of a primitive.

    Terminals are named "<component>.in" and "<component>.out".

Author
    David Kelly, UCD.

SourceFiles
    circulationBuilder.C

\*---------------------------------------------------------------------------*/

#ifndef circulationBuilder_H
#define circulationBuilder_H

#include "circulationComponent.H"
#include "PtrList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "HashSet.H"
#include "Pair.H"

namespace Foam
{

class circulationBuilder
{
    const circulationContext& ctx_;

    PtrList<circulationComponent> primitives_;

    //- Every name used so far, composite or primitive
    HashSet<word> names_;

    DynamicList<Pair<word>> links_;

    HashTable<word> aliases_;

public:

    explicit circulationBuilder(const circulationContext& ctx);

    const circulationContext& ctx() const { return ctx_; }

    //- Add a component; a composite is expanded immediately
    void add(autoPtr<circulationComponent> component);

    //- Join two terminals ("A.out", "B.in")
    void connect(const word& from, const word& to);

    //- Declare terminal "<composite>.in" to be the inner terminal
    void alias(const word& moduleTerminal, const word& innerTerminal);

    //- Follow aliases down to a primitive terminal name
    word resolve(const word& terminal) const;

    PtrList<circulationComponent>& primitives() { return primitives_; }

    const DynamicList<Pair<word>>& links() const { return links_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
