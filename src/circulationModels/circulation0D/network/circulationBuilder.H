/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

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
