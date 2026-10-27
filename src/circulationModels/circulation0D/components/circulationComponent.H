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
    Foam::circulationComponent

Description
    Base class of everything that can appear in the "components" dictionary
    of a circulation network, selected at run time by its "type".

    A component has two terminals, "in" and "out". Flow is positive from in
    to out and pressure is shared by all terminals joined by a connection (a
    "net"). The role of a component says what it contributes to the network:

    \table
        role          | meaning                                    | examples
        composite     | expands into smaller components            |
                        CR, CRL, valvedChamber, coupledDomain
        branch        | flow is an algebraic function of dP        |
                        resistor, orificeValve
        inertia       | flow is a state, dq/dt = dP/L              | inductor
        fixedFlow     | flow prescribed in time (in -> out)        | flowSource
        storage       | its net holds a volume state V, P = P(V,t) |
                        compliance, elastanceChamber
        fixedPressure | its net pressure is prescribed             |
                        ground, pressureSource
        portSource    | flow injected into its net from outside    |
                        externalPort
    \endtable

    storage, fixedPressure and portSource are one-node elements: their two
    terminals are the same node.

    Sign conventions: dP = Pin - Pout is positive
    for forward flow, and dV/dt of a storage net is the sum of the flows into
    it minus the flows out of it.

    All quantities handed to and returned from the numerical members are SI.

Author
    David Kelly, UCD.

SourceFiles
    circulationComponent.C

\*---------------------------------------------------------------------------*/

#ifndef circulationComponent_H
#define circulationComponent_H

#include "circulationContext.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "DynamicList.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "typeInfo.H"

namespace Foam
{

class circulationBuilder;

class circulationComponent
{
public:

    enum class role
    {
        composite,
        branch,
        inertia,
        fixedFlow,
        storage,
        fixedPressure,
        portSource
    };


protected:

    const word name_;
    const dictionary dict_;
    const circulationContext& ctx_;

    //- Fatal error with the component name and its dictionary location
    [[noreturn]] void fail(const string& message) const;


public:

    //- A physical parameter of a component, for the check mode
    struct parameter
    {
        word name;
        circulationContext::quantity quantity;
        scalar value;       //!< SI
    };

    TypeName("circulationComponent");

    declareRunTimeSelectionTable
    (
        autoPtr,
        circulationComponent,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const circulationContext& ctx
        ),
        (name, dict, ctx)
    );

    //- Construct from components
    circulationComponent
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    //- Select by dict.get<word>("type")
    static autoPtr<circulationComponent> New
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual ~circulationComponent() = default;


    // Member Functions

        const word& name() const { return name_; }

        virtual role kind() const = 0;

        //- Terminals "in" and "out" are one node
        bool singleNode() const
        {
            const role r = kind();
            return
                r == role::storage
             || r == role::fixedPressure
             || r == role::portSource;
        }

        //- Report the physical parameters (SI), used by the setup check
        virtual void parameters(DynamicList<parameter>&) const {}

        // composite

            //- Add sub-components, connections and terminal aliases
            virtual void expand(circulationBuilder&) const {}

        // branch

            //- Flow in -> out for dP = Pin - Pout [m3/s]
            virtual scalar flow(const scalar t, const scalar dP) const;

            //- d flow / d dP, used only by the nodal solve [m3/s/Pa]
            virtual scalar dFlowdDP(const scalar t, const scalar dP) const;

            //- True if flow() is linear in dP, so one Newton step is exact
            virtual bool linear() const { return false; }

        // inertia

            //- Inertance L [Pa s2/m3]
            virtual scalar inertance() const;

            //- Initial flow [m3/s]
            virtual scalar initialFlow() const { return 0; }

        // fixedFlow

            //- Prescribed flow in -> out [m3/s]
            virtual scalar fixedFlow(const scalar t) const;

        // storage

            //- Pressure for volume V [Pa]
            virtual scalar pressure(const scalar t, const scalar V) const;

            virtual scalar initialVolume() const;

        // fixedPressure

            virtual scalar fixedPressure(const scalar t) const;

        // portSource

            //- Flow injected into the node [m3/s]
            virtual scalar injectedFlow(const scalar t) const;

            //- Flow through the port in its natural direction [m3/s]
            virtual scalar portFlow(const scalar t) const;

            //- Pressure to impose at the coupled patch for node pressure
            //  Pnode and the port flow at time t [Pa]
            virtual scalar patchPressure
            (
                const scalar Pnode,
                const scalar t
            ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
