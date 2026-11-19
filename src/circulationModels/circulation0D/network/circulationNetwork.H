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
    Foam::circulationNetwork

Description
    A lumped-parameter model of the circulation assembled from a dictionary
    and solved with an OpenFOAM ODESolver.

    The dictionary lists named components and how they are connected:

        units       mmHg_mL;
        Tcycle      1.0;

        components
        {
            SAS { type CRL; C 0.08; R 0.003; L 6.2e-5; initialPressure 100; }
            ...
        }

        connections
        (
            (LV SAS SAT SAR SCP SVN RA)        // a chain: A.out -> B.in ...
            (SVN.out RA.in)                    // or with explicit terminals
        );

        solver { type RKDP45; absTol 1e-12; relTol 1e-8; }

        probes ( Qin.out );    // optional: report the pressure at a terminal

    Components that are connected share one pressure (a net). A net holding
    a compliance or chamber has a volume state; every other net is solved
    algebraically for its pressure (Newton, exact in one step when only
    resistors meet there) so that the flows into it balance. Inductors add a
    flow state. The result is an explicit ODE system, dy/dt = f(t, y), for
    the volumes and inductor flows, handed to the ODESolver selected in
    "solver".

    Coupling to a 3D domain: replace a valvedChamber by a coupledDomain. The
    network then exposes its ports (portNames()); the coupling code supplies
    the port flows over each step with setPortForcing(), calls advance(),
    and reads back the pressures to impose with portPressure().

Author
    David Kelly, UCD.

SourceFiles
    circulationNetwork.C

\*---------------------------------------------------------------------------*/

#ifndef circulationNetwork_H
#define circulationNetwork_H

#include "circulationBuilder.H"
#include "externalPort.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "HashTable.H"
#include "stringList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class circulationNetwork Declaration
\*---------------------------------------------------------------------------*/

class circulationNetwork
:
    public ODESystem
{
    //- Role of a component in the network
    typedef circulationComponent::role role;

    // Private Data

        //- Unit system and cardiac period of the case
        circulationContext ctx_;

        //- Primitive components, composites already expanded
        PtrList<circulationComponent> comps_;

        //- Role of each component
        List<role> roles_;

        //- Net of the "in" terminal of each component
        labelList netIn_;

        //- Net of the "out" terminal of each component
        labelList netOut_;

        //- Number of nets
        label nNets_;

        //- Net -> component holding its volume state, -1 if none
        labelList netStorage_;

        //- Net -> component fixing its pressure, -1 if none
        labelList netFixed_;

        //- Net -> index in the algebraic unknowns, -1 if not algebraic
        labelList algIndex_;

        //- Number of algebraic nets
        label nAlg_;

        //- True if every branch touching an algebraic net is linear
        bool algLinear_;

        //- An inductor whose flow is not a state because it carries the whole
        //  flow of a group of algebraic nodes that has no pressure reference
        //  (for example a flow source feeding r -> L -> C: q_L = Q(t)). Its
        //  pressure drop L dq/dt pins the group instead.
        struct forcedInertia
        {
            //- The inductor
            label comp;

            //- Its end inside the group
            label pinnedNet;

            //- Its end outside
            label otherNet;

            //- +1 if the pinned end is "in", else -1
            scalar sign;

            //- The resistive group it carries the flow of
            label group;
        };

        //- The forced inductors
        List<forcedInertia> forced_;

        //- Component -> index in forced_, or -1
        labelList forcedOfComp_;

        //- Algebraic net -> resistive group
        labelList groupOfNet_;

        //- Net whose nodal row is replaced
        boolList pinned_;

        //- Component -> its state index (storage: volume, inertia: flow)
        labelList stateOfComp_;

        //- Number of ODE states
        label nStates_;

        //- Port name -> component index
        HashTable<label> portIndex_;

        //- ODE solver settings, "solver" sub-dictionary with defaults added
        dictionary odeDict_;

        //- The ODE solver
        autoPtr<ODESolver> solver_;

        //- Current time
        scalar t_;

        //- Current state: volumes and inductor flows [SI]
        scalarField y_;

        //- Step-size estimate of the ODE solver
        scalar dxEst_;

        //- Work space; derivatives() is const but warm-starts the nodal solve
        mutable scalarField pAlg_;

        //- Work space for the pressure of every net
        mutable scalarField Pwork_;

        //- Column headings of the output
        wordList columns_;

        //- Nets whose pressure is reported (the optional "probes" list)
        labelList probeNet_;


    // Private Member Functions

        //- Expand the components, classify the nets, number the states and
        //  create the ODE solver
        void build(const dictionary& dict);

        //- Flow imbalance R at every algebraic node, and optionally its
        //  Jacobian, for algebraic-node pressures pa (also fills P at those
        //  nodes)
        void nodalSystem
        (
            const scalar t,
            const scalarField& y,
            scalarField& P,
            const scalarField& pa,
            scalarField& R,
            scalarSquareMatrix* J,
            scalar& flowScale
        ) const;

        //- Pressure of every net at (t, y), solving the algebraic ones
        void netPressures
        (
            const scalar t,
            const scalarField& y,
            scalarField& P
        ) const;

        //- Flow in -> out of component k at (t, y) with net pressures P
        scalar branchFlow
        (
            const label k,
            const scalar t,
            const scalarField& y,
            const scalarField& P
        ) const;

        //- Names of the components attached to a net, for messages
        string netMembers(const label net) const;

        //- Component index of a port; fatal error if there is no such port
        label portLabel(const word& port) const;

        //- Flow in -> out of a forced inductor, from the sources feeding its
        //  group [m3/s]
        scalar forcedFlow(const label k, const scalar t) const;

        //- Its time derivative, from those of the sources [m3/s2]
        scalar forcedFlowRate(const label k, const scalar t) const;


public:

    // Constructors

        //- Construct from circulationProperties, with time starting at t0
        circulationNetwork(const dictionary& dict, const scalar t0);

        //- No copy construct
        circulationNetwork(const circulationNetwork&) = delete;

        //- No copy assignment
        void operator=(const circulationNetwork&) = delete;


    // Member Functions

        // ODESystem

            //- Number of ODE states
            virtual label nEqns() const { return nStates_; }

            //- Time derivative of the volumes and inductor flows at (x, y)
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Finite-difference Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;


        // Time stepping

            //- Current time
            scalar time() const { return t_; }

            //- Current state [SI]
            const scalarField& state() const { return y_; }

            //- Integrate from the current time to t1
            void advance(const scalar t1);

            //- Replace the state (SI), e.g. when restarting
            void setState(const scalar t, const scalarField& y);

            //- The ODE solver's current step-size estimate. Part of what a
            //  restart must carry for the resumed run to take the very same
            //  steps (and so find hard-valve switching instants where the
            //  continuous run did).
            scalar stepEstimate() const { return dxEst_; }

            //- Set the ODE solver's step-size estimate
            void setStepEstimate(const scalar dx) { dxEst_ = dx; }


        // Results

            //- Unit system and cardiac period of the case
            const circulationContext& context() const { return ctx_; }

            //- Column headings, "<component>.P" etc, with units in the case's
            //  system
            const wordList& columns() const { return columns_; }

            //- Unit label of each column
            stringList columnUnits() const;

            //- Current values in the case's units, in the order of columns()
            void outputValues(scalarField& values) const;


        // Checkpoint

            //- State by component name [SI]
            dictionary stateDict() const;

            //- Restore from stateDict(); also sets the time from "time"
            void readState(const dictionary& d);


        // Coupling ports

            //- Sorted names of the coupling ports
            wordList portNames() const { return portIndex_.sortedToc(); }

            //- Supply the flow of an externalPort over [t0, t1], SI
            void setPortForcing
            (
                const word& port,
                const scalar t0,
                const scalar t1,
                const externalPort::forcingCoeffs& c
            );

            //- Pressure to impose at the patch of a port, at the current
            //  state and using the port flow at the current time [Pa]
            scalar portPressure(const word& port) const;

            //- True if the port is an inlet (flow positive into the domain)
            bool portInlet(const word& port) const;

            //- Name of the CFD patch the port is coupled to (may be empty)
            const word& portPatch(const word& port) const;

            //- True if the flow of the port is supplied from outside, by the
            //  coupling code, and not by a table
            bool portExternal(const word& port) const;

            //- Flow through a port at time t [m3/s]
            scalar portFlow(const word& port, const scalar t) const;


        // Info

            //- Print a short summary: element, node and state counts, units,
            //  ODE solver and coupling ports
            void summary(Ostream& os) const;

            //- Print the whole set-up: units, elements, nodes, states,
            //  coupling ports and solver, and warn about parameter values
            //  outside the human-scale ranges (the usual sign of a unit
            //  mix-up). If patch names are given, each port's patch is
            //  checked against them. Returns the number of warnings.
            label check
            (
                Ostream& os,
                const wordList* meshPatches = nullptr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
