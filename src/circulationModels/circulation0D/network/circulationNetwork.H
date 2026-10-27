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

namespace Foam
{

class circulationNetwork
:
    public ODESystem
{
    typedef circulationComponent::role role;

    // Private Data

        circulationContext ctx_;

        //- Primitive components, composites already expanded
        PtrList<circulationComponent> comps_;

        List<role> roles_;

        //- Net of the "in" and "out" terminal of each component
        labelList netIn_;
        labelList netOut_;

        label nNets_;

        //- Net -> component holding its volume state / fixing its pressure
        labelList netStorage_;
        labelList netFixed_;

        //- Net -> index in the algebraic unknowns, -1 if not algebraic
        labelList algIndex_;
        label nAlg_;

        //- True if every branch touching an algebraic net is linear
        bool algLinear_;

        //- An inductor whose flow is not a state because it carries the whole
        //  flow of a group of algebraic nodes that has no pressure reference
        //  (for example a flow source feeding r -> L -> C: q_L = Q(t)). Its
        //  pressure drop L dq/dt pins the group instead.
        struct forcedInertia
        {
            label comp;         //!< the inductor
            label pinnedNet;    //!< its end inside the group
            label otherNet;     //!< its end outside
            scalar sign;        //!< +1 if the pinned end is "in", else -1
            label group;
        };

        List<forcedInertia> forced_;
        labelList forcedOfComp_;    //!< component -> index in forced_, or -1
        labelList groupOfNet_;      //!< algebraic net -> resistive group
        boolList pinned_;           //!< net whose nodal row is replaced

        //- Component -> its state index (storage: volume, inertia: flow)
        labelList stateOfComp_;
        label nStates_;

        //- Port name -> component index
        HashTable<label> portIndex_;

        dictionary odeDict_;
        autoPtr<ODESolver> solver_;

        scalar t_;
        scalarField y_;
        scalar dxEst_;

        //- Work space; derivatives() is const but warm-starts the nodal solve
        mutable scalarField pAlg_;
        mutable scalarField Pwork_;

        wordList columns_;

        //- Nets whose pressure is reported (the optional "probes" list)
        labelList probeNet_;


    // Private Member Functions

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

        scalar branchFlow
        (
            const label k,
            const scalar t,
            const scalarField& y,
            const scalarField& P
        ) const;

        //- Names of the components attached to a net, for messages
        string netMembers(const label net) const;

        label portLabel(const word& port) const;

        //- Flow in -> out of a forced inductor, from the sources feeding its
        //  group [m3/s]
        scalar forcedFlow(const label k, const scalar t) const;

        //- Its time derivative, by central difference [m3/s2]
        scalar forcedFlowRate(const label k, const scalar t) const;


public:

    // Constructors

        //- Construct from circulationProperties, with time starting at t0
        circulationNetwork(const dictionary& dict, const scalar t0);

        circulationNetwork(const circulationNetwork&) = delete;
        void operator=(const circulationNetwork&) = delete;


    // ODESystem

        virtual label nEqns() const { return nStates_; }

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

        scalar time() const { return t_; }

        const scalarField& state() const { return y_; }

        //- Integrate from the current time to t1
        void advance(const scalar t1);

        //- Replace the state (SI), e.g. when restarting
        void setState(const scalar t, const scalarField& y);

        //- The ODE solver's current step-size estimate. Part of what a restart
        //  must carry for the resumed run to take the very same steps (and so
        //  find hard-valve switching instants where the continuous run did).
        scalar stepEstimate() const { return dxEst_; }
        void setStepEstimate(const scalar dx) { dxEst_ = dx; }


    // Results

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

        wordList portNames() const { return portIndex_.sortedToc(); }

        //- Supply the flow of an externalPort over [t0, t1], SI
        void setPortForcing
        (
            const word& port,
            const scalar t0,
            const scalar t1,
            const externalPort::forcingCoeffs& c
        );

        //- Pressure to impose at the patch of a port, at the current state
        //  and using the port flow at the current time [Pa]
        scalar portPressure(const word& port) const;

        //- True if the port is an inlet (flow positive into the domain)
        bool portInlet(const word& port) const;

        //- Name of the CFD patch the port is coupled to (may be empty)
        const word& portPatch(const word& port) const;

        //- Flow through a port at time t [m3/s]
        scalar portFlow(const word& port, const scalar t) const;


    // Info

        void summary(Ostream& os) const;

        //- Print the whole set-up: units, elements, nodes, states, coupling
        //  ports and solver, and warn about parameter values outside the
        //  human-scale ranges (the usual sign of a unit mix-up). If patch
        //  names are given, each port's patch is checked against them.
        //  Returns the number of warnings.
        label check(Ostream& os, const wordList* meshPatches = nullptr) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
