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
    Foam::circulationCoupler

Description
    Couples the patches of a CFD case to a circulationNetwork, with the
    coupling algorithm chosen in the case dictionary instead of by the type of
    boundary condition. One coupler serves every coupled patch of a case and
    is shared by their circulationPressure boundary conditions.

    In constant/circulationProperties:

        coupling
        {
            scheme      implicit;   // explicit | implicit
            order       2;          // 1..4: degree of the forcing + 1
            rho         1060;       // kg/m3: converts Pa to a kinematic p
            maxPressure 1e8;        // optional [Pa]: stop if a port pressure
                                    // exceeds this
        }

    The network sees each coupled port's flow over a step [t0, t1] as a
    polynomial in alpha = (t - t0)/(t1 - t0), through samples of the flow:

    \table
        scheme    | order | samples (alpha)
        explicit  | 1     | Q(t0) held constant
        explicit  | 2     | Q(t0), Q(t0-dt)             extrapolated to t1
        explicit  | 3, 4  | Q(t0) and two / three earlier flows, extrapolated
        implicit  | 1     | the current outer iterate, held constant
        implicit  | 2     | Q(t0) at alpha 0, the iterate at alpha 1
        implicit  | 3, 4  | the above and one / two earlier flows
    \endtable

    Explicit coupling integrates the network once per time step. Implicit
    coupling re-integrates from the start of the step each time the CFD flux is
    updated, which makes the outer correctors of PIMPLE a fixed-point iteration
    on the interface flow. Orders 3 and 4 need earlier flows; a fresh run has
    fewer, and the order is reduced for those first steps (a message says so).
    They may be seeded:

        seed        { outlet ( Q(-dt) Q(-2dt) Q(-3dt) ); }  // most recent first
        seedDeltaT  0.01;

    or, when the system starts at rest so that every earlier flow is zero,

        seed        zero;

    At start-up the coupler also checks the CFD time scheme: a scheme of lower
    order than the coupling (for example "backward" with order 3) limits the
    coupled solution, and a warning says so.

    Ports are bound to patches by the "patch" entry of the port in the network
    description (see externalPort). Sign convention: Q is positive in the
    natural direction of the port (outlet: leaving the CFD domain, so the sum
    of phi over the patch; inlet: entering it, minus that sum). In parallel the
    flows are summed over all ranks, the network is integrated on the master,
    and its state and port pressures are broadcast.

    Results: postProcessing/circulationCoupling/<start>/coupling.dat, one row
    per accepted time step, written as soon as the step is known to be
    complete: with explicit coupling straight after its exchange; with implicit
    coupling, whose state can change at every outer corrector, when the next
    step starts, at every write time and at the end of the run, so that the
    file always holds every step up to the last time written. With implicit
    coupling, the relative change of the port flows |Q^k - Q^(k-1)|/|Q^k| is
    logged at every exchange, and the last one of each step is the "residual"
    column.

    Restart: the state of the network, by component name, and the flow history
    are written, with the fields, to
    <time>/uniform/circulationState-<properties name> at every output time and
    read back when the case starts from a time that has one, so the circuit
    and the CFD resume together and the run continues exactly as it would
    have. (A "seed" is then ignored.) The master writes it, and decomposePar
    and reconstructPar carry it with the uniform directory. A case that starts
    after an earlier time without one is an error.

    circulationFoam writes the same file without the flow history, so a 0D
    spin-up can start a coupled run: copy it to <start time>/uniform/ of the
    CFD case. Components are matched by name, and the state of a component
    that the coupled circuit does not have (one that the CFD domain replaces)
    is ignored. The state is taken at the start time, so spin up whole cycles
    of a periodic circuit. The flow history starts as in a fresh run.

Author
    David Kelly, UCD.

SourceFiles
    circulationCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef circulationCoupler_H
#define circulationCoupler_H

#include "regIOobject.H"
#include "circulationNetwork.H"
#include "circulationForcing.H"
#include "autoPtr.H"
#include "OFstream.H"
#include "circulationMetrics.H"
#include "IOdictionary.H"
#include "IOstreamOption.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class Time;

/*---------------------------------------------------------------------------*\
                    Class circulationCheckpoint Declaration
\*---------------------------------------------------------------------------*/

//- The restart checkpoint: a dictionary kept up to date at every exchange and
//  written automatically by the master, at full precision, with every output
//  time, to <time>/uniform/ so that decomposePar and reconstructPar carry it
class circulationCheckpoint
:
    public IOdictionary
{
public:

    // Constructors

        //- Construct for the mesh, reading the checkpoint "name" from
        //  <start time>/uniform/ if present; in parallel the master's copy
        //  is broadcast to every rank
        circulationCheckpoint(const fvMesh& mesh, const word& name);


    // Member Functions

        //- True if a checkpoint was read at construction
        bool resumed() const { return resumed_; }

        //- Write the checkpoint, on the master only, in ASCII at full
        //  precision whatever the case's write format and precision
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;


private:

    // Private Data

        //- True if a checkpoint (with a "state" entry) was read at
        //  construction
        bool resumed_;
};


/*---------------------------------------------------------------------------*\
                     Class circulationCoupler Declaration
\*---------------------------------------------------------------------------*/

class circulationCoupler
:
    public regIOobject
{
public:

    //- Coupling scheme: explicit (one exchange per time step) or implicit
    //  (one exchange per new flux field, re-integrating the step)
    enum class scheme { explicitScheme, implicitScheme };

private:

    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the volumetric flux sampled at the coupled patches
        word phiName_;

        //- Coupling scheme
        scheme scheme_;

        //- Coupling order (1 to 4): degree of the forcing polynomial + 1
        label order_;

        //- Density used to convert Pa to a kinematic pressure [kg/m3]
        scalar rho_;

        //- Pressure at a port beyond which the coupling is taken to have
        //  diverged [Pa]
        scalar maxPressure_;

        //- The circulation network
        autoPtr<circulationNetwork> network_;

        //- One coupled patch
        struct portData
        {
            //- Name of the port
            word name;

            //- Index of the patch coupled to the port
            label patchID;

            //- True for an inlet: its flow is minus the sum of phi
            bool inlet;

            //- Flow most recently sampled
            scalar Q;

            //- Accepted flow at the step start
            scalar stepStartQ;

            //- Earlier flows, most recent first
            FixedList<scalar, 3> history;

            //- Pressure to impose [Pa]
            scalar P;
        };

        //- The coupled ports
        List<portData> ports_;

        //- Interval belonging to each history entry
        FixedList<scalar, 3> historyDt_;

        //- Number of valid history entries
        label historyDepth_;

        //- Accepted network state
        scalarField y_;

        //- Time of the accepted network state
        scalar lastT_;

        //- Accepted network state at the start of this step
        scalarField stepStartY_;

        //- Time at the start of this step
        scalar stepStartT_;

        //- Time index of the step in progress (-1 before the first)
        label activeTimeIndex_;

        //- Event number of the flux field consumed by the last exchange
        label lastConsumedFluxEvent_;

        //- Number of exchanges in this step
        label stepUpdateCount_;

        //- Number of exchanges in the run
        label totalUpdateCount_;

        //- The state accepted last has not been logged yet
        bool rowPending_;

        //- Implicit coupling: |Q^k - Q^(k-1)|/|Q^k| over the ports at the
        //  last exchange, Q^0 being the flow at the start of the step
        scalar interfaceResidual_;

        //- Restart checkpoint, <time>/uniform/circulationState-<name>
        autoPtr<circulationCheckpoint> checkpoint_;

        //- Results file, coupling.dat (master only)
        autoPtr<OFstream> log_;

        //- Cycle metrics, if a "metrics" dictionary is given (master only)
        autoPtr<circulationMetrics> metrics_;


    // Private Member Functions

        //- Flow of a port from the flux, positive in its natural direction
        scalar portFlow
        (
            const surfaceScalarField& phi,
            const portData& p
        ) const;

        //- Start a new time step if timeIndex is new: roll the previous
        //  step into the flow history, log it, and freeze the accepted
        //  state and flows as the start of this step
        void beginTimeIndex(const label timeIndex);

        //- Build the forcing polynomial of each port and integrate the
        //  network from the start of the step to the current time, on the
        //  master; broadcast the state and port pressures, and stop if a
        //  pressure exceeds maxPressure
        void advance(const Time& runTime);

        //- Write a row of the results file (and add it to the metrics)
        void writeRow();

        //- Log the state accepted last, if it has not been
        void flushRow();

        //- Store what a restart from the end of the step just accepted needs
        void updateCheckpoint();

        //- Resume from a checkpoint found at the start time
        void restore();

        //- No copy construct
        circulationCoupler(const circulationCoupler&) = delete;

        //- No copy assignment
        void operator=(const circulationCoupler&) = delete;


public:

    //- Runtime type information
    TypeName("circulationCoupler");


    // Constructors

        //- Construct reading the circulation description "propertiesName"
        //  from constant/, sampling the flux "phiName"
        circulationCoupler
        (
            const word& propertiesName,
            const fvMesh& mesh,
            const word& phiName = "phi"
        );


    // Destructor

        //- Destructor: logs the state accepted last, if it has not been
        virtual ~circulationCoupler();


    // Selectors

        //- The coupler registered with the mesh, constructed on first use.
        //  Every patch of one coupler must name the same flux.
        static circulationCoupler& lookupOrCreate
        (
            const fvMesh& mesh,
            const word& propertiesName = "circulationProperties",
            const word& phiName = "phi"
        );


    // Member Functions

        //- Name of the port coupled to a patch; fatal if there is none
        const word& portOfPatch(const label patchID) const;

        //- Sample the flows and, on a new exchange, advance the network
        void sampleAndAdvance(const fvMesh& mesh);

        //- Pressure to impose at a port [Pa]
        scalar patchPressure(const word& port) const;

        //- Density used to convert Pa to a kinematic pressure [kg/m3]
        scalar rho() const { return rho_; }

        //- Const access to the circulation network
        const circulationNetwork& network() const { return network_(); }

        //- Called at every write time, when the step is complete: logs its
        //  row. Writes no file.
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        //- Write data: nothing to write
        virtual bool writeData(Ostream&) const { return true; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
