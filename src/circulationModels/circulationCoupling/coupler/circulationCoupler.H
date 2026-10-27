/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

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
            rho         1060;       // kg/m3: the CFD solves p/rho
            maxPressure 1e8;        // optional [Pa]: stop if a port pressure exceeds this
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

        seed        { outlet ( Q(-dt) Q(-2dt) Q(-3dt) ); }    // most recent first
        seedDeltaT  0.01;

    or, when the system starts at rest so that every earlier flow is zero,

        seed        zero;

    At start-up the coupler also checks the CFD time scheme: a scheme of lower
    order than the coupling (for example "backward" with order 3) limits the
    coupled solution, and a warning says so.

    Ports are bound to patches by the "patch" entry of the port in the network
    description (see externalPort). Sign convention: Q is positive in the
    natural direction of the port (outlet: leaving the CFD domain, so the sum
    of phi over the patch; inlet: entering it, minus that sum).

    Results: postProcessing/circulationCoupling/<start>/coupling.dat, one row
    per accepted time step.

    Restart: the state of the network and the flow history are written, with
    the fields, to <time>/circulationState-<properties name> at every output
    time and read back when the case starts from a time that has one, so the
    circuit and the CFD resume together and the run continues exactly as it
    would have. (A "seed" is then ignored.)

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
#include "localIOdictionary.H"
#include "IOstreamOption.H"

namespace Foam
{

class fvMesh;
class Time;

//- The restart checkpoint: a dictionary kept up to date at every exchange and
//  written automatically, at full precision, with every output time
class circulationCheckpoint
:
    public localIOdictionary
{
public:

    circulationCheckpoint(const fvMesh& mesh, const word& name);

    //- True if a checkpoint was read at construction
    bool resumed() const { return resumed_; }

    virtual bool writeObject(IOstreamOption streamOpt, const bool writeOnProc) const;

private:

    bool resumed_;
};


class circulationCoupler
:
    public regIOobject
{
public:

    enum class scheme { explicitScheme, implicitScheme };

private:

    // Private Data

        const fvMesh& mesh_;

        scheme scheme_;
        label order_;
        scalar rho_;

        //- Pressure at a port beyond which the coupling is taken to have
        //  diverged [Pa]
        scalar maxPressure_;

        autoPtr<circulationNetwork> network_;

        //- One coupled patch
        struct portData
        {
            word name;
            label patchID;
            bool inlet;
            scalar Q;                       //!< flow most recently sampled
            scalar stepStartQ;              //!< accepted flow at the step start
            FixedList<scalar, 3> history;   //!< earlier flows, most recent first
            scalar P;                       //!< pressure to impose [Pa]
        };

        List<portData> ports_;

        //- Interval belonging to each history entry, and how many are valid
        FixedList<scalar, 3> historyDt_;
        label historyDepth_;

        //- Accepted network state and its time; and at the start of this step
        scalarField y_;
        scalar lastT_;
        scalarField stepStartY_;
        scalar stepStartT_;

        label activeTimeIndex_;
        label lastConsumedFluxEvent_;
        label stepUpdateCount_;
        label totalUpdateCount_;

        //- Restart checkpoint, <time>/circulationState-<name>
        autoPtr<circulationCheckpoint> checkpoint_;

        autoPtr<OFstream> log_;
        autoPtr<circulationMetrics> metrics_;


    // Private Member Functions

        void beginTimeIndex(const label timeIndex);

        void advance(const Time& runTime);

        void writeRow();

        //- Store what a restart from the end of the step just accepted needs
        void updateCheckpoint();

        //- Resume from a checkpoint found at the start time
        void restore();

        circulationCoupler(const circulationCoupler&) = delete;
        void operator=(const circulationCoupler&) = delete;


public:

    TypeName("circulationCoupler");

    //- Construct reading the circulation description "propertiesName"
    //  from constant/
    circulationCoupler
    (
        const word& propertiesName,
        const fvMesh& mesh
    );

    virtual ~circulationCoupler();

    //- The coupler registered with the mesh, constructed on first use
    static circulationCoupler& lookupOrCreate
    (
        const fvMesh& mesh,
        const word& propertiesName = "circulationProperties"
    );

    //- Name of the port coupled to a patch; fatal if there is none
    const word& portOfPatch(const label patchID) const;

    //- Sample the flows and, on a new exchange, advance the network
    void sampleAndAdvance(const fvMesh& mesh);

    //- Pressure to impose at a port [Pa]
    scalar patchPressure(const word& port) const;

    scalar rho() const { return rho_; }

    const circulationNetwork& network() const { return network_(); }

    virtual bool writeData(Ostream&) const { return true; }
};

} // End namespace Foam

#endif
