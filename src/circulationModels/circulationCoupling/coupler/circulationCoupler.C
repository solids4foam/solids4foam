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
    Samples the flow at each coupled patch, advances the circulation network
    with the explicit or implicit coupling scheme of the chosen order, and
    returns the pressure to impose at each patch. It also writes the results
    and the restart checkpoint.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationCoupler.H"
#include "fvMesh.H"
#include "Time.H"
#include "IOdictionary.H"
#include "surfaceFields.H"
#include "Pstream.H"

#include <map>
#include <cmath>
#include <csignal>
#include <limits>
#include <cstdio>
#include <cstring>
#include <unistd.h>

// * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * * //

namespace
{

// The signal a floating point trap raises. OpenFOAM does the same on Apple
// Silicon, where the trap arrives as SIGILL rather than SIGFPE.
#if defined(__APPLE__) && defined(__arm64__)
constexpr int trapSignal = SIGILL;
#else
constexpr int trapSignal = SIGFPE;
#endif

// What the signal handler needs to describe the failure. Plain data, updated
// after every exchange.
struct fpeInfo
{
    double t = 0, dt = 0, Pmax = 0;
    int order = 0;
    bool implicitScheme = false;
    bool installed = false;
} fpe;

struct sigaction previousFpe;

void fpeAdvice(int sig, siginfo_t* info, void* context)
{
    char buf[1400];
    const int n = std::snprintf
    (
        buf, sizeof(buf),
        "\n--> circulation coupling: floating point exception in the CFD "
        "solution\n"
        "    after t = %g (last accepted step), %s coupling, order %d, "
        "deltaT %g s.\n"
        "    The largest port pressure at the last exchange was %g Pa. "
        "If port pressures\n"
        "    or flows were growing before this, the coupling is probably "
        "unstable here:\n"
        "      - explicit coupling is only conditionally stable: "
        "use  scheme implicit;\n"
        "        or a smaller deltaT\n"
        "      - a high coupling order amplifies noise at a large deltaT: "
        "try order 2\n"
        "      - check the circuit with  circulationFoam -check\n"
        "    (the original CFD error follows)\n\n",
        fpe.t,
        fpe.implicitScheme ? "implicit" : "explicit",
        fpe.order,
        fpe.dt,
        fpe.Pmax
    );
    if (n > 0)
    {
        const ssize_t w = ::write(2, buf, static_cast<size_t>(n));
        (void)w;
    }

    // Hand on to the handler that was installed before ours (OpenFOAM's, which
    // prints the stack trace)
    if (previousFpe.sa_flags & SA_SIGINFO)
    {
        if (previousFpe.sa_sigaction)
        {
            previousFpe.sa_sigaction(sig, info, context);
        }
    }
    else if
    (
        previousFpe.sa_handler
     && previousFpe.sa_handler != SIG_IGN
     && previousFpe.sa_handler != SIG_DFL
    )
    {
        previousFpe.sa_handler(sig);
    }
    ::signal(trapSignal, SIG_DFL);
    ::raise(trapSignal);
}

//- Make sure our advice handler is the active one. OpenFOAM re-arms its own
//  handler at times (sigFpe::ignore guards), so this is checked every exchange.
void installFpeAdvice()
{
    struct sigaction cur;
    if (sigaction(trapSignal, nullptr, &cur) != 0)
    {
        return;
    }
    if ((cur.sa_flags & SA_SIGINFO) && cur.sa_sigaction == fpeAdvice)
    {
        return;
    }

    struct sigaction act;
    std::memset(&act, 0, sizeof(act));
    act.sa_sigaction = fpeAdvice;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;

    struct sigaction before;
    if (sigaction(trapSignal, &act, &before) == 0)
    {
        previousFpe = before;
    }
    fpe.installed = true;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(circulationCoupler, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationCheckpoint::circulationCheckpoint
(
    const fvMesh& mesh,
    const word& name
)
:
    localIOdictionary
    (
        IOobject
        (
            name,
            mesh.time().timeName(mesh.time().startTime().value()),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null
    ),
    resumed_(found("y"))
{}


bool Foam::circulationCheckpoint::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    streamOpt.format(IOstreamOption::ASCII);

    // Full double precision, whatever writePrecision the case uses: this is the
    // state a restart resumes from, so rounding it to the field precision
    // would offset the resumed run by that much
    const unsigned int oldPrecision = IOstream::defaultPrecision
    (
        std::numeric_limits<scalar>::max_digits10
    );
    const bool ok = localIOdictionary::writeObject(streamOpt, writeOnProc);
    IOstream::defaultPrecision(oldPrecision);
    return ok;
}


Foam::circulationCoupler::circulationCoupler
(
    const word& propertiesName,
    const fvMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            propertiesName + ".coupler",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    scheme_(scheme::implicitScheme),
    order_(1),
    rho_(1),
    maxPressure_(1e8),
    historyDt_(0.0),
    historyDepth_(0),
    lastT_(mesh.time().startTime().value()),
    stepStartT_(lastT_),
    activeTimeIndex_(-1),
    lastConsumedFluxEvent_(-1),
    stepUpdateCount_(0),
    totalUpdateCount_(0)
{
    const IOdictionary props
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& cd = props.subDict("coupling");

    const word sch = cd.get<word>("scheme");
    if (sch == "explicit")
    {
        scheme_ = scheme::explicitScheme;
    }
    else if (sch == "implicit")
    {
        scheme_ = scheme::implicitScheme;
    }
    else
    {
        FatalIOErrorInFunction(cd)
            << "scheme must be explicit or implicit, not '" << sch << "'"
            << exit(FatalIOError);
    }

    order_ = cd.get<label>("order");
    if (order_ < 1 || order_ > 4)
    {
        FatalIOErrorInFunction(cd)
            << "order must be 1, 2, 3 or 4, not " << order_
            << exit(FatalIOError);
    }

    rho_ = cd.get<scalar>("rho");
    maxPressure_ = cd.getOrDefault<scalar>("maxPressure", 1e8);

    // The CFD time scheme should be at least as accurate as the coupling
    {
        const std::map<std::string, label> schemeOrder =
        {
            {"Euler", 1}, {"localEuler", 1}, {"backward", 2},
            {"CrankNicolson", 2}
        };

        for (const entry& e : mesh.schemes().ddtSchemes())
        {
            const ITstream& is = e.stream();
            if (is.empty() || !is[0].isWord())
            {
                continue;
            }

            std::string name = is[0].wordToken();
            label i = 0;
            if (name == "bounded" && is.size() > 1 && is[1].isWord())
            {
                i = 1;
                name = is[1].wordToken();
            }

            // The BDF scheme takes its order as an argument, e.g. "BDF 3"
            label cfdOrder = -1;
            if (name == "BDF" && is.size() > i + 1 && is[i + 1].isLabel())
            {
                cfdOrder = is[i + 1].labelToken();
                name += " " + std::to_string(cfdOrder);
            }
            else
            {
                const auto it = schemeOrder.find(name);
                if (it != schemeOrder.end())
                {
                    cfdOrder = it->second;
                }
            }

            if (cfdOrder > 0 && cfdOrder < order_ && Pstream::master())
            {
                WarningInFunction
                    << "The coupling is order " << order_ << " but the CFD time"
                    << " scheme for '" << e.keyword() << "' is " << name
                    << ", order " << cfdOrder << ". The coupled solution is"
                    << " limited to order " << cfdOrder << "; use a CFD time"
                    << " scheme of order " << order_ << " (e.g. BDF " << order_
                    << ", if available) or lower the coupling order." << endl;
            }
        }
    }
    network_.reset(new circulationNetwork(props, lastT_));

    // Bind every port to its patch
    for (const word& pn : network_().portNames())
    {
        const word& patchName = network_().portPatch(pn);

        if (patchName.empty())
        {
            FatalIOErrorInFunction(props)
                << "Port '" << pn << "' has no \"patch\" entry, so it is not"
                << " coupled to a CFD patch" << exit(FatalIOError);
        }

        const label patchID = mesh.boundaryMesh().findPatchID(patchName);
        if (patchID < 0)
        {
            FatalIOErrorInFunction(props)
                << "Port '" << pn << "' names patch '" << patchName
                << "', which is not in the mesh. Patches: "
                << mesh.boundaryMesh().names() << exit(FatalIOError);
        }

        portData p;
        p.name = pn;
        p.patchID = patchID;
        p.inlet = network_().portInlet(pn);
        p.Q = 0;
        p.stepStartQ = 0;
        p.history = 0.0;
        p.P = 0;
        ports_.append(p);
    }

    if (ports_.empty())
    {
        FatalIOErrorInFunction(props)
            << "The network has no coupling ports (externalPort or"
            << " coupledDomain)" << exit(FatalIOError);
    }

    // Optional history seeds, most recent first. "seed zero;" says the system
    // starts at rest, so every earlier flow is zero and the coupling is at full
    // order from the first step.
    if (cd.found("seed") && !cd.isDict("seed"))
    {
        const word how = cd.get<word>("seed");

        if (how != "zero")
        {
            FatalIOErrorInFunction(cd)
                << "seed must be \"zero\" or a dictionary of earlier flows per"
                << " port, not '" << how << "'" << exit(FatalIOError);
        }

        for (portData& p : ports_)
        {
            p.history = 0.0;
        }
        historyDepth_ = 3;
        historyDt_ = mesh.time().deltaTValue();
    }
    else if (cd.isDict("seed"))
    {
        const dictionary& sd = cd.subDict("seed");
        const scalar sdt = cd.get<scalar>("seedDeltaT");

        label depth = -1;
        for (portData& p : ports_)
        {
            scalarList q(sd.get<scalarList>(p.name));
            if (q.size() > 3)
            {
                FatalIOErrorInFunction(sd)
                    << "At most three seed flows, port " << p.name
                    << exit(FatalIOError);
            }
            if (depth >= 0 && q.size() != depth)
            {
                FatalIOErrorInFunction(sd)
                    << "Every port needs the same number of seed flows"
                    << exit(FatalIOError);
            }
            depth = q.size();
            forAll(q, i)
            {
                p.history[i] = q[i];
            }
        }
        historyDepth_ = depth;
        for (label i = 0; i < historyDepth_; ++i)
        {
            historyDt_[i] = sdt;
        }
    }

    y_ = network_().state();
    stepStartY_ = y_;

    // On a floating point exception, say what it probably means
    installFpeAdvice();
    fpe.order = order_;
    fpe.implicitScheme = (scheme_ == scheme::implicitScheme);
    fpe.dt = mesh.time().deltaTValue();

    // Pressure before the first exchange: the initial state with no flow
    for (portData& p : ports_)
    {
        p.P = network_().portPressure(p.name);
    }

    // Resume from the checkpoint of the start time, if there is one
    checkpoint_.reset
    (
        new circulationCheckpoint(mesh, "circulationState-" + propertiesName)
    );
    if (checkpoint_().resumed())
    {
        restore();
    }

    if (Pstream::master())
    {
        const fileName dir
        (
            mesh.time().globalPath()/"postProcessing"/"circulationCoupling"
           /mesh.time().timeName(mesh.time().startTime().value())
        );
        mkDir(dir);
        log_.reset(new OFstream(dir/"coupling.dat"));
        log_().precision(12);

        log_() << "# circulationCoupler: " << sch << " order " << order_ << nl
               << "# time updates";
        for (const word& c : network_().columns())
        {
            log_() << ' ' << c;
        }
        for (const portData& p : ports_)
        {
            log_() << ' ' << p.name << ".Qsampled";
        }
        log_() << nl;

        if (props.isDict("metrics"))
        {
            metrics_.reset
            (
                new circulationMetrics
                (
                    props.subDict("metrics"),
                    network_().columns(),
                    network_().context(),
                    dir
                )
            );
        }

        Info<< type() << ": " << sch << " coupling, order " << order_
            << ", ports:";
        for (const portData& p : ports_)
        {
            Info<< ' ' << p.name << "->" << network_().portPatch(p.name);
        }
        Info<< nl;
        network_().summary(Info);
    }
}


Foam::circulationCoupler::~circulationCoupler()
{
    // The last accepted step is only logged when the next one begins
    if (activeTimeIndex_ >= 0 && stepUpdateCount_ > 0)
    {
        writeRow();
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::circulationCoupler& Foam::circulationCoupler::lookupOrCreate
(
    const fvMesh& mesh,
    const word& propertiesName
)
{
    objectRegistry& obr = const_cast<objectRegistry&>(mesh.thisDb());
    const word key = propertiesName + ".coupler";

    if (obr.foundObject<circulationCoupler>(key))
    {
        return obr.lookupObjectRef<circulationCoupler>(key);
    }

    circulationCoupler* ptr = new circulationCoupler(propertiesName, mesh);
    obr.store(ptr);
    return *ptr;
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::circulationCoupler::writeRow()
{
    if (!log_.valid())
    {
        return;
    }

    scalarField values;
    network_().outputValues(values);

    log_() << lastT_ << ' ' << stepUpdateCount_;
    for (const scalar v : values)
    {
        log_() << ' ' << v;
    }

    // The flows the CFD delivered, as opposed to the forcing the scheme built
    // from them, in the units of the network description
    const scalar fQ =
        network_().context().factor(circulationContext::quantity::flow);
    for (const portData& p : ports_)
    {
        log_() << ' ' << p.Q/fQ;
    }

    log_() << nl;
    log_().flush();

    if (metrics_.valid())
    {
        metrics_().add(lastT_, values);
    }
}


void Foam::circulationCoupler::updateCheckpoint()
{
    if (!checkpoint_.valid())
    {
        return;
    }

    // What the next step expects to find: the history as it will be once the
    // step just accepted has been rolled into it (see beginTimeIndex)
    const bool imp = scheme_ == scheme::implicitScheme;

    const unsigned int oldPrecision = IOstream::defaultPrecision
    (
        std::numeric_limits<scalar>::max_digits10
    );

    circulationCheckpoint& cp = checkpoint_();
    cp.set("time", lastT_);
    cp.set("y", y_);
    cp.set("stepEstimate", network_().stepEstimate());
    cp.set("historyDepth", min(historyDepth_ + 1, label(3)));

    FixedList<scalar, 3> dt;
    dt[0] = lastT_ - stepStartT_;
    dt[1] = historyDt_[0];
    dt[2] = historyDt_[1];
    cp.set("historyDt", List<scalar>(dt));

    dictionary ports;
    for (const portData& p : ports_)
    {
        dictionary pd;
        pd.add("Q", p.Q);
        FixedList<scalar, 3> h;
        h[0] = imp ? p.stepStartQ : p.Q;
        h[1] = p.history[0];
        h[2] = p.history[1];
        pd.add("history", List<scalar>(h));
        ports.add(p.name, pd);
    }
    cp.set("ports", ports);

    IOstream::defaultPrecision(oldPrecision);
}


void Foam::circulationCoupler::restore()
{
    const circulationCheckpoint& cp = checkpoint_();

    const scalarField y(cp.get<scalarField>("y"));
    if (y.size() != network_().nEqns())
    {
        FatalErrorInFunction
            << "Checkpoint " << cp.objectRelPath() << " holds " << y.size()
            << " states but the circuit has " << network_().nEqns()
            << ": the circuit was changed since it was written"
            << exit(FatalError);
    }

    lastT_ = cp.get<scalar>("time");
    y_ = y;
    network_().setState(lastT_, y_);
    network_().setStepEstimate
    (
        cp.getOrDefault<scalar>("stepEstimate", network_().stepEstimate())
    );
    stepStartY_ = y_;
    stepStartT_ = lastT_;

    historyDepth_ = cp.get<label>("historyDepth");
    const List<scalar> hdt(cp.get<List<scalar>>("historyDt"));
    for (label i = 0; i < 3; ++i)
    {
        historyDt_[i] = hdt[i];
    }

    const dictionary& pd = cp.subDict("ports");
    for (portData& p : ports_)
    {
        const dictionary& d = pd.subDict(p.name);
        p.Q = d.get<scalar>("Q");
        p.stepStartQ = p.Q;
        const List<scalar> h(d.get<List<scalar>>("history"));
        for (label i = 0; i < 3; ++i)
        {
            p.history[i] = h[i];
        }
        p.P = network_().portPressure(p.name);
    }

    Info<< type() << ": resumed from " << cp.objectRelPath() << " at t = "
        << lastT_ << endl;
}


void Foam::circulationCoupler::beginTimeIndex(const label timeIndex)
{
    if (timeIndex == activeTimeIndex_)
    {
        return;
    }

    if (activeTimeIndex_ >= 0)
    {
        // The previous step is complete: its start flow (implicit) or the flow
        // it sampled (explicit) becomes the newest history entry
        const bool imp = scheme_ == scheme::implicitScheme;

        for (portData& p : ports_)
        {
            p.history[2] = p.history[1];
            p.history[1] = p.history[0];
            p.history[0] = imp ? p.stepStartQ : p.Q;
        }
        historyDt_[2] = historyDt_[1];
        historyDt_[1] = historyDt_[0];
        historyDt_[0] = lastT_ - stepStartT_;
        historyDepth_ = min(historyDepth_ + 1, label(3));

        // One row per accepted step
        network_().setState(lastT_, y_);
        writeRow();
    }

    activeTimeIndex_ = timeIndex;

    // Freeze the accepted state and flows: every attempt in this step
    // integrates from here, so a repeated attempt cannot march the network
    // forward twice
    stepStartY_ = y_;
    stepStartT_ = lastT_;
    for (portData& p : ports_)
    {
        p.stepStartQ = p.Q;
    }

    stepUpdateCount_ = 0;
    lastConsumedFluxEvent_ = -1;
}


void Foam::circulationCoupler::advance(const Time& runTime)
{
    const scalar t0 = stepStartT_;
    const scalar t1 = runTime.value();
    const scalar dt = t1 - t0;

    if (dt <= SMALL)
    {
        return;
    }

    const bool imp = scheme_ == scheme::implicitScheme;

    // ---- Forcing polynomial of each port ------------------------------------
    for (portData& p : ports_)
    {
        scalarList alpha, Q;

        if (!imp)
        {
            // Q(t0) and up to order-1 earlier flows, extrapolated across the
            // step
            const label n = min(order_, 1 + historyDepth_);
            alpha.setSize(n);
            Q.setSize(n);
            alpha[0] = 0;
            Q[0] = p.Q;
            for (label j = 1; j < n; ++j)
            {
                const scalar h =
                    historyDt_[j-1] > VSMALL ? historyDt_[j-1] : dt;
                alpha[j] = alpha[j-1] - h/dt;
                Q[j] = p.history[j-1];
            }
        }
        else if (order_ == 1)
        {
            // The current iterate, held constant
            alpha.setSize(1, scalar(1));
            Q.setSize(1, p.Q);
        }
        else
        {
            // The iterate at alpha 1, the accepted start flow at 0, and up to
            // order-2 earlier flows
            const label n = 2 + min(order_ - 2, historyDepth_);
            alpha.setSize(n);
            Q.setSize(n);
            alpha[0] = 1;
            Q[0] = p.Q;
            alpha[1] = 0;
            Q[1] = p.stepStartQ;
            for (label j = 2; j < n; ++j)
            {
                const scalar h =
                    historyDt_[j-2] > VSMALL ? historyDt_[j-2] : dt;
                alpha[j] = alpha[j-1] - h/dt;
                Q[j] = p.history[j-2];
            }
        }

        if (stepUpdateCount_ == 0 && Pstream::master() && alpha.size() < order_)
        {
            Info<< type() << ": only " << historyDepth_ << " earlier flows at"
                << " time index " << runTime.timeIndex() << ", port " << p.name
                << "; using order-" << alpha.size() << " forcing."
                << (historyDepth_ == 0
                  ? "  If the system starts at rest, \"seed zero;\""
                    " in the coupling dictionary starts at full order."
                  : "")
                << endl;
        }

        network_().setPortForcing
        (
            p.name, t0, t1, circulationForcing::polynomialThrough(alpha, Q)
        );
    }

    // ---- Integrate from the start of the step ---------------------------
    network_().setState(t0, stepStartY_);
    network_().advance(t1);
    const scalarField candidate = network_().state();

    // Both schemes accept the integrated state whole
    y_ = candidate;

    lastT_ = t1;
    network_().setState(t1, y_);

    for (portData& p : ports_)
    {
        p.P = network_().portPressure(p.name);
    }

    // Divergence guard: stop at a clear threshold with advice, rather than
    // when the CFD solver finally overflows
    for (const portData& p : ports_)
    {
        if (!std::isfinite(p.P) || Foam::mag(p.P) > maxPressure_)
        {
            FatalErrorInFunction
                << "The coupling has diverged: the pressure at port '"
                << p.name
                << "' is " << p.P << " Pa at t = " << t1 << " (limit "
                << maxPressure_ << " Pa)." << nl
                << "    Scheme " << (imp ? "implicit" : "explicit")
                << ", order "
                << order_ << ", deltaT " << dt << " s. Likely causes:" << nl
                << "      - "
                << (imp ? "the outer iteration is not converging"
                        : "explicit coupling is only conditionally stable")
                << ": " << (imp ? "raise nOuterCorrectors or reduce deltaT"
                                : "use  scheme implicit;  or reduce deltaT")
                << nl
                << "      - a high coupling order amplifies flow noise at a "
                   "large deltaT: try order 2" << nl
                << "      - the circuit parameters are unphysical "
                   "(see circulationFoam -check)" << nl
                << "    If this pressure is real, raise  maxPressure  in the "
                   "coupling dictionary."
                << exit(FatalError);
        }
    }

    installFpeAdvice();
    updateCheckpoint();

    fpe.t = t1;
    fpe.Pmax = 0;
    for (const portData& p : ports_)
    {
        fpe.Pmax = std::max(fpe.Pmax, std::abs(double(p.P)));
    }

    ++stepUpdateCount_;
    ++totalUpdateCount_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::circulationCoupler::portOfPatch
(
    const label patchID
) const
{
    for (const portData& p : ports_)
    {
        if (p.patchID == patchID)
        {
            return p.name;
        }
    }

    FatalErrorInFunction
        << "Patch " << mesh_.boundaryMesh()[patchID].name() << " is not"
        << " coupled to any port of the circulation network" << nl
        << "    (a port needs a \"patch "
        << mesh_.boundaryMesh()[patchID].name()
        << ";\" entry)" << exit(FatalError);

    return word::null;
}


void Foam::circulationCoupler::sampleAndAdvance(const fvMesh& mesh)
{
    const Time& runTime = mesh.time();

    beginTimeIndex(runTime.timeIndex());

    if (scheme_ == scheme::explicitScheme)
    {
        // One exchange per time step; later evaluations (outer, inner and
        // non-orthogonal correctors) reuse the pressure already computed
        if (stepUpdateCount_ > 0)
        {
            return;
        }
    }
    else
    {
        // One exchange per new flux field. phi's event number identifies it;
        // the repeated boundary evaluations of one outer corrector all see the
        // same number. Every rank must reach the same decision, so the
        // master's is broadcast.
        const label event =
            mesh.lookupObject<surfaceScalarField>("phi").eventNo();
        bool fresh = (event != lastConsumedFluxEvent_);
        Pstream::broadcast(fresh);
        if (!fresh)
        {
            return;
        }
    }

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() != dimVolume/dimTime)
    {
        FatalErrorInFunction
            << "A volumetric flux phi is required, found " << phi.dimensions()
            << exit(FatalError);
    }

    // All ranks take part in the reductions, including those without faces
    for (portData& p : ports_)
    {
        const scalar leaving = gSum(phi.boundaryField()[p.patchID]);
        p.Q = p.inlet ? -leaving : leaving;
    }

    lastConsumedFluxEvent_ = phi.eventNo();

    advance(runTime);
}


Foam::scalar Foam::circulationCoupler::patchPressure(const word& port) const
{
    for (const portData& p : ports_)
    {
        if (p.name == port)
        {
            return p.P;
        }
    }

    FatalErrorInFunction << "No port '" << port << "'" << exit(FatalError);
    return 0;
}

// ************************************************************************* //
