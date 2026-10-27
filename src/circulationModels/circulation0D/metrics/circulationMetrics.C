/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "circulationMetrics.H"
#include "OSspecific.H"

#include <cmath>
#include <iomanip>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationMetrics::circulationMetrics
(
    const dictionary& dict,
    const wordList& columns,
    const circulationContext& ctx,
    const fileName& outDir
)
:
    period_(ctx.period()),
    dir_(outDir),
    headerWritten_(false),
    cycle_(-1),
    cycleStartedOnBoundary_(false),
    volumeLabel_(ctx.label(circulationContext::quantity::volume)),
    pressureLabel_(ctx.label(circulationContext::quantity::pressure)),
    flowLabel_(ctx.label(circulationContext::quantity::flow))
{
    auto columnIndex = [&](const dictionary& d, const word& key) -> label
    {
        const word name = d.get<word>(key);
        const label i = columns.find(name);

        if (i < 0)
        {
            FatalIOErrorInFunction(d)
                << "Metrics: there is no output column '" << name << "'. Columns:"
                << nl << columns << exit(FatalIOError);
        }
        return i;
    };

    DynamicList<item> items;
    for (const entry& e : dict)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Metric '" << e.keyword() << "' must be a dictionary"
                << exit(FatalIOError);
        }

        const dictionary& d = e.dict();
        item it;
        it.name = e.keyword();
        it.type = d.get<word>("type");

        if (it.type == "chamber")
        {
            it.cVolume = columnIndex(d, "volume");
            it.cPressure = columnIndex(d, "pressure");
        }
        else if (it.type == "pressure")
        {
            it.cPressure = columnIndex(d, "column");
        }
        else if (it.type == "flow")
        {
            it.cFlow = columnIndex(d, "column");
        }
        else
        {
            FatalIOErrorInFunction(d)
                << "Metric '" << it.name << "': unknown type '" << it.type
                << "'. Valid types: chamber pressure flow" << exit(FatalIOError);
        }
        items.append(it);
    }
    items_ = items;

    // The distinct columns that have to be sampled
    DynamicList<label> mon;
    for (const item& it : items_)
    {
        for (const label c : {it.cVolume, it.cPressure, it.cFlow})
        {
            if (c >= 0 && !mon.found(c)) mon.append(c);
        }
    }
    monitored_ = mon;
    v_.setSize(monitored_.size());

    if (!isDir(dir_)) mkDir(dir_);
    out_.open(dir_/"metrics.dat");
    out_ << std::setprecision(10);
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

Foam::label Foam::circulationMetrics::slot(const label column) const
{
    return monitored_.find(column);
}


void Foam::circulationMetrics::startCycle(const label k)
{
    cycle_ = k;
    t_.clear();
    for (auto& l : v_) l.clear();
}


void Foam::circulationMetrics::finishCycle()
{
    if (t_.size() < 2 || !cycleStartedOnBoundary_) return;

    const scalar span = t_.last() - t_.first();

    auto integral = [&](const DynamicList<scalar>& y, const bool positive, const bool negative)
    {
        scalar s = 0;
        for (label i = 1; i < t_.size(); ++i)
        {
            auto f = [&](const scalar q)
            {
                return positive ? Foam::max(q, scalar(0)) : negative ? Foam::min(q, scalar(0)) : q;
            };
            s += 0.5*(f(y[i]) + f(y[i-1]))*(t_[i] - t_[i-1]);
        }
        return s;
    };
    auto maxOf = [](const DynamicList<scalar>& y) { return Foam::max(y); };
    auto minOf = [](const DynamicList<scalar>& y) { return Foam::min(y); };

    DynamicList<std::string> names, units;
    DynamicList<scalar> vals;
    auto put = [&](const std::string& n, const std::string& u, const scalar v)
    {
        names.append(n); units.append(u); vals.append(v);
    };

    for (const item& it : items_)
    {
        const std::string n = it.name + ".";

        if (it.type == "chamber")
        {
            const DynamicList<scalar>& V = v_[slot(it.cVolume)];
            const DynamicList<scalar>& P = v_[slot(it.cPressure)];
            const scalar edv = maxOf(V), esv = minOf(V);
            put(n + "EDV", volumeLabel_, edv);
            put(n + "ESV", volumeLabel_, esv);
            put(n + "SV", volumeLabel_, edv - esv);
            put(n + "EF", "-", edv > 0 ? (edv - esv)/edv : 0);
            put(n + "CO", volumeLabel_ + "/min", (edv - esv)*60/period_);
            put(n + "Pmax", pressureLabel_, maxOf(P));
            put(n + "Pmin", pressureLabel_, minOf(P));

            // Pressure-volume loop of this cycle
            std::ofstream pv(dir_/(it.name + ".pvLoop.dat"));
            pv << std::setprecision(10) << "# cycle " << cycle_ << ": t " << t_.first() << " .. " << t_.last()
               << "\n# V P\n# " << volumeLabel_ << ' ' << pressureLabel_ << '\n';
            for (label i = 0; i < t_.size(); ++i) pv << V[i] << ' ' << P[i] << '\n';
        }
        else if (it.type == "pressure")
        {
            const DynamicList<scalar>& P = v_[slot(it.cPressure)];
            put(n + "mean", pressureLabel_, integral(P, false, false)/span);
            put(n + "systolic", pressureLabel_, maxOf(P));
            put(n + "diastolic", pressureLabel_, minOf(P));
            put(n + "pulse", pressureLabel_, maxOf(P) - minOf(P));
        }
        else
        {
            const DynamicList<scalar>& Q = v_[slot(it.cFlow)];
            const scalar fwd = integral(Q, true, false);
            const scalar back = integral(Q, false, true);
            put(n + "mean", flowLabel_, integral(Q, false, false)/span);
            put(n + "max", flowLabel_, maxOf(Q));
            put(n + "min", flowLabel_, minOf(Q));
            put(n + "netVolume", volumeLabel_, fwd + back);
            put(n + "fwdVolume", volumeLabel_, fwd);
            put(n + "backVolume", volumeLabel_, back);
            put(n + "regurgitantFraction", "-", fwd > 0 ? -back/fwd : 0);
        }
    }

    if (!headerWritten_)
    {
        out_ << "# per-cycle metrics, units of the case; cycle k is [k T, (k+1) T]\n# cycle tEnd";
        for (const auto& n : names) out_ << ' ' << n;
        out_ << "\n# - s";
        for (const auto& u : units) out_ << ' ' << (u.empty() ? "-" : u);
        out_ << '\n';
        headerWritten_ = true;
    }

    out_ << cycle_ << ' ' << t_.last();
    for (const scalar v : vals) out_ << ' ' << v;
    out_ << '\n';
    out_.flush();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circulationMetrics::add(const scalar t, const scalarField& values)
{
    const scalar tol = 1e-9*period_;
    const label k = static_cast<label>(std::floor((t + tol)/period_));

    if (cycle_ < 0)
    {
        startCycle(k);
        cycleStartedOnBoundary_ = std::abs(t - k*period_) <= tol;
    }
    else if (k > cycle_)
    {
        // The cycle ends with the sample on its closing boundary, which also
        // opens the next cycle
        if (std::abs(t - (cycle_ + 1)*period_) <= tol)
        {
            t_.append(t);
            forAll(monitored_, i) v_[i].append(values[monitored_[i]]);
        }

        finishCycle();
        startCycle(k);
        cycleStartedOnBoundary_ = std::abs(t - k*period_) <= tol;
    }

    t_.append(t);
    forAll(monitored_, i) v_[i].append(values[monitored_[i]]);
}

// ************************************************************************* //
