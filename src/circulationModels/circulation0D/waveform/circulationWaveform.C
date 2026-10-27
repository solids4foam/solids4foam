/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "circulationWaveform.H"
#include "fileName.H"
#include "mathematicalConstants.H"
#include <cmath>

#include <fstream>
#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationWaveform::circulationWaveform
(
    const dictionary& cd,
    const circulationContext& ctx,
    const circulationContext::quantity q
)
:
    kind_(kind::constant),
    value_(0),
    amplitude_(0),
    period_(ctx.period()),
    phase_(0),
    systolicDuration_(0),
    periodic_(false)
{
    if (!cd.isDict("waveform"))
    {
        value_ = ctx.read(cd, "value", q);
        return;
    }

    const dictionary& d = cd.subDict("waveform");
    const word type = d.get<word>("type");

    if (type == "constant")
    {
        kind_ = kind::constant;
        value_ = ctx.read(d, "value", q);
    }
    else if (type == "sine")
    {
        kind_ = kind::sine;
        value_ = ctx.read(d, "mean", q, 0);
        amplitude_ = ctx.read(d, "amplitude", q);
        period_ = d.getOrDefault<scalar>("period", ctx.period());
        phase_ = d.getOrDefault<scalar>("phase", 0);
    }
    else if (type == "halfSine")
    {
        kind_ = kind::halfSine;
        amplitude_ = ctx.read(d, "amplitude", q);
        period_ = d.getOrDefault<scalar>("period", ctx.period());
        systolicDuration_ = d.get<scalar>("systolicDuration");
    }
    else if (type == "table")
    {
        kind_ = kind::table;
        periodic_ = d.getOrDefault<bool>("periodic", false);
        const label column = d.getOrDefault<label>("column", 1);

        fileName path(d.get<fileName>("file"));
        path.expand();

        std::ifstream is(path);
        if (!is)
        {
            FatalIOErrorInFunction(d)
                << "Cannot open waveform table '" << path << "'"
                << exit(FatalIOError);
        }

        DynamicList<scalar> tt, vv;
        std::string line;
        while (std::getline(is, line))
        {
            const auto hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);

            std::istringstream ls(line);
            scalarList row;
            double x;
            while (ls >> x) row.append(x);
            if (row.empty()) continue;

            if (row.size() <= column)
            {
                FatalIOErrorInFunction(d)
                    << "Table '" << path << "' has a row with " << row.size()
                    << " columns but column " << column << " was requested"
                    << exit(FatalIOError);
            }

            if (!tt.empty() && row[0] <= tt.last())
            {
                FatalIOErrorInFunction(d)
                    << "Times in table '" << path << "' must increase"
                    << exit(FatalIOError);
            }
            tt.append(row[0]);
            vv.append(row[column]*ctx.factor(q));
        }

        if (tt.size() < 2)
        {
            FatalIOErrorInFunction(d)
                << "Table '" << path << "' needs at least two rows"
                << exit(FatalIOError);
        }

        tTable_ = tt;
        vTable_ = vv;
    }
    else
    {
        FatalIOErrorInFunction(d)
            << "Unknown waveform type '" << type << "'. Valid types: constant"
            << " sine halfSine table" << exit(FatalIOError);
    }

    if (period_ <= 0)
    {
        FatalIOErrorInFunction(d) << "period must be positive"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::circulationWaveform::value(const scalar t) const
{
    switch (kind_)
    {
        case kind::constant:
            return value_;

        case kind::sine:
            return
                value_
              + amplitude_
               *Foam::sin(constant::mathematical::twoPi*t/period_ + phase_);

        case kind::halfSine:
        {
            const scalar tc = t - period_*std::floor(t/period_);
            return
                tc < systolicDuration_
              ? amplitude_*Foam::sin(constant::mathematical::pi*tc/systolicDuration_)
              : scalar(0);
        }

        case kind::table:
        {
            scalar tq = t;
            const scalar t0 = tTable_.first();
            const scalar t1 = tTable_.last();

            if (periodic_)
            {
                const scalar span = t1 - t0;
                tq = t0 + (t - t0) - span*std::floor((t - t0)/span);
            }

            if (tq <= t0) return vTable_.first();
            if (tq >= t1) return vTable_.last();

            // Binary search for the bracketing interval
            label lo = 0, hi = tTable_.size() - 1;
            while (hi - lo > 1)
            {
                const label mid = (lo + hi)/2;
                (tTable_[mid] <= tq ? lo : hi) = mid;
            }

            const scalar a = (tq - tTable_[lo])/(tTable_[hi] - tTable_[lo]);
            return (1 - a)*vTable_[lo] + a*vTable_[hi];
        }
    }

    return 0;
}

// ************************************************************************* //
