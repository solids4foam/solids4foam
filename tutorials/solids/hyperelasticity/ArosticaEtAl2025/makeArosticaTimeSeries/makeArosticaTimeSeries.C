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

Application
    makeArosticaTimeSeries

Description
    Tabulate the time series of the Arostica et al. (2025) cardiac
    elastodynamics benchmarks, as OpenFOAM interpolation tables:

      - the active stress tau(t), eq. (5)-(6) of the paper (Bestel model);
      - a cavity pressure p(t), eq. (7)-(8), with the parameters of Table 4
        (benchmark 1) or Table 9 (benchmark 2, LV or RV).

    Each is a scalar ODE in time only, y' = -c(t) y + s(t), y(0) = 0. Over a
    sub-step h with c and s frozen at its mid-point the exact update is
        y(t + h) = y_eq + (y(t) - y_eq) exp(-c h),   y_eq = s/c,
    which is unconditionally stable. With h = 1e-5 s it agrees with a
    tightly toleranced Radau solve to better than 1e-7 relative.

Usage
    makeArosticaTimeSeries -active -sigma0 1.5e5 [-tSys 0.16] [-tDias 0.484]
        -out constant/activeStress
    makeArosticaTimeSeries -pressure benchmark1|bivLV|bivRV
        -out constant/timeVsPressure
    Add -zero for an all-zero table, and -dt to set the table spacing
    (default 1e-4 s). The tables cover t in [0, 1] s.

Author
    Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.
    A port of the earlier makeTimeSeries.py of this tutorial.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "scalarField.H"

using namespace Foam;

// OpenFOAM.com (which defines OPENFOAM) and OpenFOAM.org name the option
// accessors of argList differently
#ifdef OPENFOAM
    #define ARG_FOUND(args, opt) args.found(opt)
    #define ARG_GET(T, args, opt) args.get<T>(opt)
    #define ARG_GET_OR(T, args, opt, def) args.getOrDefault<T>(opt, def)
#else
    #define ARG_FOUND(args, opt) args.optionFound(opt)
    #define ARG_GET(T, args, opt) args.optionRead<T>(opt)
    #define ARG_GET_OR(T, args, opt, def) args.optionLookupOrDefault<T>(opt, def)
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

// Smoothed step, 0.5*(1 + sign*tanh(dt/gamma))
scalar S(const scalar dt, const scalar gamma, const scalar sign)
{
    return 0.5*(1.0 + sign*::tanh(dt/gamma));
}

// Integrate y' = -c(t) y + s(t), y(0) = 0, over [0, tEnd] with sub-step h
// and sample it at tOut, by linear interpolation of the sub-step values
template<class CoeffFunction>
scalarField integrate
(
    const CoeffFunction& coeffs,
    const scalar tEnd,
    const scalar h,
    const scalarField& tOut
)
{
    const label n = label(::round(tEnd/h));
    scalarField y(n + 1, 0.0);

    for (label i = 0; i < n; ++i)
    {
        // Mid-point of the sub-step
        const scalar tm = (i + 0.5)*h;
        scalar c = 0;
        scalar s = 0;
        coeffs(tm, c, s);

        if (c > 0)
        {
            const scalar yEq = s/c;
            y[i + 1] = yEq + (y[i] - yEq)*::exp(-c*h);
        }
        else
        {
            // y' = s when c == 0 (not reached for these parameters)
            y[i + 1] = y[i] + s*h;
        }
    }

    scalarField result(tOut.size(), 0.0);
    forAll(tOut, j)
    {
        const scalar r = min(max(tOut[j]/h, scalar(0)), scalar(n));
        const label i = min(label(r), n - 1);
        const scalar w = r - i;
        result[j] = (1.0 - w)*y[i] + w*y[i + 1];
    }

    return result;
}

// Active stress, eq. (5)-(6), Table 3 (benchmark 1) or Table 8
struct activeStressCoeffs
{
    scalar sigma0, tSys, tDias, gamma, aMax, aMin;

    void operator()(const scalar t, scalar& c, scalar& s) const
    {
        const scalar f = S(t - tSys, gamma, 1)*S(t - tDias, gamma, -1);
        const scalar a = aMax*f + aMin*(1.0 - f);
        c = mag(a);
        s = sigma0*max(a, scalar(0));
    }
};

// Cavity pressure, eq. (7)-(8), Table 4 (benchmark 1) or Table 9
struct pressureCoeffs
{
    scalar sigmaPre, sigmaMid, alphaPre, alphaMid, tSysPre, tDiasPre,
        gamma, aMax, aMin;

    void operator()(const scalar t, scalar& c, scalar& s) const
    {
        const scalar f = S(t - tSysPre, gamma, 1)*S(t - tDiasPre, gamma, -1);
        const scalar a = aMax*f + aMin*(1.0 - f);
        const scalar g = S(t - tDiasPre, gamma, -1);
        const scalar b = a + alphaPre*g + alphaMid;
        c = mag(b);
        s = sigmaMid*max(b, scalar(0)) + sigmaPre*max(g, scalar(0));
    }
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::addBoolOption("active", "tabulate the active stress tau(t)");
    argList::addOption
    (
        "pressure", "name", "tabulate a pressure: benchmark1, bivLV or bivRV"
    );
    argList::addOption("sigma0", "scalar", "contractility sigma_0 [Pa]");
    argList::addOption("tSys", "scalar", "t_sys of the active stress [s]");
    argList::addOption("tDias", "scalar", "t_dias of the active stress [s]");
    argList::addOption("dt", "scalar", "table spacing [s] (default 1e-4)");
    argList::addBoolOption("zero", "write an all-zero table");
    argList::addOption("out", "file", "output file");

    argList args(argc, argv, false, false, false);

    if (ARG_FOUND(args, "active") == ARG_FOUND(args, "pressure"))
    {
        FatalErrorIn(args.executable())
            << "Give exactly one of -active and -pressure" << exit(FatalError);
    }
    if (!ARG_FOUND(args, "out"))
    {
        FatalErrorIn(args.executable())
            << "Give the output file with -out" << exit(FatalError);
    }

    const scalar dtOut = ARG_GET_OR(scalar, args, "dt", 1e-4);
    const label nOut = label(::round(1.0/dtOut));
    scalarField t(nOut + 1);
    forAll(t, i)
    {
        t[i] = i*dtOut;
    }

    const bool zero = ARG_FOUND(args, "zero");
    scalarField y(t.size(), 0.0);
    word what;
    string params;

    // Integration sub-step and end time
    const scalar h = 1e-5;
    const scalar tEnd = 1.0;

    if (ARG_FOUND(args, "active"))
    {
        activeStressCoeffs a;
        a.sigma0 = ARG_GET_OR(scalar, args, "sigma0", 1.5e5);
        a.tSys = ARG_GET_OR(scalar, args, "tSys", 0.16);
        a.tDias = ARG_GET_OR(scalar, args, "tDias", 0.484);
        a.gamma = 0.005;
        a.aMax = 5.0;
        a.aMin = -30.0;

        what = "Active stress tau(t), eq. (5)";
        params =
            zero ? string("zero")
          : string
            (
                "sigma0 = " + Foam::name(a.sigma0) + " Pa, t_sys = "
              + Foam::name(a.tSys) + " s, t_dias = " + Foam::name(a.tDias)
              + " s"
            );

        if (!zero)
        {
            y = integrate(a, tEnd, h, t);
        }
    }
    else
    {
        const word preset(ARG_GET(word, args, "pressure"));

        // Table 4 (benchmark 1) and Table 9 (benchmark 2)
        pressureCoeffs p;
        p.tSysPre = 0.17;
        p.tDiasPre = 0.484;
        p.gamma = 0.005;
        p.aMax = 5.0;
        p.aMin = -30.0;

        if (preset == "benchmark1")
        {
            p.sigmaPre = 7000; p.sigmaMid = 16000;
            p.alphaPre = 5; p.alphaMid = 1;
        }
        else if (preset == "bivLV")
        {
            p.sigmaPre = 12000; p.sigmaMid = 16000;
            p.alphaPre = 5; p.alphaMid = 15;
        }
        else if (preset == "bivRV")
        {
            p.sigmaPre = 3000; p.sigmaMid = 4000;
            p.alphaPre = 1; p.alphaMid = 10;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Unknown pressure " << preset << ": valid are benchmark1, "
                << "bivLV and bivRV" << exit(FatalError);
        }

        what = "Cavity pressure p(t), eq. (7)";
        params = zero ? string("zero") : string("preset " + preset);

        if (!zero)
        {
            y = integrate(p, tEnd, h, t);
        }
    }

    // Write the table
    const fileName out(ARG_GET(fileName, args, "out"));
    {
        OFstream os(out);
        os  << "// " << what.c_str() << nl
            << "// Generated by makeArosticaTimeSeries (" << params.c_str()
            << ")" << nl
            << "// Columns: time [s]   value [Pa]" << nl
            << "(" << nl;

        forAll(t, i)
        {
            os  << "    (" << setprecision(6) << t[i] << " "
                << setprecision(10) << y[i] << ")" << nl;
        }

        os  << ")" << nl;
    }

    label iMax = 0;
    forAll(y, i)
    {
        if (y[i] > y[iMax])
        {
            iMax = i;
        }
    }

    Info<< out << ": max " << y[iMax] << " Pa at t = " << t[iMax] << " s"
        << endl;

    return 0;
}


// ************************************************************************* //
