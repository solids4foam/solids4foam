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
    Test-circulationNetworkChecks

Description
    Checks of the lumped circulation network (circulation0D) on small
    circuits with known answers. No mesh or case is needed; the program
    writes the waveform table it reads into the working directory.

    1. A coupling port feeding an inductor: at the end of a step the port
       pressure is P_C + L dQ/dt, with dQ/dt from the forcing polynomial, for
       steps down to 1e-7 s.
    2. The derivatives of the sine, halfSine and table waveforms against
       central differences.
    3. A flow source feeding an inductor through a valve: the pressure of the
       node the inductor holds is P_C + L dQ/dt, for a constant and a sine
       flow.
    4. A flow source feeding a compliance through a hard valve (epsilon 0,
       reverseFraction 0), from an unpressurised and a pressurised compliance,
       with and without an inductor after the valve: the valve opens to
       dP = (Q/CQ)^2. A valve that should stay shut does, and a circuit with
       no solution is an error.

    5. A smooth valve (epsilon > 0) with a regurgitant leak: far from the
       switch its flow is that of the hard law, forward and backward.

    A line "PASS: ..." or "FAIL: ..." is printed for each check, and a
    summary at the end; the exit code is 1 if any check failed.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationNetwork.H"
#include "circulationWaveform.H"
#include "IStringStream.H"
#include "OFstream.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

label nPass = 0;
label nFail = 0;


void report(const bool ok, const string& what, const scalar value)
{
    if (ok)
    {
        ++nPass;
        Info<< "PASS: ";
    }
    else
    {
        ++nFail;
        Info<< "FAIL: ";
    }
    Info<< what.c_str() << " (" << value << ")" << endl;
}


dictionary readDict(const string& text)
{
    return dictionary(IStringStream(text)());
}


// 1. Port -> inductor -> compliance -> resistor -> ground
void portInductor()
{
    const scalar L = 2;
    const dictionary dict
    (
        readDict
        (
            "units SI;"
            "components {"
            "  ground { type ground; }"
            "  outlet { type externalPort; direction outlet; }"
            "  L { type inductor; L 2; }"
            "  C { type compliance; C 0.01; initialPressure 0; }"
            "  R { type resistor; R 10; }"
            "}"
            "connections ((outlet L C R ground));"
            "solver { type RKDP45; absTol 1e-14; relTol 1e-12;"
            "  maxSteps 100000000; }"
        )
    );

    // Q = c0 + c1 a + c2 a^2 + c3 a^3, a = t/dt
    externalPort::forcingCoeffs c;
    c[0] = 1;
    c[1] = 2;
    c[2] = -3;
    c[3] = 0.5;

    for (const scalar dt : {1e-2, 1e-3, 1e-6, 1e-7})
    {
        circulationNetwork net(dict, 0);
        net.setPortForcing("outlet", 0, dt, c);
        net.advance(dt);

        scalarField v;
        net.outputValues(v);
        const scalar PC = v[net.columns().find("C.P")];

        const scalar LdQ = L*(c[1] + 2*c[2] + 3*c[3])/dt;
        const scalar err = (net.portPressure("outlet") - PC - LdQ)/LdQ;

        report
        (
            mag(err) < 1e-8,
            "port -> inductor, dt " + Foam::name(dt)
          + ": (P_port - P_C - L dQ/dt)/(L dQ/dt)",
            err
        );
    }
}


// 2. Waveform derivatives against central differences
void waveformDerivatives()
{
    {
        OFstream os("circulationNetworkChecks.table");
        os  << "0 0" << nl << "0.1 2" << nl << "0.3 1" << nl << "0.6 4" << nl
            << "1 0" << nl;
    }

    const dictionary units(readDict("units SI;"));
    const circulationContext ctx(units);
    const dictionary wd
    (
        readDict
        (
            "sine { type sine; mean 1; amplitude 3; period 0.8; phase 0.4; }"
            "halfSine { type halfSine; amplitude 4; period 1;"
            "  systolicDuration 0.3; }"
            "table { type table; file \"circulationNetworkChecks.table\";"
            "  column 1; periodic yes; }"
        )
    );

    for (const word w : {"sine", "halfSine", "table"})
    {
        dictionary d;
        d.add("waveform", wd.subDict(w));
        const circulationWaveform f
        (
            d, ctx, circulationContext::quantity::flow
        );

        // Away from the kinks of halfSine and table
        scalar worst = 0, scale = 0;
        for (const scalar t : {0.05, 0.13, 0.27, 0.42, 0.71, 0.93, 1.37, 2.2})
        {
            const scalar h = 1e-7;
            const scalar fd = (f.value(t + h) - f.value(t - h))/(2*h);
            worst = max(worst, mag(f.derivative(t) - fd));
            scale = max(scale, mag(fd));
        }

        report
        (
            worst/scale < 1e-6,
            "waveform " + w + ": |derivative - central difference|/max",
            worst/scale
        );
    }
}


// 3. ground -> flowSource -> valve -> inductor -> C -> R -> ground
void pinnedRow(const word& title, const string& wave, const bool sine)
{
    const scalar L = 2;
    const dictionary dict
    (
        readDict
        (
            "units SI;"
            "components {"
            "  ground { type ground; }"
            "  S { type flowSource; waveform { " + wave + " } }"
            "  V { type orificeValve; CQ 2; epsilon 0.01; }"
            "  L { type inductor; L 2; }"
            "  C { type compliance; C 0.5; initialPressure 1; }"
            "  R { type resistor; R 3; }"
            "}"
            "connections ((ground S V L C R ground));"
            "probes (L.in);"
            "solver { type RKDP45; absTol 1e-12; relTol 1e-10; }"
        )
    );

    circulationNetwork net(dict, 0);
    const label iPC = net.columns().find("C.P");
    const label iPB = net.columns().find("L.in.P");

    scalar worst = 0;
    for (label i = 1; i <= 10; ++i)
    {
        const scalar t = 0.1*i;
        net.advance(t);

        scalarField v;
        net.outputValues(v);

        const scalar w = constant::mathematical::twoPi;
        const scalar dQ = sine ? 0.5*w*Foam::cos(w*t) : 0;
        worst = max(worst, mag(v[iPB] - v[iPC] - L*dQ));
    }

    report
    (
        worst < 1e-9,
        "valve -> inductor, " + title + ": max |P_B - P_C - L dQ/dt|",
        worst
    );
}


// 4. ground -> flowSource -> hard valve [-> inductor] -> C -> R -> ground
void valveOpening
(
    const bool inductor,
    const scalar P0,
    const scalar Q,
    const scalar CQ,
    const scalar C
)
{
    const string chain =
        inductor ? "ground S V L C R ground" : "ground S V C R ground";

    const dictionary dict
    (
        readDict
        (
            "units SI;"
            "components {"
            "  ground { type ground; }"
            "  S { type flowSource; waveform { type constant; value "
          + Foam::name(Q) + "; } }"
            "  V { type orificeValve; CQ " + Foam::name(CQ) + "; }"
          + (inductor ? "  L { type inductor; L 2; }" : "")
          + "  C { type compliance; C " + Foam::name(C)
          + "; initialPressure " + Foam::name(P0) + "; }"
            "  R { type resistor; R " + Foam::name(1/C) + "; }"
            "}"
            "connections ((" + chain + "));"
            "probes (V.in V.out);"
            "solver { type RKDP45; absTol 1e-14; relTol 1e-10; }"
        )
    );

    const string what =
        "hard valve opening, " + chain + ", P0 " + Foam::name(P0) + ", Q "
      + Foam::name(Q) + ", CQ " + Foam::name(CQ) + ", C " + Foam::name(C)
      + ": dP/(Q/CQ)^2 - 1";

    try
    {
        circulationNetwork net(dict, 0);
        net.advance(0.01);

        scalarField v;
        net.outputValues(v);
        const scalar dP =
            v[net.columns().find("V.in.P")] - v[net.columns().find("V.out.P")];
        const scalar err = dP/sqr(Q/CQ) - 1;

        report(mag(err) < 1e-8, what, err);
    }
    catch (const Foam::error& e)
    {
        report(false, what + ": " + e.message(), 0);
    }
}


// 4, continued: a source pulling flow back through the valve has no solution
void valveNoSolution()
{
    const dictionary dict
    (
        readDict
        (
            "units SI;"
            "components {"
            "  ground { type ground; }"
            "  S { type flowSource; waveform { type constant; value -1; } }"
            "  V { type orificeValve; CQ 2; }"
            "  C { type compliance; C 1; initialPressure 0; }"
            "  R { type resistor; R 1; }"
            "}"
            "connections ((ground S V C R ground));"
            "solver { type RKDP45; absTol 1e-14; relTol 1e-10; }"
        )
    );

    bool failed = false;
    try
    {
        circulationNetwork net(dict, 0);
        net.advance(0.01);
    }
    catch (const Foam::error&)
    {
        failed = true;
    }

    report
    (
        failed,
        "hard valve, a source pulling flow back through it: an error",
        failed
    );
}


// 4, continued: node A drains through Ra to C1 and through the valve to C2
void valveShut(const scalar P2, const bool shut)
{
    const dictionary dict
    (
        readDict
        (
            "units SI;"
            "components {"
            "  ground { type ground; }"
            "  S { type flowSource; waveform { type constant; value 1; } }"
            "  Ra { type resistor; R 2; }"
            "  C1 { type compliance; C 1; initialPressure 0; }"
            "  R1 { type resistor; R 1; }"
            "  V { type orificeValve; CQ 2; }"
            "  C2 { type compliance; C 1; initialPressure "
          + Foam::name(P2) + "; }"
            "  R2 { type resistor; R 1e6; }"
            "}"
            "connections ((ground S Ra C1 R1 ground) (S.out V C2 R2 ground));"
            "probes (S.out);"
            "solver { type RKDP45; absTol 1e-14; relTol 1e-10; }"
        )
    );

    circulationNetwork net(dict, 0);
    net.advance(0.01);

    scalarField v;
    net.outputValues(v);
    const scalar PA = v[net.columns().find("S.out.P")];
    const scalar QV = v[net.columns().find("V.Q")];
    const scalar QRa = (PA - v[net.columns().find("C1.P")])/2;

    if (shut)
    {
        report
        (
            QV == 0 && mag(QRa - 1) < 1e-8,
            "hard valve that stays shut, P_C2 " + Foam::name(P2)
          + ": valve flow 0, the source flow through Ra; Ra flow - 1",
            QRa - 1
        );
    }
    else
    {
        report
        (
            QV > 0 && mag(QV + QRa - 1) < 1e-8,
            "hard valve that opens, P_C2 " + Foam::name(P2)
          + ": the two flows add up to the source; sum - 1",
            QV + QRa - 1
        );
    }
}


// 5. pressureSource -> smooth leaky valve -> compliance -> R -> ground: the
//    valve flow at a fixed pressure drop, far from the switch
void smoothValveLeak(const scalar dP, const scalar fr)
{
    const scalar CQ = 350;
    const dictionary dict
    (
        readDict
        (
            "units SI;"
            "components {"
            "  ground { type ground; }"
            "  src { type pressureSource; waveform { type constant; value "
          + Foam::name(max(dP, scalar(0))) + "; } }"
            "  V { type orificeValve; CQ 350; epsilon 1; reverseFraction "
          + Foam::name(fr) + "; }"
            "  C { type compliance; C 1e6; initialPressure "
          + Foam::name(max(-dP, scalar(0))) + "; }"
            "  R { type resistor; R 1; }"
            "}"
            "connections ((src V C R ground));"
            "solver { type RKDP45; absTol 1e-14; relTol 1e-10; }"
        )
    );

    circulationNetwork net(dict, 0);
    net.advance(1e-9);

    scalarField v;
    net.outputValues(v);
    const scalar Q = v[net.columns().find("V.Q")];

    // The hard law at the same pressure drop
    const scalar Qhard =
        dP >= 0 ? CQ*Foam::sqrt(dP) : -fr*CQ*Foam::sqrt(-dP);
    const scalar err = mag(Q - Qhard)/(CQ*Foam::sqrt(mag(dP)));

    report
    (
        err < 1e-4,
        "smooth valve, epsilon 1, reverseFraction " + Foam::name(fr)
      + ", dP " + Foam::name(dP) + ": |Q - Q_hard|/(CQ sqrt|dP|)",
        err
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    FatalError.throwing(true);
    FatalIOError.throwing(true);

    portInductor();
    waveformDerivatives();

    pinnedRow("constant flow", "type constant; value 1;", false);
    pinnedRow
    (
        "sine flow",
        "type sine; mean 1; amplitude 0.5; period 1; phase 0;",
        true
    );

    for (const bool inductor : {false, true})
    {
        valveOpening(inductor, 0, 1, 2, 1);
        valveOpening(inductor, 1, 1, 2, 1);
        valveOpening(inductor, 0, 1e-4, 4.3e-6, 1e-8);
        valveOpening(inductor, 1.3e4, 1e-4, 4.3e-6, 1e-8);
        valveOpening(inductor, 0, 1, 2, 1e-8);
    }
    valveNoSolution();
    valveShut(100, true);
    valveShut(1.5, false);

    smoothValveLeak(100, 0.05);
    smoothValveLeak(-100, 0.05);
    smoothValveLeak(-100, 0);

    Info<< nl << "Circulation network checks: " << nPass << " passed, "
        << nFail << " failed" << endl;

    return (nFail == 0 ? 0 : 1);
}


// ************************************************************************* //
