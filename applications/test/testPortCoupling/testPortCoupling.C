/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Application
    testPortCoupling

Description
    Drives a circulation network through its external-coupling interface the
    way a CFD solver does, without a CFD solver. The network must have a
    coupledDomain with provider external. Each step the "CFD" supplies the
    port flows, the network is advanced, and the pressures to impose at the
    patches are read back:

        for each time step [t0, t1]
            network.setPortForcing(port, t0, t1, c)   // Q over the step
            network.advance(t1)                       // 0D integration
            P = network.portPressure(port)            // BC for the 3D patch

    The flows are taken from a table (constant/portFlows.dat: t Qinlet
    Qoutlet, as written by the case's Allrun) and presented as the straight
    line between their values at t0 and t1, which is what a first-order
    coupling supplies. The result is written like circulationFoam's, so it can
    be compared with the table-driven run of the same circuit.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "circulationNetwork.H"
#include "circulationWaveform.H"

#include <fstream>
#include <iomanip>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary props
    (
        IOobject
        (
            "circulationProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    circulationNetwork network(props, runTime.value());
    network.summary(Info);

    // "CFD" flows: reuse the waveform reader for the two table columns
    circulationContext ctx(props);
    auto column = [&](const label c)
    {
        dictionary d, w;
        w.add("type", word("table"));
        w.add("file", fileName("constant/portFlows.dat"));
        w.add("column", c);
        d.add("waveform", w);
        return circulationWaveform(d, ctx, circulationContext::quantity::flow);
    };
    const circulationWaveform QinTable = column(1);
    const circulationWaveform QoutTable = column(2);

    const fileName outDir
    (
        runTime.globalPath()/"postProcessing"/"circulation"/runTime.timeName()
    );
    mkDir(outDir);
    std::ofstream log(outDir/"circulation.dat");
    log << std::setprecision(12) << "# testPortCoupling\n# t";
    for (const word& c : network.columns()) log << ' ' << c;
    log << '\n';

    scalarField values;
    auto writeRow = [&]()
    {
        network.outputValues(values);
        log << network.time();
        for (const scalar v : values) log << ' ' << v;
        log << '\n';
    };
    writeRow();

    while (runTime.run())
    {
        const scalar t0 = runTime.value();
        runTime++;
        const scalar t1 = runTime.value();

        const scalar qi0 = QinTable.value(t0), qi1 = QinTable.value(t1);
        const scalar qo0 = QoutTable.value(t0), qo1 = QoutTable.value(t1);

        network.setPortForcing("LV.inlet", t0, t1, {qi0, qi1 - qi0, 0, 0});
        network.setPortForcing("LV.outlet", t0, t1, {qo0, qo1 - qo0, 0, 0});

        network.advance(t1);

        // What the CFD solver would impose on its patches next
        const scalar Pmitral = network.portPressure("LV.inlet");
        const scalar Paortic = network.portPressure("LV.outlet");
        (void)Pmitral; (void)Paortic;

        writeRow();
    }

    Info<< "End: final patch pressures [Pa] mitral "
        << network.portPressure("LV.inlet") << ", aortic "
        << network.portPressure("LV.outlet") << endl;

    return 0;
}

// ************************************************************************* //
