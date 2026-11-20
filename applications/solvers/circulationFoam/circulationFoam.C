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
    circulationFoam

Description
    Solves a lumped-parameter (0D) model of the circulation described in
    constant/circulationProperties: named components joined by connections,
    integrated with an OpenFOAM ODE solver. No mesh is needed.

    controlDict supplies startTime, endTime, deltaT and the write times.
    deltaT is the output interval; the ODE solver takes whatever internal
    steps its tolerances need inside it.

    With -check nothing is solved: the circuit is described (elements, nodes,
    states, coupling ports, solver) and its set-up checked, including parameter
    magnitudes and, when the case has a mesh, that each coupled port's patch
    exists. Works on any case that has constant/circulationProperties, CFD
    cases included.

    Every deltaT a row is appended to
        postProcessing/circulation/<startTime>/circulation.dat
    with the pressure, volume and flow of every element in the units the case
    is written in. At write times the state is stored in
    <time>/uniform/circulationState-circulationProperties, and is read from the
    start time if present, so a long spin-up can be continued. The file is the
    restart checkpoint of the CFD coupler too, so a spin-up can start a coupled
    run (see circulationCoupler).

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "polyBoundaryMeshEntries.H"
#include "circulationNetwork.H"
#include "circulationMetrics.H"
#include "OSspecific.H"
#include "OFstream.H"

#include <limits>

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//- Patch names of constant/polyMesh/boundary, read without building the
//  mesh. Empty if the case has no mesh.
static wordList readPatchNames(const Time& runTime)
{
    IOobject boundaryIO
    (
        "boundary",
        runTime.constant(),
        polyMesh::meshSubDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    wordList names;
    if (boundaryIO.typeHeaderOk<polyBoundaryMeshEntries>(true))
    {
        const polyBoundaryMeshEntries entries(boundaryIO);
        names.setSize(entries.size());
        forAll(entries, patchi)
        {
            names[patchi] = entries[patchi].keyword();
        }
    }
    return names;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Lumped-parameter (0D) circulation network solver");
    argList::addBoolOption
    (
        "check",
        "Describe the circuit and check the set-up, then exit"
    );

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

    if (args.found("check"))
    {
        const wordList patches = readPatchNames(runTime);

        network.check(Info, patches.empty() ? nullptr : &patches);

        if (props.isDict("coupling"))
        {
            const dictionary& cd = props.subDict("coupling");
            Info<< nl << "CFD coupling:  scheme " << cd.get<word>("scheme")
                << ", order " << cd.get<label>("order")
                << ", rho " << cd.get<scalar>("rho") << nl;
        }

        return 0;
    }

    network.summary(Info);

    // The flow of an external port is supplied by a CFD solver, which
    // circulationFoam is not: warn that those ports carry no flow
    {
        DynamicList<word> externalPorts;

        for (const word& port : network.portNames())
        {
            if (network.portExternal(port))
            {
                externalPorts.append(port);
            }
        }

        if (!externalPorts.empty())
        {
            WarningInFunction
                << "The circuit has " << externalPorts.size()
                << " external coupling port(s) " << flatOutput(externalPorts)
                << ", whose flow is supplied by a CFD solver." << nl
                << "    circulationFoam supplies none, so each port carries"
                << " zero flow." << nl
                << "    Run the coupled case with the CFD solver, or give the"
                << " port a \"provider table\"."
                << endl;
        }
    }

    // The restart state, in the format of the CFD coupler's checkpoint
    const word stateName("circulationState-circulationProperties");

    // Continue from a stored state if there is one at the start time
    {
        IOobject stateIO
        (
            stateName,
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (stateIO.typeHeaderOk<IOdictionary>(true))
        {
            const IOdictionary state(stateIO);
            network.readState(state);

            // A state written at another time is taken at the start time
            if
            (
                mag(network.time() - runTime.value())
              > 1e-6*runTime.deltaTValue()
            )
            {
                Info<< "The circuit state of t = " << network.time()
                    << " is taken at the start time t = " << runTime.value()
                    << endl;
                network.setState(runTime.value(), network.state());
            }

            Info<< "Restarted from " << stateIO.objectPath() << nl << endl;
        }
    }

    const fileName outDir
    (
        runTime.globalPath()/"postProcessing"/"circulation"/runTime.timeName()
    );
    mkDir(outDir);

    OFstream log(outDir/"circulation.dat");
    log.precision(12);

    {
        const wordList& cols = network.columns();
        const stringList units = network.columnUnits();

        log << "# circulationFoam, units " << network.context().name() << nl;
        log << "# t";
        for (const word& c : cols)
        {
            log << ' ' << c;
        }
        log << "\n# s";
        for (const string& u : units)
        {
            log << ' ' << (u.empty() ? "-" : u.c_str());
        }
        log << '\n';
    }

    // Per-cycle metrics, if the case asks for them
    autoPtr<circulationMetrics> metrics;
    if (props.isDict("metrics"))
    {
        metrics.reset
        (
            new circulationMetrics
            (
                props.subDict("metrics"),
                network.columns(),
                network.context(),
                outDir
            )
        );
    }

    scalarField values;

    auto writeRow = [&]()
    {
        network.outputValues(values);
        log << network.time();
        for (const scalar v : values)
        {
            log << ' ' << v;
        }
        log << '\n';

        if (metrics.valid())
        {
            metrics().add(network.time(), values);
        }
    };

    writeRow();

    Info<< "Starting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        network.advance(runTime.value());
        writeRow();

        if (runTime.writeTime())
        {
            Info<< "Time = " << runTime.timeName() << endl;

            // The restart state at full double precision whatever
            // writePrecision the case uses. The values are converted to text
            // when the dictionary is built, so the precision is raised before
            // that.
            const unsigned int oldPrecision = IOstream::defaultPrecision
            (
                std::numeric_limits<scalar>::max_digits10
            );

            IOdictionary state
            (
                IOobject
                (
                    stateName,
                    runTime.timeName(),
                    "uniform",
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                network.stateDict()
            );
            state.regIOobject::writeObject
            (
                IOstreamOption(IOstreamOption::ASCII),
                true
            );

            IOstream::defaultPrecision(oldPrecision);

            log.flush();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
