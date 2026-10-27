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
    <time>/circulationState, and is read from the start time if present, so a
    long spin-up can be continued.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "circulationNetwork.H"
#include "circulationMetrics.H"
#include "OSspecific.H"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace Foam;

//- Patch names of constant/polyMesh/boundary, read without building the mesh
//  (a flat scan: a word followed by "{" inside the outer list). Empty if the
//  case has no mesh.
static wordList readPatchNames(const fileName& path)
{
    wordList names;
    std::ifstream is(path);
    if (!is)
    {
        return names;
    }

    std::stringstream buf;
    buf << is.rdbuf();
    const std::string text = buf.str();

    // Skip the FoamFile header block
    std::size_t pos = text.find("FoamFile");
    if (pos != std::string::npos)
    {
        pos = text.find('}', pos);
        if (pos != std::string::npos)
        {
            ++pos;
        }
        else
        {
            pos = 0;
        }
    }
    else
    {
        pos = 0;
    }

    // Tokens: whitespace separated, with ( ) { } as their own tokens
    std::vector<std::string> tok;
    std::string cur;
    for (std::size_t i = pos; i < text.size(); ++i)
    {
        const char c = text[i];
        if (c == '/' && i + 1 < text.size() && text[i+1] == '/')
        {
            while (i < text.size() && text[i] != '\n')
            {
                ++i;
            }
            continue;
        }
        if
        (
            std::isspace(static_cast<unsigned char>(c))
         || c == '(' || c == ')' || c == '{' || c == '}'
        )
        {
            if (!cur.empty())
            {
                tok.push_back(cur);
                cur.clear();
            }
            if (!std::isspace(static_cast<unsigned char>(c)))
            {
                tok.push_back(std::string(1, c));
            }
        }
        else
        {
            cur += c;
        }
    }
    if (!cur.empty())
    {
        tok.push_back(cur);
    }

    int depth = 0;
    for (std::size_t i = 0; i + 1 < tok.size(); ++i)
    {
        if (tok[i] == "(" )
        {
            ++depth;
            continue;
        }
        if (tok[i] == ")" )
        {
            --depth;
            continue;
        }
        if (tok[i] == "{" )
        {
            ++depth;
            continue;
        }
        if (tok[i] == "}" )
        {
            --depth;
            continue;
        }
        if (depth == 1 && tok[i+1] == "{")
        {
            names.append(word(tok[i], false));
        }
    }
    return names;
}

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
        const wordList patches =
            readPatchNames(runTime.path()/"constant"/"polyMesh"/"boundary");

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

    // Continue from a stored state if there is one at the start time
    {
        IOobject stateIO
        (
            "circulationState",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (stateIO.typeHeaderOk<IOdictionary>(true))
        {
            const IOdictionary state(stateIO);
            network.readState(state);
            Info<< "Restarted from " << stateIO.objectPath() << nl << endl;
        }
    }

    const fileName outDir
    (
        runTime.globalPath()/"postProcessing"/"circulation"/runTime.timeName()
    );
    mkDir(outDir);

    std::ofstream log(outDir/"circulation.dat");
    log << std::setprecision(12);

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
            log << ' ' << (u.empty() ? "-" : u);
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
                    "circulationState",
                    runTime.timeName(),
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

    log.close();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
