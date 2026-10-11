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
    immersedBodyRefinementCells

Description
    Writes the cell set of the cells within a given distance of the immersed
    bodies of the immersedBoundaryForce finite volume options, over the motion
    of the bodies between two times, for refineMesh.

    The immersed boundary methods need about three cells or more across the
    thickness of a body, so that the cells whose centre is inside the body
    represent it: thin moving bodies, such as the leaflets of a heart valve,
    are resolved by refining the mesh around the region they sweep, one level
    at a time:

    \verbatim
    immersedBodyRefinementCells -distance 0.002
    refineMesh -overwrite
    immersedBodyRefinementCells -distance 0.001
    refineMesh -overwrite
    \endverbatim

    The surfaces of the bodies are sampled at the times where a vertex has
    moved by half the distance since the previous sample, so that the bands
    around the samples overlap. The cells whose centre is within the distance
    plus half the cell width of a sampled surface are selected.

Usage
    immersedBodyRefinementCells [OPTIONS]

    -distance <d>     Distance from the surfaces [m] (required)
    -set <name>       Name of the cell set (default refineCells)
    -bodies <list>    Names of the bodies (default all)
    -dict <file>      The fvOptions file (default constant/fvOptions, or
                      constant/fluid/fvOptions)
    -startTime <t0>   Start of the motion (default the start time)
    -endTime <t1>     End of the motion (default the end time)

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "immersedBody.H"
#include "cellSet.H"
#include "triSurfaceSearch.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("distance", "scalar", "distance from the surfaces");
    argList::addOption("set", "name", "name of the cell set");
    argList::addOption("bodies", "(name ...)", "names of the bodies");
    argList::addOption("dict", "file", "the fvOptions file");
    argList::addOption("startTime", "scalar", "start of the motion");
    argList::addOption("endTime", "scalar", "end of the motion");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    runTime.functionObjects().off();

    const scalar distance = args.get<scalar>("distance");
    const word setName(args.getOrDefault<word>("set", "refineCells"));
    wordList bodyNames;
    args.readListIfPresent("bodies", bodyNames);
    const scalar t0 =
        args.getOrDefault<scalar>("startTime", runTime.startTime().value());
    const scalar t1 =
        args.getOrDefault<scalar>("endTime", runTime.endTime().value());

    // The fvOptions file
    fileName dictFile;
    if (args.found("dict"))
    {
        dictFile = args.get<fileName>("dict");
    }
    else
    {
        dictFile = runTime.constant()/"fvOptions";
        if (!isFile(runTime.path()/dictFile))
        {
            dictFile = runTime.constant()/"fluid"/"fvOptions";
        }
    }

    const IOdictionary fvOptions
    (
        IOobject
        (
            dictFile,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // The immersed bodies of the immersedBoundaryForce options
    PtrList<immersedBody> bodies;
    for (const entry& e : fvOptions)
    {
        if
        (
            e.isDict()
         && e.dict().getOrDefault<word>("type", "") == "immersedBoundaryForce"
        )
        {
            for (const entry& be : e.dict().subDict("bodies"))
            {
                if
                (
                    be.isDict()
                 && (bodyNames.empty() || bodyNames.found(be.keyword()))
                )
                {
                    bodies.append
                    (
                        new immersedBody(be.keyword(), be.dict(), mesh)
                    );
                }
            }
        }
    }

    if (bodies.empty())
    {
        FatalErrorInFunction
            << "No immersed bodies found in " << dictFile << exit(FatalError);
    }

    const pointField& C = mesh.C();
    const scalarField halfWidth(0.5*cbrt(mesh.V().field()));

    cellSet cells(mesh, setName, mesh.nCells()/10 + 1);

    // Number of trial steps of the motion between the start and end times
    const label nSteps = 10000;

    for (immersedBody& body : bodies)
    {
        pointField lastPoints;
        label nSamples = 0;

        for (label stepi = 0; stepi <= nSteps; ++stepi)
        {
            const scalar t = t0 + (t1 - t0)*stepi/nSteps;

            body.move(t);
            const triSurface& surf = body.surface();

            // Sample the first and last configurations, and those where a
            // vertex has moved by half the distance since the last sample
            if (!lastPoints.empty() && stepi < nSteps)
            {
                if (max(mag(surf.points() - lastPoints)) < 0.5*distance)
                {
                    continue;
                }
            }

            lastPoints = surf.points();
            ++nSamples;

            const triSurfaceSearch search(surf);
            List<pointIndexHit> hits;
            search.findNearest
            (
                C,
                sqr(distance + halfWidth),
                hits
            );

            forAll(hits, celli)
            {
                if (hits[celli].hit())
                {
                    cells.insert(celli);
                }
            }

            if (!body.moving())
            {
                break;
            }
        }

        Info<< "Immersed body " << body.name() << ": " << nSamples
            << " samples of the surface between t = " << t0 << " and "
            << t1 << " s" << endl;
    }

    Info<< "Writing the " << cells.size() << " cells within " << distance
        << " m of the immersed bodies to the cell set " << setName
        << endl;

    cells.instance() = mesh.facesInstance();
    cells.write();

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
