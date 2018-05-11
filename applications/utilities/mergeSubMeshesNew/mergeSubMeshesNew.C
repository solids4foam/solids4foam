/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Merge multiple region meshes in the current case to create one default
    region mesh.

    The mesh is overwritten.

    A cellZone is created for each sub-region.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "mergePolyMesh.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::noParallel();
    argList::validArgs.append("regions");
    argList::validOptions.insert("noFunctionObjects", "");

    argList args(argc, argv, false, true);

    if (args.additionalArgs().size() < 2)
    {
        FatalError
            << "At least two regions should be specified"
            << abort(FatalError);
    }

    Info<< "Reading regions" << endl;

    List<word> regionNames(args.additionalArgs().size());
    forAll(regionNames, argI)
    {
        regionNames[argI] = word(args.additionalArgs()[argI]);
        Info<< "    " << regionNames[argI] << endl;
    }
    Info<< endl;

#   include "createTime.H"

    Info<< "Reading region meshes for time = " << runTime.timeName() << endl;

    PtrList<polyMesh> meshes(regionNames.size());
    forAll(regionNames, meshI)
    {
        meshes.set
        (
            meshI,
            new polyMesh
            (
                IOobject
                (
                    regionNames[meshI],
                    runTime.timeName(),
                    runTime
                )
            )
        );

        Info<< "    " << regionNames[meshI]
            << ", nCells: " << meshes[meshI].nCells() << endl;

        // If the subMesh has no cellZones then we will add one cellZone which
        // contains all the cells in the subMesh

        polyMesh& mesh = meshes[meshI];

        if (mesh.cellZones().size() == 0)
        {
            mesh.cellZones().setSize(1);

            labelList addr = labelList(mesh.nCells(), -1);

            forAll(addr, cI)
            {
                addr[cI] = cI;
            }

            // Remove trailing 'Subset' from the regionName, if it exists
            fileName czName = regionNames[meshI];
            const int s = czName.size();
            if (s > 6)
            {
                if
                (
                    czName[s - 1] == 't'
                 && czName[s - 2] == 'e'
                 && czName[s - 3] == 's'
                 && czName[s - 4] == 'b'
                 && czName[s - 5] == 'u'
                 && czName[s - 6] == 'S'
                )
                {
                    czName.resize(s - 6);
                }
            }

            Info<< "    Adding cellZone: " << czName << endl;

            // Add new cellZone
            mesh.cellZones().set
            (
                0,
                new cellZone
                (
                    czName,
                    addr,
                    0,
                    mesh.cellZones()
                )
            );

            mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
        }
    }

    // Write the first subMesh to the base mesh
    {
        Info<< nl << "Copying region " << regionNames[0]
            << " into the base mesh" << nl << endl;

        polyMesh& mesh = meshes[0];

        mesh.rename(polyMesh::defaultRegion);

        mesh.path() = fileName(mesh.rootPath()/mesh.caseName());
        mesh.setInstance("constant");
        mesh.polyMesh::write();
    }

    // Merge the other region meshes into the base mesh
    for (int meshI = 1; meshI < regionNames.size(); meshI++)
    {
        Info<< "Merging region " << regionNames[meshI]
            << " into base mesh" << nl << endl;

        mergePolyMesh mesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.timeName(),
                runTime
            )
        );

        mesh.addMesh(meshes[meshI]);
        mesh.merge();

        // Write the new merged mesh
        mesh.path() = fileName(mesh.rootPath()/mesh.caseName());
        mesh.setInstance("constant");
        mesh.polyMesh::write();
    }

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
