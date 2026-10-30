/*---------------------------------------------------------------------------*\
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Summary of mesh information (eg, after blockMesh)

    Copied unchanged from printMeshSummary.H of PDRblockMesh in
    OpenFOAM-v2306.

\*---------------------------------------------------------------------------*/

{
    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
        << "  " << "nPoints: " << mesh.nPoints() << nl
        << "  " << "nCells: " << mesh.nCells() << nl
        << "  " << "nFaces: " << mesh.nFaces() << nl
        << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

    const auto printZone =
        [](const Foam::zone& zn)
        {
            Info<< "  " << "zone " << zn.index()
                << " (size: " << zn.size()
                << ") name: " << zn.name() << nl;
        };

    if (mesh.cellZones().size())
    {
        Info<< "----------------" << nl
            << "Cell Zones" << nl
            << "----------------" << nl;

        for (const cellZone& zn : mesh.cellZones())
        {
            printZone(zn);
        }
    }
    if (mesh.faceZones().size())
    {
        Info<< "----------------" << nl
            << "Face Zones" << nl
            << "----------------" << nl;

        for (const faceZone& zn : mesh.faceZones())
        {
            printZone(zn);
        }
    }
    if (mesh.pointZones().size())
    {
        Info<< "----------------" << nl
            << "Point Zones" << nl
            << "----------------" << nl;

        for (const pointZone& zn : mesh.pointZones())
        {
            printZone(zn);
        }
    }

    Info<< "----------------" << nl
        << "Patches" << nl
        << "----------------" << nl;

    for (const polyPatch& p : mesh.boundaryMesh())
    {
        Info<< "  " << "patch " << p.index()
            << " (start: " << p.start()
            << " size: " << p.size()
            << ") name: " << p.name()
            << nl;
    }
}


// ************************************************************************* //
