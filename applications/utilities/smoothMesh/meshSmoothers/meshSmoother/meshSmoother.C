/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    meshSmoother

\*---------------------------------------------------------------------------*/

#include "meshSmoother.H"
#include "volFields.H"
#include "fvc.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(meshSmoother, 0);
defineRunTimeSelectionTable(meshSmoother, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSmoother::meshSmoother
(
    const word& name,
    fvMesh& mesh
)
:
    mesh_(mesh),
    meshSmootherDict_     //("meshSmootherDict")
    (
            IOobject
            (
                "meshSmootherDict",
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
    ),
    oldInstance_("constant")
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void meshSmoother::writeSmoothedMesh(Time& runTime)
{
    word writeOptionTag
    (
        dict().lookupOrDefault<word>
        (
            "writeOption",
            "overWrite"
        )
    );

    Info<< nl << "Writing mesh with option: " << writeOptionTag << nl << endl;

    bool writeMesh = false;

    if (writeOptionTag == "lastTimeStep")
    {
        Info<< "Writing the mesh to the previous time-step" << endl;
        mesh().setInstance(runTime.timeName());
        writeMesh = true;
    }
    else if (writeOptionTag == "overWrite")
    {
        Info<< "Overwriting the mesh" << endl;
        mesh().setInstance(oldInstance());
        writeMesh = true;
    }
    else if (writeOptionTag == "nextTimeStep")
    {
        Info<< "Writing the mesh to the next time-step" << endl;
        runTime++;
        mesh().setInstance(runTime.timeName());
        writeMesh = true;
    }
    else
    {
        FatalError
            << "Invalid option for meshSmoother::writeSmoothedMesh" << nl
            << "Possible valid options are: overWrite, lastTimeStep "
            << "and nextTimeStep. Default is overWrite." << nl
            << exit(FatalError);
    }

    if (writeMesh)
    {
        Info << "Writing mesh..." << nl << endl;
        mesh().write();
    }

}

} // End namespace Foam

// ************************************************************************* //
