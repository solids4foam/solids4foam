/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    smoothMesh

Description
    Smoothing mesh based on meshSmootherDict in system directory

Author
    Philip Cardiff UCD
    Peter De Jaeger
    Zeljko Tukovic

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshSmoother.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("latestTime", "");
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    runTime.functionObjects().off();

    if (args.optionFound("latestTime"))
    {
        instantList Times = runTime.times();
        runTime.setTime(Times[Times.size() - 1], Times.size() - 1);
        Info<< "Reading mesh from time: " << runTime.value() << endl;
    }

#   include "createMesh.H"

    IOdictionary meshSmootherDict
    (
        IOobject
        (
            "meshSmootherDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word smootherType =
        meshSmootherDict.lookupOrDefault<word>("type", "taubinSmoother");

    meshSmoother* smootherPtr = meshSmoother::New(smootherType, mesh).ptr();
    meshSmoother& smoother = *smootherPtr;

    smoother.smooth(false);
    smoother.writeSmoothedMesh(runTime);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
