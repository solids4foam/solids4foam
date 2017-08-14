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
    setInitialTemperature

Description
    Simple utility to set the initial temperature field based on a mathematical
    funtion.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Read the temperature field
    Info<< "Reading T field from " << runTime.timeName() << nl << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Set temperature field
    Info<< "Setting T field to (500 + 6*y/sqrt(x^2+y^2))" << nl << endl;

    // X-coordinates
    const volScalarField x = mesh.C().component(vector::X);

    // Y-coordinates
    const volScalarField y = mesh.C().component(vector::Y);

    // Temperature parameters
    const dimensionedScalar a("a", dimTemperature, 500.0);
    const dimensionedScalar b("b", dimTemperature, 6.0);

    T = a + b*y/Foam::sqrt(pow(x, 2) + pow(y, 2));

    // Write the temperature field
    Info<< "Writing T field to " << runTime.timeName() << nl << endl;
    T.write();

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
