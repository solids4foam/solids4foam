/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    initialiseTaylorGreenVortex

Description
    Utility that initialises the U and p field for the a Taylor-Green vortex
    case given by

        u(x, y, t) = e^{-2pi^2 t/Re} sin(pi x) cos(pi y)
        v(x, y, t) = -e^{-2pi^2 t/Re} cos(pi x) sin(pi y)
        p(x, y, t) = e^{-4pi^2 t/Re} (1/4) [cos(2pi x) + sin(2pi y)]

    Assuming time (t) is 0.0, the intial fields are

        u(x, y, 0) = sin(pi x) cos(pi y)
        v(x, y, 0) = -cos(pi x) sin(pi y)
        p(x, y, 0) = (1/4) [cos(2pi x) + sin(2pi y)] + pOffset

 Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const scalar pi = constant::mathematical::pi;

    const scalar pOffset = 0.5;
    //const scalar pOffset = 0;

    // Initialise the velocity field
    Info<< "Creating U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    );

    // Initialise the pressure field
    Info<< "Creating p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, 0.0)
    );

    // Set the internal field
    {
        const scalarField x(mesh.C().component(vector::X));
        const scalarField y(mesh.C().component(vector::Y));

        U.primitiveFieldRef().replace
        (
            vector::X,
            Foam::sin(pi*x)*Foam::cos(pi*y)
        );

        U.primitiveFieldRef().replace
        (
            vector::Y,
          - Foam::cos(pi*x)*Foam::sin(pi*y)
        );

        p.primitiveFieldRef() =
            0.25*(Foam::cos(2*pi*x) + Foam::cos(2*pi*y)) + pOffset;
    }

    // Set the boundary fields

    forAll(mesh.boundary(), patchI)
    {
        const scalarField x
        (
            mesh.C().boundaryField()[patchI].component(vector::X)
        );
        const scalarField y
        (
            mesh.C().boundaryField()[patchI].component(vector::Y)
        );

        vectorField patchU
        (
            Foam::sin(pi*x)*Foam::cos(pi*y)*vector(1, 0, 0)
          - Foam::cos(pi*x)*Foam::sin(pi*y)*vector(0, 1, 0)
        );

        U.boundaryFieldRef()[patchI] == patchU;

        p.boundaryFieldRef()[patchI] ==
            0.25*(Foam::cos(2*pi*x) + Foam::cos(2*pi*y)) + pOffset;
    }

    // Write U and p
    Info<< nl << "Writing U and p to " << runTime.timeName() << endl;
    U.write();
    p.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
