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
    setPoiseuilleFlow

Description
    Sets the initial fluid fields of the collapsibleChannel tutorial to the
    steady Poiseuille flow of unit mean velocity in the unit-height channel,

        u = 6 y (1 - y),    v = 0,
        p = 12 nu (25 - x) = 0.24 (25 - x),

    where p is the kinematic pressure for nu = 0.02, and 12 nu L = 6 is the
    inlet pressure that drives the flow through the channel of length
    L = 25.

    Only the internal fields are set; the boundary conditions are kept as
    read. It replaces runtime-compiled initial conditions, which OpenFOAM.org
    and foam-extend cannot compile or do not allow by default.

Usage
    setPoiseuilleFlow -region fluid

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    Info<< "Reading U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const vectorField& C = mesh.cellCentres();

    vectorField& UI = primitiveFieldRef(U);
    scalarField& pI = primitiveFieldRef(p);

    forAll(C, cellI)
    {
        const scalar x = C[cellI].x();
        const scalar y = C[cellI].y();

        UI[cellI] = vector(6.0*y*(1.0 - y), 0, 0);
        pI[cellI] = 0.24*(25.0 - x);
    }

    Info<< nl << "Writing U and p to " << runTime.timeName() << endl;
    U.write();
    p.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
