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
    setPointPressure

Description
    Initialise the point pressure field (pointP) to have a constant gradient
    in the x direction.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Hard-code the pressure gradient
const scalar PRESSURE_GRADIENT = -4.2e-4;
//const scalar PRESSURE_GRADIENT = -4.2e-1;

// Hard-code the pressure origin: x location where pressure is zero
const scalar PRESSURE_ORIGIN = 0.2;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
    argList::noParallel();
#   include "createTime.H"
#   include "createMesh.H"

    // Create point mesh
    const pointMesh& pMesh(pointMesh::New(mesh));

    // Read point pressure field
    Info<< "Reading pointP" << nl << endl;
    pointScalarField pointP
    (
        IOobject
        (
            "pointP",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        pMesh
    );

    // Set pressure
    scalarField& pointPI = pointP;
    const pointField& points = mesh.points();
    forAll(pointPI, pI)
    {
        pointPI[pI] = PRESSURE_GRADIENT*(points[pI].x() - PRESSURE_ORIGIN);
    }
    pointP.correctBoundaryConditions();

    // Write pointP
    Info<< "Writing pointP to " << runTime.timeName() << endl;
    pointP.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
