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
    initialiseTaylorGreenVortex

Description
    Utility that initialises the U and p field for the a Taylor-Green vortex
    case given by

        u(x, y, t) = e^{-2pi^2 t/Re} sin(pi x) cos(pi y)
        v(x, y, t) = -e^{-2pi^2 t/Re} cos(pi x) sin(pi y)
        p(x, y, t) = e^{-4pi^2 t/Re} (1/4) [cos(2pi x) + cos(2pi y)]

    Assuming time (t) is 0.0, the intial fields are

        u(x, y, 0) = sin(pi x) cos(pi y)
        v(x, y, 0) = -cos(pi x) sin(pi y)
        p(x, y, 0) = (1/4) [cos(2pi x) + cos(2pi y)] + pOffset

    By default, the cell values are the exact solution at the cell centres.
    With the -cellAverage option, they are the exact cell averages, computed
    with a high-order cell quadrature. A high-order (k-exact) finite volume
    scheme needs cell averages, because the gap between the average and the
    centre value is of second order in the cell size. The boundary values are
    the exact solution at the face centres in both cases.

 Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "fvMeshQuadrature.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

const scalar pi = constant::mathematical::pi;

const scalar pOffset = 0.5;
//const scalar pOffset = 0;

vector exactU(const point& c)
{
    return vector
    (
        Foam::sin(pi*c.x())*Foam::cos(pi*c.y()),
      - Foam::cos(pi*c.x())*Foam::sin(pi*c.y()),
        0
    );
}

scalar exactP(const point& c)
{
    return 0.25*(Foam::cos(2*pi*c.x()) + Foam::cos(2*pi*c.y())) + pOffset;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "cellAverage",
        "set the cell values to the exact cell averages"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

#ifdef OPENFOAM_COM
    const bool cellAverage = args.found("cellAverage");
#else
    const bool cellAverage = args.optionFound("cellAverage");
#endif

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
    if (cellAverage)
    {
        Info<< "Setting the cell values to the exact cell averages" << endl;

        // The cell quadrature is exact for polynomials of degree 6, so its
        // error is far below the discretisation errors of interest
        const fvMeshQuadrature quad(mesh, 6, 6);
        auto& cellQuadPoints = compactListListCRef(quad.cellQuadPoints());
        auto& cellQuadWeights = compactListListCRef(quad.cellQuadWeights());
        const scalarField& V = mesh.V();

        forAll(U, cellI)
        {
            vector UIntegral = Zero;
            scalar pIntegral = 0;

            forAll(cellQuadPoints[cellI], qI)
            {
                const scalar w = cellQuadWeights[cellI][qI];
                UIntegral += w*exactU(cellQuadPoints[cellI][qI]);
                pIntegral += w*exactP(cellQuadPoints[cellI][qI]);
            }

            U[cellI] = UIntegral/V[cellI];
            p[cellI] = pIntegral/V[cellI];
        }
    }
    else
    {
        const vectorField& C = mesh.C();

        forAll(U, cellI)
        {
            U[cellI] = exactU(C[cellI]);
            p[cellI] = exactP(C[cellI]);
        }
    }

    // Set the boundary fields

    forAll(mesh.boundary(), patchI)
    {
        const vectorField& Cf = mesh.C().boundaryField()[patchI];

        vectorField patchU(Cf.size());
        scalarField patchP(Cf.size());

        forAll(Cf, faceI)
        {
            patchU[faceI] = exactU(Cf[faceI]);
            patchP[faceI] = exactP(Cf[faceI]);
        }

        U.boundaryFieldRef()[patchI] == patchU;

        p.boundaryFieldRef()[patchI] == patchP;
    }

    // Write U and p
    Info<< nl << "Writing U and p to " << runTime.timeName() << endl;
    U.write();
    p.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
