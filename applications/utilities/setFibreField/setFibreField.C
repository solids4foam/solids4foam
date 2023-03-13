/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    setFibreField

Description
    Set the initial fibre field for cardiac simulations.

    The utility currently sets the fibre field for the Land et al. (2015)
    ellipsoidal ventricle benchmarks.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Read the transmural distance field
    // This will be calculated using the given boundary conditions
    Info<< "Reading t" << endl;
    volScalarField t
    (
        IOobject
        (
            "t",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Solve diffusion equation to calculate the transmural distance (t)
    // Set fixedValue = 0.0 on the endocardium and fixedValue = 1.0 on the
    // epicardium, with zeroGradient on other patches
    Info<< "Solving the diffusion equation for t" << endl;
    fvScalarMatrix tEqn(fvm::laplacian(t));
    tEqn.solve();
    Info<< "Writing t" << endl;
    t.write();

    const dimensionedScalar dimL("1", dimLength, 1.0);

    const volScalarField rs("rs", 0.001*(7.0*dimL + 3.0*t));
    Info<< "Writing rs" << endl;
    rs.write();

    const volScalarField rl("rl", 0.001*(17.0*dimL + 3.0*t));
    Info<< "Writing rl" << endl;
    rl.write();

    const volScalarField alpha("alpha", 90.0*dimL - 180.0*t);
    Info<< "Writing alpha" << endl;
    alpha.write();

    // Convert alpha to radians
    const volScalarField alphaRadians
    (
        "alphaRadians", alpha*constant::mathematical::pi/180.0
    );
    Info<< "Writing alphaRadians" << endl;
    alphaRadians.write();

    // Unit Cartesian vectors
    const vector iHat(1, 0, 0);
    const vector jHat(0, 1, 0);
    const vector kHat(0, 0, 1);

    // Coordinate fields
    const volScalarField x("x", mesh.C().component(vector::X));
    const volScalarField y("y", mesh.C().component(vector::Y));
    const volScalarField z("z", mesh.C().component(vector::Z));
    z.write();

    //const volScalarField posX(pos(x));
    const volScalarField posY(pos(y));
    const volScalarField negY(1.0 - posY);
    const volScalarField posZ(pos(z));
    const volScalarField negZ(1.0 - posZ);

    // Ellipsoid parameters
    const volScalarField u
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -acos(z/rl)
    );
    Info<< "Writing u" << endl;
    u.write();

    const volScalarField v
    (
        IOobject
        (
            "v",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (negY - posY)*acos(x/(rs*sin(u)))
        //asin(y/(rs*sin(u)))
    );
    Info<< "Writing v" << endl;
    v.write();

    const volVectorField dxdu
    (
        IOobject
        (
            "dxdu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
      //   (1.0 - posY)*rs*cos(u)*cos(v)*iHat
      // - posY*rs*cos(u)*cos(v)*iHat
      // + (1.0 - posY)*rs*cos(u)*sin(v)*jHat
      // - posY*rs*cos(u)*sin(v)*jHat
        rs*cos(u)*cos(v)*iHat
      + rs*cos(u)*sin(v)*jHat
      - rl*sin(u)*kHat
    );
    Info<< "Writing dxdu" << endl;
    dxdu.write();

    const volVectorField dxduHat
    (
        IOobject
        (
            "dxduHat",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dxdu/mag(dxdu)
    );
    Info<< "Writing dxduHat" << endl;
    dxduHat.write();

    const volVectorField dxdv
    (
        IOobject
        (
            "dxdv",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       - rs*sin(u)*sin(v)*iHat
       + rs*sin(u)*cos(v)*jHat
       //  - (1.0 - posY)*rs*sin(u)*sin(v)*iHat
       // + posY*rs*sin(u)*sin(v)*iHat
       // + (1.0 - posY)*rs*sin(u)*cos(v)*jHat
       // - posY*rs*sin(u)*cos(v)*jHat
       // + (1.0 - posY)*rs*sin(u)*cos(v)*jHat
       // - posY*rs*sin(u)*cos(v)*jHat
    );
    Info<< "Writing dxdv" << endl;
    dxdv.write();

    const volVectorField dxdvHat
    (
        IOobject
        (
            "dxdvHat",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dxdv/mag(dxdv)
    );
    Info<< "Writing dxdvHat" << endl;
    dxdvHat.write();

    // Initialise the fibre field
    const word fibreFieldName("f0");
    volVectorField fibres
    (
        IOobject
        (
            fibreFieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dxduHat*sin(alphaRadians/dimL) + dxdvHat*cos(alphaRadians/dimL)
    );
    fibres /= mag(fibres);

    Info<< "Writing " << fibres.name() << endl;
    fibres.write();

    // Create sheet direction field
    // This is not defined for Lund et al. because the Guccione model is
    // transversly isotropic, so it does not matter what direction it points
    // as long as it is normal to the fibres
    volVectorField sheet
    (
        IOobject
        (
            "s0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless, vector::zero)
    );

    // Arbitrarily select z direction and remove the component the thre fibre
    // direction
    sheet = (I - sqr(fibres)) & kHat;
    sheet /= mag(sheet);

    // Remove component

    sheet /= mag(sheet);

    Info<< "Writing " << sheet.name() << endl;
    sheet.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
