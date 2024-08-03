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
#include "transform.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

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

    // Base direction as per Rossi-Lassila et al. approach
    const vector k(0, 0, 1);

    // Transmural direction
    volVectorField et("et", fvc::grad(t));
    et /= mag(et);

    // Normal direction
    volVectorField en("en", k - (k & et)*et);
    en /= mag(en);

    // Longitudinal direction
    volVectorField el("el", en ^ et);
    el /= mag(el);

    Info<< "Writing et, en, el" << endl;
    et.write();
    en.write();
    el.write();


    // alpha at the endocardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEndo = degToRad(60.0);

    // alpha at the epicardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEpi = degToRad(-60.0);

    // Calculate alphaRadians
    const volScalarField alphaRadians
    (
        "alphaRadians",
        alphaEndo*(1.0 - t) + alphaEpi*t
    );
    Info<< "Writing alphaRadians" << endl;
    alphaRadians.write();

    // Rotate el about et by alphaRadians to get f
    volVectorField f0("f0", el);
    forAll(et, cellI)
    {
        const tensor rotT = Ra(et[cellI], alphaRadians[cellI]);
        f0[cellI] = transform(-rotT, f0[cellI]);
    }
    forAll(et.boundaryField(), patchI)
    {
        forAll(et.boundaryField()[patchI], faceI)
        {
            const tensor rotT = Ra
            (
                et.boundaryField()[patchI][faceI],
                alphaRadians.boundaryField()[patchI][faceI]
            );
            f0.boundaryFieldRef()[patchI][faceI] = transform
            (
                -rotT,
                f0.boundaryField()[patchI][faceI]
            );
        }
    }

    f0.write();


    // Calculate surface fields
    const surfaceScalarField tf(fvc::interpolate(t));

    // Transmural direction
    surfaceVectorField etf("etf", fvc::interpolate(fvc::grad(t)));
    //surfaceVectorField n(mesh.Sf()/mesh.magSf());
    // etf += fvc::snGrad(t)*n - ((I - sqr(n)) & etf);
    etf /= mag(etf);

    // Normal direction
    surfaceVectorField enf("enf", k - (k & etf)*etf);
    enf /= mag(enf);

    // Longitudinal direction
    surfaceVectorField elf("elf", enf ^ etf);
    elf /= mag(elf);

    etf.write();
    enf.write();
    elf.write();

    // Calculate alphaRadians
    const surfaceScalarField alphaRadiansf
    (
        "alphaRadiansf",
        alphaEndo*(1.0 - tf) + alphaEpi*tf
    );
    Info<< "Writing alphaRadiansf" << endl;
    alphaRadiansf.write();

    // Rotate el about et by alphaRadians to get f
    surfaceVectorField f0f("f0f", elf);
    forAll(etf, faceI)
    {
        const tensor rotT = Ra(etf[faceI], alphaRadiansf[faceI]);
        f0f[faceI] = transform(-rotT, f0f[faceI]);
    }
    forAll(etf.boundaryField(), patchI)
    {
        forAll(etf.boundaryField()[patchI], faceI)
        {
            const tensor rotT = Ra
            (
                etf.boundaryField()[patchI][faceI],
                alphaRadiansf.boundaryField()[patchI][faceI]
            );
            f0f.boundaryFieldRef()[patchI][faceI] = transform
            (
                -rotT,
                f0f.boundaryField()[patchI][faceI]
            );
        }
    }

    f0f.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
