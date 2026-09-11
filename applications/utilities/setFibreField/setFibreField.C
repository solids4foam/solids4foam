/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

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

// Rotation tensor for a rotation of omega about the axis a, by Rodrigues'
// formula. OpenFOAM.com and .org provide this as rotationAboutAxis() in transform.H;
// foam-extend does not, so it is spelled out here rather than the utility
// being excluded from that build
namespace Foam
{
    static inline tensor rotationAboutAxis
    (
        const vector& a,
        const scalar omega
    )
    {
        const vector n = a/(mag(a) + VSMALL);
        const scalar c = cos(omega);
        const scalar s = sin(omega);

        return
            c*tensor::I
          + s*tensor(0, -n.z(), n.y(), n.z(), 0, -n.x(), -n.y(), n.x(), 0)
          + (1 - c)*(n*n);
    }
}

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
    const scalar alphaEndo = degToRad(90.0);

    // alpha at the epicardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEpi = degToRad(-90.0);

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
        const tensor rotT = rotationAboutAxis(et[cellI], alphaRadians[cellI]);
        f0[cellI] = transform(-rotT, f0[cellI]);
    }
    // foam-extend writes through boundaryField(); the others require
    // boundaryFieldRef()
#ifdef OPENFOAM_NOT_EXTEND
    auto& f0BoundaryRef = f0.boundaryFieldRef();
#else
    auto& f0BoundaryRef = f0.boundaryField();
#endif

    forAll(et.boundaryField(), patchI)
    {
        forAll(et.boundaryField()[patchI], faceI)
        {
            const tensor rotT = rotationAboutAxis
            (
                et.boundaryField()[patchI][faceI],
                alphaRadians.boundaryField()[patchI][faceI]
            );
            f0BoundaryRef[patchI][faceI] = transform
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
        const tensor rotT = rotationAboutAxis(etf[faceI], alphaRadiansf[faceI]);
        f0f[faceI] = transform(-rotT, f0f[faceI]);
    }
    // foam-extend writes through boundaryField(); the others require
    // boundaryFieldRef()
#ifdef OPENFOAM_NOT_EXTEND
    auto& f0fBoundaryRef = f0f.boundaryFieldRef();
#else
    auto& f0fBoundaryRef = f0f.boundaryField();
#endif

    forAll(etf.boundaryField(), patchI)
    {
        forAll(etf.boundaryField()[patchI], faceI)
        {
            const tensor rotT = rotationAboutAxis
            (
                etf.boundaryField()[patchI][faceI],
                alphaRadiansf.boundaryField()[patchI][faceI]
            );
            f0fBoundaryRef[patchI][faceI] = transform
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
