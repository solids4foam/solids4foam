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

    Both methods first solve laplacian(t) = 0 for the transmural coordinate
    t, which needs fixedValue 0 (endocardium) and 1 (epicardium) conditions
    in 0/t and zeroGradient elsewhere. The method is selected in the
    optional dictionary system/setFibreFieldDict:

    \verbatim
    method      rossiLassila;   // default when the dictionary is absent
    nCorrectors 0;              // optional, non-orthogonal correctors for t
    \endverbatim

    rossiLassila (default)
        The fibre field f0 (cell and face values) of the Land et al. (2015)
        ellipsoidal ventricle benchmarks: the longitudinal direction is
        rotated about the transmural direction by an angle varying linearly
        from 90 degrees at the endocardium to -90 degrees at the epicardium,
        after Rossi, Lassila et al. The base normal is z.

    ellipsoid
        The fibre (f0), sheet (s0) and sheet-normal (n0) fields of a
        truncated ellipsoid with its long axis along x, following the
        rule-based definition of Arostica et al. (2025), Section 3.2,
        eq. (12)-(16), a port of the benchmark's reference implementation
        (cardiac_benchmark, compute_system()):

        1. with r_l(t), r_s(t) and alpha(t) interpolated linearly between
           their endocardial and epicardial values, the ellipsoidal
           coordinates of a point x = (x, y, z) are
               mu    = atan2(sqrt(y^2 + z^2)/r_s, x/r_l)
               theta = pi - atan2(z, -y),  or 0 where mu < 1e-7,
           and the local basis is
               e_mu    = dx/dmu    = (-r_l sin mu, r_s cos mu cos theta,
                                       r_s cos mu sin theta)
               e_theta = dx/dtheta = (0, -r_s sin mu sin theta,
                                      r_s sin mu cos theta),
           both normalised;
        2. f0 = sin(alpha) e_mu + cos(alpha) e_theta,
           n0 = e_mu x e_theta / |e_mu x e_theta|,
           s0 = f0 x n0 / |f0 x n0|.

        It needs these further entries (lengths in m, angles in degrees):

        \verbatim
        method      ellipsoid;
        rShortEndo  0.025;
        rLongEndo   0.09;
        rShortEpi   0.035;
        rLongEpi    0.097;
        alphaEndo   -60;
        alphaEpi    60;
        nCorrectors 3;
        \endverbatim

Author
    Philip Cardiff, UCD.
    ellipsoid method: Created by Claude Opus 5.5, supervised by Chanikya
    Valeti, UCD.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "transform.H"
#include "unitConversion.H"

// Rotation tensor for a rotation of omega about the axis a, by Rodrigues'
// formula. OpenFOAM.com and .org provide this as rotationAboutAxis() in
// transform.H; foam-extend does not, so it is spelled out here rather than
// the utility being excluded from that build
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


    // Unit vector, or the fallback where v vanishes (on the long axis)
    static inline vector unitOr(const vector& v, const vector& fallback)
    {
        const scalar m = mag(v);
        return (m > VSMALL) ? v/m : fallback;
    }

    // Geometry and fibre angles of the ellipsoid method
    struct ellipsoidParams
    {
        scalar rsEndo, rlEndo, rsEpi, rlEpi, alphaEndo, alphaEpi;
    };

    // Fibre, sheet and sheet-normal at a point with transmural coordinate t,
    // for the ellipsoid method
    static void ellipsoidBasis
    (
        const ellipsoidParams& p,
        const point& x,
        const scalar tIn,
        vector& f,
        vector& s,
        vector& n
    )
    {
#ifdef OPENFOAM_NOT_EXTEND
        const scalar pi = constant::mathematical::pi;
#else
        const scalar pi = mathematicalConstant::pi;
#endif

        // Bound t so that a slightly over- or under-shooting Laplace
        // solution near the boundaries does not extrapolate the radii
        const scalar t = min(max(tIn, scalar(0)), scalar(1));

        const scalar rl = p.rlEndo + (p.rlEpi - p.rlEndo)*t;
        const scalar rs = p.rsEndo + (p.rsEpi - p.rsEndo)*t;
        const scalar al = p.alphaEndo + (p.alphaEpi - p.alphaEndo)*t;

        const scalar a = sqrt(sqr(x.y()) + sqr(x.z()))/rs;
        const scalar b = x.x()/rl;
        const scalar mu = atan2(a, b);
        const scalar theta = (mu < 1e-7) ? 0.0 : pi - atan2(x.z(), -x.y());

        const scalar sMu = sin(mu);
        const scalar cMu = cos(mu);
        const scalar sTh = sin(theta);
        const scalar cTh = cos(theta);

        const vector eMu
        (
            unitOr(vector(-rl*sMu, rs*cMu*cTh, rs*cMu*sTh), vector(0, 1, 0))
        );
        const vector eTheta
        (
            unitOr(vector(0, -rs*sMu*sTh, rs*sMu*cTh), vector(0, 0, 1))
        );

        f = unitOr(sin(al)*eMu + cos(al)*eTheta, eTheta);
        n = unitOr(eMu ^ eTheta, vector(1, 0, 0));
        s = unitOr(f ^ n, eMu);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Optional settings: when system/setFibreFieldDict is absent, the
    // rossiLassila method is used, as before the dictionary was introduced
    IOdictionary dict
    (
        IOobject
        (
            "setFibreFieldDict",
            runTime.system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    const word method(dict.lookupOrDefault<word>("method", "rossiLassila"));

    if (method != "rossiLassila" && method != "ellipsoid")
    {
        FatalErrorIn(args.executable())
            << "Unknown method '" << method << "' in setFibreFieldDict. "
            << "Valid methods are rossiLassila and ellipsoid."
            << exit(FatalError);
    }

    Info<< "Method: " << method << nl << endl;

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
    const label nCorr = dict.lookupOrDefault<label>("nCorrectors", 0);
    for (label corr = 0; corr <= nCorr; ++corr)
    {
        fvScalarMatrix tEqn(fvm::laplacian(t));
        tEqn.solve();
    }
    Info<< "Writing t" << endl;
    t.write();

    if (method == "ellipsoid")
    {
        #include "setEllipsoidFibres.H"

        Info<< "\nEnd\n" << endl;

        return 0;
    }

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
