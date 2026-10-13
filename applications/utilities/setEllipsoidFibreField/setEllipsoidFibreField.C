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
    setEllipsoidFibreField

Description
    Set the fibre (f0), sheet (s0) and sheet-normal (n0) fields of a truncated
    ellipsoidal ventricle, following the rule-based definition of
    Arostica et al. (2025), Section 3.2, eq. (12)-(16):

    1. The transmural coordinate t is the solution of
           laplacian(t) = 0,  t = 0 on the endocardium, t = 1 on the
           epicardium, zero gradient elsewhere (the base),
       which needs fixedValue 0 and 1 conditions on those patches in 0/t.

    2. With r_l(t), r_s(t) and alpha(t) interpolated linearly between their
       endocardial and epicardial values, the ellipsoidal coordinates of a
       point x = (x, y, z) (long axis along x) are
           mu    = atan2(sqrt(y^2 + z^2)/r_s, x/r_l)
           theta = pi - atan2(z, -y),  or 0 where mu < 1e-7,
       and the local basis is
           e_mu    = dx/dmu    = (-r_l sin mu, r_s cos mu cos theta,
                                   r_s cos mu sin theta)
           e_theta = dx/dtheta = (0, -r_s sin mu sin theta,
                                  r_s sin mu cos theta),
       both normalised.

    3. f0 = sin(alpha) e_mu + cos(alpha) e_theta,
       n0 = e_mu x e_theta / |e_mu x e_theta|,
       s0 = f0 x n0 / |f0 x n0|.

    This is a direct port of the reference implementation used by the
    benchmark (cardiac_benchmark, compute_system()).

    The parameters are read from system/setEllipsoidFibreFieldDict:
    \verbatim
    rShortEndo  0.025;
    rLongEndo   0.09;
    rShortEpi   0.035;
    rLongEpi    0.097;
    alphaEndo   -60;    // degrees
    alphaEpi    60;     // degrees
    nCorrectors 3;      // optional, non-orthogonal correctors for t
    \endverbatim

Author
    Written for solids4foam, after the Arostica et al. (2025) benchmark.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Unit vector, or the fallback where v vanishes (on the long axis)
    static inline vector unitOr(const vector& v, const vector& fallback)
    {
        const scalar m = mag(v);
        return (m > VSMALL) ? v/m : fallback;
    }

    struct ellipsoidParams
    {
        scalar rsEndo, rlEndo, rsEpi, rlEpi, alphaEndo, alphaEpi;
    };

    // Fibre, sheet and sheet-normal at a point with transmural coordinate t
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
        // Bound t so that a slightly over- or under-shooting Laplace
        // solution near the boundaries does not extrapolate the radii
        const scalar t = min(max(tIn, scalar(0)), scalar(1));

        const scalar rl = p.rlEndo + (p.rlEpi - p.rlEndo)*t;
        const scalar rs = p.rsEndo + (p.rsEpi - p.rsEndo)*t;
        const scalar al = p.alphaEndo + (p.alphaEpi - p.alphaEndo)*t;

        const scalar a = sqrt(sqr(x.y()) + sqr(x.z()))/rs;
        const scalar b = x.x()/rl;
        const scalar mu = atan2(a, b);
        const scalar theta =
            (mu < 1e-7) ? 0.0 : constant::mathematical::pi - atan2(x.z(), -x.y());

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


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set the fibre, sheet and sheet-normal fields of a truncated "
        "ellipsoid (Arostica et al. 2025)"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "setEllipsoidFibreFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    ellipsoidParams p;
    p.rsEndo = readScalar(dict.lookup("rShortEndo"));
    p.rlEndo = readScalar(dict.lookup("rLongEndo"));
    p.rsEpi = readScalar(dict.lookup("rShortEpi"));
    p.rlEpi = readScalar(dict.lookup("rLongEpi"));
    p.alphaEndo = degToRad(readScalar(dict.lookup("alphaEndo")));
    p.alphaEpi = degToRad(readScalar(dict.lookup("alphaEpi")));
    const label nCorr = dict.lookupOrDefault<label>("nCorrectors", 3);

    Info<< "Ellipsoid: endo (" << p.rsEndo << ", " << p.rlEndo << "), epi ("
        << p.rsEpi << ", " << p.rlEpi << "), alpha from "
        << radToDeg(p.alphaEndo) << " to " << radToDeg(p.alphaEpi)
        << " degrees" << nl << endl;

    // Transmural coordinate
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

    Info<< "Solving the Laplace equation for t" << endl;
    for (label corr = 0; corr <= nCorr; ++corr)
    {
        fvScalarMatrix tEqn(fvm::laplacian(t));
        tEqn.solve();
    }

    Info<< "    t in [" << gMin(t.primitiveField()) << ", "
        << gMax(t.primitiveField()) << "]" << endl;

    t.write();

    // Fibre, sheet and sheet-normal directions
    volVectorField f0
    (
        IOobject
        (
            "f0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero),
        "calculated"
    );
    volVectorField s0("s0", f0);
    volVectorField n0("n0", f0);

    const volVectorField& C = mesh.C();

    forAll(C, cellI)
    {
        ellipsoidBasis
        (
            p, C[cellI], t[cellI], f0[cellI], s0[cellI], n0[cellI]
        );
    }

#ifdef OPENFOAM_NOT_EXTEND
    auto& f0B = f0.boundaryFieldRef();
    auto& s0B = s0.boundaryFieldRef();
    auto& n0B = n0.boundaryFieldRef();
#else
    auto& f0B = f0.boundaryField();
    auto& s0B = s0.boundaryField();
    auto& n0B = n0.boundaryField();
#endif

    forAll(mesh.boundary(), patchI)
    {
        const vectorField& Cf = mesh.C().boundaryField()[patchI];
        const scalarField& tP = t.boundaryField()[patchI];

        forAll(Cf, faceI)
        {
            ellipsoidBasis
            (
                p,
                Cf[faceI],
                tP[faceI],
                f0B[patchI][faceI],
                s0B[patchI][faceI],
                n0B[patchI][faceI]
            );
        }
    }

    Info<< "Writing f0, s0 and n0" << endl;
    f0.write();
    s0.write();
    n0.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
