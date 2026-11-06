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
    Test-fvcDdt

Description
    Checks the order of accuracy and the time-level handling of the BDF ddt
    scheme, for orders 1 to 6, on the mesh of the case:

    - Start-up and exactness: from a fresh start the order ramps 1, 2, ...,
      q. On a step of order m, fvc::ddt is exact for a time polynomial of
      degree m but not m + 1.
    - Convergence: dy/dt = lambda y, integrated with fvm::ddt from a supplied
      history, has an observed order of q as the step is halved.
    - Restart: stopping, writing and reading back the old-time fields
      reproduces an uninterrupted run.
    - A rho without stored history: fvc::ddt(rho, U) is exact for a rho
      linear in time, rather than weighting its copied levels.
    - Fluxes without stored history: the ddtCorr flux corrections keep the
      order of U and reconstruct the missing flux levels from U.

    The scheme is built for OpenFOAM.com only, so the test is skipped for
    other versions.

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_COM
    #include "fvCFD.H"
#else
    #include "IOstreams.H"
#endif

using namespace Foam;

#ifdef OPENFOAM_COM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

label failures = 0;

void check(const bool pass, const string& message)
{
    Info<< (pass ? "PASS: " : "FAIL: ") << message.c_str() << endl;

    if (!pass)
    {
        ++failures;
    }
}


//- The BDF scheme of order q, with an optional ddtPhiCoeff
template<class Type>
tmp<fv::ddtScheme<Type>> bdf
(
    const fvMesh& mesh,
    const label q,
    const scalar ddtPhiCoeff = -1
)
{
    OStringStream os;
    os  << "BDF " << q;
    if (ddtPhiCoeff >= 0)
    {
        os  << ' ' << ddtPhiCoeff;
    }

    IStringStream is(os.str());

    return fv::ddtScheme<Type>::New(mesh, is);
}


//- Start a new run at t0 with step dt. The time index keeps increasing, so
//  that fields created now start a fresh history.
void startRun(Time& runTime, const scalar t0, const scalar dt)
{
    runTime.setTime(t0, runTime.timeIndex());
    runTime.setDeltaT(dt, false);
}


scalar relativeError(const scalarField& a, const scalarField& b)
{
    return gMax(mag(a - b))/max(gMax(mag(b)), VSMALL);
}


volScalarField scalarField0(const fvMesh& mesh, const word& name)
{
    return volScalarField
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
}


//- Set y = exp(lambda t), the solution of dy/dt = lambda y with y(0) = 1
void setDecay(volScalarField& y, const scalar lambda)
{
    y.oldTime();
    y.primitiveFieldRef() = Foam::exp(lambda*y.time().value());
    y.correctBoundaryConditions();
}


//- Advance dy/dt = lambda y with fvm::ddt to tEnd. The matrix is diagonal,
//  so it is solved directly.
void advanceDecay
(
    Time& runTime,
    volScalarField& y,
    fv::ddtScheme<scalar>& scheme,
    const scalar lambda,
    const scalar tEnd
)
{
    const dimensionedScalar lambdaDim(inv(dimTime), lambda);

    while (runTime.value() < tEnd - 0.5*runTime.deltaTValue())
    {
        ++runTime;

        fvScalarMatrix eqn(scheme.fvmDdt(y));
        eqn -= fvm::Sp(lambdaDim, y);

        y.primitiveFieldRef() = eqn.source()/eqn.diag();
        y.correctBoundaryConditions();
    }
}


//- Start dy/dt = lambda y at t0 - (q - 1) dt and prescribe the exact
//  solution up to t0, so that the scheme has a full history at t0
void prescribeDecay
(
    Time& runTime,
    volScalarField& y,
    fv::ddtScheme<scalar>& scheme,
    const scalar lambda,
    const label q
)
{
    setDecay(y, lambda);

    for (label i = 1; i < q; ++i)
    {
        ++runTime;
        setDecay(y, lambda);

        // Create and fill the old levels
        scheme.fvcDdt(y);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void testStartUpAndExactness(Time& runTime, const fvMesh& mesh)
{
    Info<< nl << "Start-up and exactness" << endl;

    const scalarField magC(mag(mesh.C().primitiveField()));
    const scalarField f(1 + 0.5*magC/gMax(magC));
    const scalar dt = 0.1;

    for (label q = 1; q <= 6; ++q)
    {
        tmp<fv::ddtScheme<scalar>> scheme(bdf<scalar>(mesh, q));

        startRun(runTime, 0, dt);

        // y_d = (1 + t)^d f, for d = 1 to q + 1
        PtrList<volScalarField> y(q + 1);
        forAll(y, i)
        {
            y.set
            (
                i,
                new volScalarField
                (
                    scalarField0
                    (
                        mesh,
                        "yq" + Foam::name(q) + "d" + Foam::name(i + 1)
                    )
                )
            );
            y[i].primitiveFieldRef() = f;
        }

        bool rampOk = true;
        bool exactOk = true;

        for (label step = 1; step <= q + 2; ++step)
        {
            ++runTime;
            const scalar t = runTime.value();

            forAll(y, i)
            {
                y[i].oldTime();
                y[i].primitiveFieldRef() = Foam::pow(1 + t, i + 1)*f;
                y[i].correctBoundaryConditions();
            }

            // On this step the order should be m: exact for degree m, not
            // for degree m + 1
            const label m = min(step, q);

            forAll(y, i)
            {
                const label d = i + 1;
                const scalar error = relativeError
                (
                    scheme.ref().fvcDdt(y[i])().primitiveField(),
                    d*Foam::pow(1 + t, d - 1)*f
                );

                if (d <= m && error > 1e-10)
                {
                    (step <= q ? rampOk : exactOk) = false;
                    Info<< "    q = " << q << ", step " << step
                        << ": degree " << d << " error " << error << endl;
                }
                else if (d == m + 1 && error < 1e-6)
                {
                    rampOk = false;
                    Info<< "    q = " << q << ", step " << step
                        << ": degree " << d << " unexpectedly exact" << endl;
                }
            }
        }

        check
        (
            rampOk,
            "BDF " + Foam::name(q) + " start-up ramps the order 1 to "
          + Foam::name(q)
        );
        check
        (
            exactOk,
            "BDF " + Foam::name(q) + " is exact for degree " + Foam::name(q)
          + " polynomials"
        );
    }
}


void testConvergence(Time& runTime, const fvMesh& mesh)
{
    Info<< nl << "Convergence of dy/dt = lambda y, lambda = -1, to t = 1"
        << endl;

    const scalar lambda = -1;
    const scalar tEnd = 1;
    const scalarList dts({0.05, 0.025, 0.0125});

    for (label q = 1; q <= 6; ++q)
    {
        tmp<fv::ddtScheme<scalar>> scheme(bdf<scalar>(mesh, q));

        scalarList errors(dts.size());

        forAll(dts, j)
        {
            const scalar dt = dts[j];

            startRun(runTime, -(q - 1)*dt, dt);

            volScalarField y
            (
                scalarField0(mesh, "yConv" + Foam::name(q) + Foam::name(j))
            );

            prescribeDecay(runTime, y, scheme.ref(), lambda, q);
            advanceDecay(runTime, y, scheme.ref(), lambda, tEnd);

            errors[j] = gMax
            (
                mag
                (
                    y.primitiveField()
                  - Foam::exp(lambda*runTime.value())
                )
            );
        }

        const scalar observed =
            Foam::log(errors[1]/errors[2])/Foam::log(dts[1]/dts[2]);

        Info<< "    BDF " << q << ": errors " << errors
            << ", observed order " << observed << endl;

        check
        (
            mag(observed - q) < 0.2,
            "BDF " + Foam::name(q) + " observed order "
          + Foam::name(observed)
        );
    }
}


void testRestart(Time& runTime, const fvMesh& mesh)
{
    Info<< nl << "Restart" << endl;

    // Away from the times of the case, so as not to overwrite its results
    const scalar t0 = 100;
    const scalar lambda = -1;
    const scalar dt = 0.05;
    const label q = 3;

    tmp<fv::ddtScheme<scalar>> scheme(bdf<scalar>(mesh, q));

    // Uninterrupted
    startRun(runTime, t0 - (q - 1)*dt, dt);
    volScalarField yCont(scalarField0(mesh, "yCont"));
    prescribeDecay(runTime, yCont, scheme.ref(), lambda, q);
    advanceDecay(runTime, yCont, scheme.ref(), lambda, t0 + 1);

    // Stopped and written half way
    startRun(runTime, t0 - (q - 1)*dt, dt);
    fileName writeDir;
    {
        volScalarField y
        (
            IOobject
            (
                "yRestart",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
        prescribeDecay(runTime, y, scheme.ref(), lambda, q);
        advanceDecay(runTime, y, scheme.ref(), lambda, t0 + 0.5);

        // Write to full precision, so that the restart can be exact
        IOstream::defaultPrecision(17);
        runTime.writeNow();
        writeDir = runTime.timePath();
    }

    // Read back, with the old-time fields, and continue
    volScalarField y
    (
        IOobject
        (
            "yRestart",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    advanceDecay(runTime, y, scheme.ref(), lambda, t0 + 1);

    check
    (
        isFile(writeDir/"yRestart_0") && isFile(writeDir/"yRestart_0_0"),
        "BDF 3 writes the old-time fields yRestart_0 and yRestart_0_0"
    );

    const scalar error =
        relativeError(y.primitiveField(), yCont.primitiveField());
    check
    (
        error < 1e-12,
        "BDF 3 restart reproduces the uninterrupted run, relative "
        "difference " + Foam::name(error)
    );

    rmDir(writeDir);
}


void testRhoWithoutHistory(Time& runTime, const fvMesh& mesh)
{
    Info<< nl << "rho without stored history" << endl;

    const scalar dt = 0.1;
    const scalar c = 2;
    const label q = 3;

    tmp<fv::ddtScheme<vector>> scheme(bdf<vector>(mesh, q));

    startRun(runTime, 0, dt);

    volVectorField U
    (
        IOobject
        (
            "URho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, vector(1, 0, 0))
    );

    autoPtr<volScalarField> rhoPtr;
    scalar maxError = 0;

    for (label step = 1; step <= 7; ++step)
    {
        ++runTime;
        const scalar t = runTime.value();

        if (!rhoPtr)
        {
            scheme.ref().fvcDdt(U);

            // rho is created at the end of step 3, as a field derived rather
            // than read on a restart, with no old levels
            if (step == 3)
            {
                rhoPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "rhoNoHistory",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar(dimDensity, 1 + c*t)
                    )
                );
            }

            continue;
        }

        volScalarField& rho = rhoPtr.ref();

        // As a solver's first ddt would, before rho is updated
        rho.oldTime();
        rho.primitiveFieldRef() = 1 + c*t;
        rho.correctBoundaryConditions();

        const scalarField ddtx
        (
            scheme.ref().fvcDdt(rho, U)().primitiveField().component(vector::X)
        );

        maxError = max(maxError, gMax(mag(ddtx - c))/c);
    }

    check
    (
        maxError < 1e-10,
        "ddt(rho, U) is exact for a rho without history, relative error "
      + Foam::name(maxError)
    );
}


void testFluxWithoutHistory(Time& runTime, const fvMesh& mesh)
{
    Info<< nl << "Fluxes without stored history" << endl;

    // With ddtPhiCoeff 1 and each genuine flux level offset from the flux of
    // U by eps |Sf|, ddtCorr/(eps |Sf|/dt) is the sum of -a_k over the
    // genuine flux levels, as reconstructed levels contribute nothing: for
    // BDF 3, 3 and 3 - 3/2, then 11/6 once all are genuine
    const scalar dt = 0.1;
    const scalar eps = 1e-3;
    const label q = 3;
    const scalarList expected({3, 1.5, 11.0/6.0, 11.0/6.0});

    tmp<fv::ddtScheme<vector>> scheme(bdf<vector>(mesh, q, 1));

    const dimensionedScalar one(dimVelocity, 1);
    surfaceVectorField nf("nfTest", mesh.Sf()/mesh.magSf());
    nf.setOriented(false);

    startRun(runTime, 0, dt);

    volVectorField U
    (
        IOobject
        (
            "UFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, vector(1, 0, 0))
    );

    autoPtr<surfaceScalarField> phiPtr;
    autoPtr<surfaceVectorField> UfPtr;

    scalar maxPhiError = 0;
    scalar maxUfError = 0;

    for (label step = 1; step <= 7; ++step)
    {
        ++runTime;
        const scalar t = runTime.value();

        U.oldTime();
        U == dimensionedVector(dimVelocity, vector(1 + t, 0.5*t*t, 0));

        if (!phiPtr)
        {
            scheme.ref().fvcDdt(U);

            // The fluxes are created at the end of step 3 with no old levels
            if (step == 3)
            {
                phiPtr.reset
                (
                    new surfaceScalarField
                    (
                        "phiNoHistory",
                        mesh.Sf() & (fvc::interpolate(U) + eps*nf*one)
                    )
                );
                UfPtr.reset
                (
                    new surfaceVectorField
                    (
                        "UfNoHistory",
                        fvc::interpolate(U) + eps*nf*one
                    )
                );
            }

            continue;
        }

        surfaceScalarField& phi = phiPtr.ref();
        surfaceVectorField& Uf = UfPtr.ref();

        // As a solver's first ddt would, before the fluxes are updated
        phi.oldTime();
        Uf.oldTime();
        phi = mesh.Sf() & (fvc::interpolate(U) + eps*nf*one);
        Uf = fvc::interpolate(U) + eps*nf*one;

        const scalar e = expected[step - 4];

        const scalar rPhi =
            gAverage
            (
                (scheme.ref().fvcDdtPhiCorr(U, phi)/mesh.magSf())()
               .primitiveField()
            )*dt/eps;
        const scalar rUf =
            gAverage
            (
                (scheme.ref().fvcDdtUfCorr(U, Uf)/mesh.magSf())()
               .primitiveField()
            )*dt/eps;

        Info<< "    step " << step << ": phi " << rPhi << ", Uf " << rUf
            << ", expected " << e << endl;

        maxPhiError = max(maxPhiError, mag(rPhi - e));
        maxUfError = max(maxUfError, mag(rUf - e));
    }

    check
    (
        maxPhiError < 1e-8,
        "ddtCorr(U, phi) reconstructs the missing phi levels"
    );
    check
    (
        maxUfError < 1e-8,
        "ddtCorr(U, Uf) reconstructs the missing Uf levels"
    );
}

} // End anonymous namespace

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#ifdef OPENFOAM_COM

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    testStartUpAndExactness(runTime, mesh);
    testConvergence(runTime, mesh);
    testRestart(runTime, mesh);
    testRhoWithoutHistory(runTime, mesh);
    testFluxWithoutHistory(runTime, mesh);

    Info<< endl;

    if (failures)
    {
        FatalErrorInFunction
            << failures << " checks of the BDF ddt scheme failed"
            << abort(FatalError);
    }

    Info<< "Test-fvcDdt: PASSED" << endl;

#else

    Info<< "Test-fvcDdt: SKIPPED, the BDF ddt scheme is only built for "
        << "OpenFOAM.com" << endl;

#endif

    return 0;
}


// ************************************************************************* //
