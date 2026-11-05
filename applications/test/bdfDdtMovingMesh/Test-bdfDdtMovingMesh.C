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
    Test-bdfDdtMovingMesh

Description
    Checks the BDF ddt scheme of orders 1 to 6 on a moving mesh against the
    exact rate of change of a cell integral.

    The mesh of the case (a box along x) is moved every step by
    x = X + eps L sin(pi X/L) g(t), with g = sin(omega t) (default) or g = t
    (-poly), so that every cell's volume is known in closed form, and the
    field psi_i = f(t) (1 + X_i) is set from the exact function f each step,
    f = 1 + sin(2 t) or, with -poly, a polynomial of degree q - 1 for the
    order q under test. The finite-volume time derivative on a moving mesh
    is that of the cell integral, (1/V) d(V psi)/dt = psi' + psi V'/V, and
    BDFq must reproduce it to order q, exactly when V psi is a polynomial of
    degree q or less. Three things are compared, as maxima over the run
    (the derivative over its second half, after the order has ramped up),
    normalised by the maxima of the quantities themselves, for every order:

      - fvc::ddt(psi) against the exact derivative (relative max error);
      - fvm::ddt(psi) evaluated at psi, (diag psi - source)/V, against
        fvc::ddt(psi): the implicit and explicit forms must agree to
        round-off;
      - fvc::ddt(1) against div(meshPhi): the geometric conservation law,
        which the scheme's mesh flux must satisfy to round-off.

    The old-time chains of psi and of the scheme's volumes are created before
    the first step, as a solver's first call to fvm::ddt does, so the order
    ramps 1, 2, ..., q over the first steps and reaches q by the final step
    (the case has many more steps than 6).

    Usage: Test-bdfDdtMovingMesh [-deltaT dt] [-endTime T] [-poly]
    prints one line per order; see Allrun for the convergence study.

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_COM
    #include "fvCFD.H"
    #include "bdfDdtScheme.H"
    #include "fvcDiv.H"
    #include "IStringStream.H"
#else
    #include "IOstreams.H"
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#ifdef OPENFOAM_COM

    argList::addOption("deltaT", "scalar", "time step (default: controlDict)");
    argList::addOption("endTime", "scalar", "end time (default: controlDict)");
    argList::addBoolOption
    (
        "poly",
        "polynomial motion and fields: BDFq must be exact for order q"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.found("deltaT"))
    {
        runTime.setDeltaT(args.get<scalar>("deltaT"));
    }
    if (args.found("endTime"))
    {
        runTime.setEndTime(args.get<scalar>("endTime"));
    }
    const bool poly = args.found("poly");

    const scalar deltaT = runTime.deltaTValue();
    const label nSteps = round(runTime.endTime().value()/deltaT);

    // Motion x = X + eps L sin(pi X/L) g(t) along x; the box spans 0..L
    const pointField X(mesh.points());
    const scalar L = max(X.component(0)) - min(X.component(0));
    const scalar eps = 0.2;
    const scalar omega = 2*constant::mathematical::pi;

    auto g = [&](const scalar t) { return poly ? t : Foam::sin(omega*t); };
    auto dg = [&](const scalar t) { return poly ? 1.0 : omega*Foam::cos(omega*t); };
    auto shape = [&](const auto& x)
    {
        return eps*L*Foam::sin(constant::mathematical::pi*x/L);
    };

    // Per cell: the x extent of its points at rest, and the cross-section
    const scalarField V0(mesh.V());
    scalarField xlo(mesh.nCells(), GREAT), xhi(mesh.nCells(), -GREAT);
    forAll(mesh.cells(), celli)
    {
        for (const label pointi : mesh.cellPoints()[celli])
        {
            xlo[celli] = min(xlo[celli], X[pointi].x());
            xhi[celli] = max(xhi[celli], X[pointi].x());
        }
    }
    const scalarField Ayz(V0/(xhi - xlo));
    const scalarField Xc(mesh.C().primitiveField().component(0));

    // Exact cell volume and its rate
    auto Vexact = [&](const scalar t)
    {
        return Ayz*((xhi + shape(xhi)*g(t)) - (xlo + shape(xlo)*g(t)));
    };
    auto dVexact = [&](const scalar t)
    {
        return Ayz*(shape(xhi) - shape(xlo))*dg(t);
    };

    // The field: psi_i = f(t) (1 + X_i); f polynomial of degree q - 1 with
    // -poly (q from the order under test, so f is set per order below)
    auto f = [&](const scalar t, const label q)
    {
        if (poly)
        {
            scalar s = 0;
            for (label k = 0; k < q; ++k)
            {
                s += Foam::pow(1 + t, scalar(k));
            }
            return s;
        }
        return 1 + Foam::sin(2*t);
    };
    auto df = [&](const scalar t, const label q)
    {
        if (poly)
        {
            scalar s = 0;
            for (label k = 1; k < q; ++k)
            {
                s += k*Foam::pow(1 + t, scalar(k - 1));
            }
            return s;
        }
        return 2*Foam::cos(2*t);
    };

    const label maxOrder = 6;

    // One field per order, so that each has its own polynomial degree
    PtrList<volScalarField> psis(maxOrder);
    for (label q = 1; q <= maxOrder; ++q)
    {
        psis.set
        (
            q - 1,
            new volScalarField
            (
                IOobject("psi" + Foam::name(q), runTime.timeName(), mesh),
                mesh,
                dimensionedScalar(dimless, Zero)
            )
        );
        psis[q - 1].primitiveFieldRef() = f(0, q)*(1 + Xc);
        psis[q - 1].correctBoundaryConditions();

        // The old-time chain, as a solver's first fvm::ddt creates it
        const volScalarField* old = &psis[q - 1];
        for (label k = 0; k < maxOrder; ++k)
        {
            old = &old->oldTime();
        }
    }

    // Absolute maxima over the run, normalised at the end by the maxima of
    // the quantities themselves (the instantaneous ones pass through zero)
    scalarList ddtErr(maxOrder, 0.0), fvmErr(maxOrder, 0.0), gclErr(maxOrder, 0.0);
    scalarList ddtScale(maxOrder, 0.0), gclScale(maxOrder, 0.0);
    scalar geomErr = 0;

    for (label step = 1; step <= nSteps; ++step)
    {
        runTime++;
        const scalar t = runTime.value();

        pointField newPoints(X);
        forAll(newPoints, pointi)
        {
            newPoints[pointi].x() = X[pointi].x() + shape(X[pointi].x())*g(t);
        }
        mesh.movePoints(newPoints);

        const scalarField V(Vexact(t));
        geomErr = max(geomErr, gMax(mag(V - mesh.V()))/gMax(V));

        for (label q = maxOrder; q >= 1; --q)
        {
            volScalarField& psi = psis[q - 1];

            // The "solution" of this step; assigning it shifts the old levels
            psi.primitiveFieldRef() = f(t, q)*(1 + Xc);
            psi.correctBoundaryConditions();

            IStringStream is(Foam::name(q));
            fv::bdfDdtScheme<scalar> scheme(mesh, is);

            const scalarField exact
            (
                (df(t, q)*(1 + Xc)*V + f(t, q)*(1 + Xc)*dVexact(t))/V
            );

            // The errors are maxima over the second half of the run (the
            // order has ramped up by then at every step size used), so that
            // a zero crossing of the leading error term at one instant does
            // not hide it
            const bool measure = t > 0.5*runTime.endTime().value();

            const volScalarField ddt(scheme.fvcDdt(psi));
            ddtScale[q - 1] = max(ddtScale[q - 1], gMax(mag(exact)));
            if (measure)
            {
                ddtErr[q - 1] = max
                (
                    ddtErr[q - 1],
                    gMax(mag(ddt.primitiveField() - exact))
                );
            }

            const fvScalarMatrix M(scheme.fvmDdt(psi));
            const scalarField fvmDdt
            (
                (M.diag()*psi.primitiveField() - M.source())/mesh.V()
            );
            fvmErr[q - 1] = max
            (
                fvmErr[q - 1],
                gMax(mag(fvmDdt - ddt.primitiveField()))
            );

            const volScalarField ddtOne
            (
                scheme.fvcDdt(dimensionedScalar(dimless, 1))
            );
            const volScalarField divMeshPhi(fvc::div(scheme.meshPhi(psi)));
            gclScale[q - 1] =
                max(gclScale[q - 1], gMax(mag(divMeshPhi.primitiveField())));
            gclErr[q - 1] = max
            (
                gclErr[q - 1],
                gMax(mag(ddtOne.primitiveField() - divMeshPhi.primitiveField()))
            );
        }
    }

    for (label q = 1; q <= maxOrder; ++q)
    {
        ddtErr[q - 1] /= max(ddtScale[q - 1], SMALL);
        fvmErr[q - 1] /= max(ddtScale[q - 1], SMALL);
        gclErr[q - 1] /= max(gclScale[q - 1], SMALL);
    }

    Info<< "deltaT " << deltaT << " steps " << nSteps
        << " motion " << (poly ? "polynomial" : "smooth")
        << " geometry error " << geomErr << nl;

    for (label q = 1; q <= maxOrder; ++q)
    {
        Info<< "BDF " << q << ": ddt error " << ddtErr[q - 1]
            << " fvm-fvc " << fvmErr[q - 1]
            << " GCL " << gclErr[q - 1] << nl;
    }

    Info<< nl << "End" << nl << endl;

#else

    Info<< "Test-bdfDdtMovingMesh: the BDF ddt scheme requires OpenFOAM.com"
        << nl << endl;

#endif

    return 0;
}


// ************************************************************************* //
