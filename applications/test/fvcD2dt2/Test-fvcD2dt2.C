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
    Test-fvcD2dt2

Description
    Checks that the explicit second time derivative used by the solid models,
    fvcD2dt2Compat(), matches the implicit fvm::d2dt2 evaluated at the same
    field, for both the D and the (rho, D) forms, using the d2dt2 scheme
    of the case (#502).

    OpenFOAM's fvc::d2dt2 takes its scheme from ddtSchemes, whereas
    fvm::d2dt2 takes it from d2dt2Schemes; the mismatch of fvc::d2dt2 is
    reported for information only.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compatibilityFunctions.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar relativeDifference
(
    const GeometricField<Type, fvPatchField, volMesh>& a,
    const GeometricField<Type, fvPatchField, volMesh>& b
)
{
    const scalar scale =
        max(gMax(mag(primitiveField(a))), gMax(mag(primitiveField(b))));

    return
        gMax(mag(primitiveField(a) - primitiveField(b)))/max(scale, VSMALL);
}


// Return the matrix residual, M & D, i.e. the matrix evaluated at D
tmp<volVectorField> residual(const fvVectorMatrix& M, const volVectorField& D)
{
    // The foam-extend steadyState d2dt2 matrix is empty, without a diagonal
    if (!M.hasDiag())
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "residual",
                    D.time().timeName(),
                    D.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                D.mesh(),
                dimensionedVector
                (
                    "zero", M.dimensions()/dimVolume, vector::zero
                )
            )
        );
    }

    return M & D;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const scalar tolerance = 1e-10;

    const scalarField magC(mag(primitiveField(mesh.C())));

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1000)
    );
    primitiveFieldRef(rho) *= 1.0 + 0.5*magC/gMax(magC);

    volVectorField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Store two old-time levels
    D.oldTime().oldTime();

    // Advance a few time-steps with a displacement that is cubic in time, so
    // that every d2dt2 scheme gives a non-trivial result
    for (label timeI = 0; timeI < 3; timeI++)
    {
        runTime++;

        primitiveFieldRef(D) =
            pow3(runTime.value())*primitiveField(mesh.C());
        D.correctBoundaryConditions();
    }

    const fvVectorMatrix DEqn(fvm::d2dt2(D));
    const fvVectorMatrix rhoDEqn(fvm::d2dt2(rho, D));
    const volVectorField fvmD2dt2(residual(DEqn, D));
    const volVectorField fvmRhoD2dt2(residual(rhoDEqn, D));

    const scalar error = relativeDifference(fvcD2dt2Compat(D)(), fvmD2dt2);
    const scalar rhoError =
        relativeDifference(fvcD2dt2Compat(rho, D)(), fvmRhoD2dt2);

    Info<< "fvcD2dt2Compat(D) vs fvm::d2dt2(D): relative difference = "
        << error << nl
        << "fvcD2dt2Compat(rho, D) vs fvm::d2dt2(rho, D): relative "
        << "difference = " << rhoError << nl
        << "fvc::d2dt2(D) vs fvm::d2dt2(D) (information only): relative "
        << "difference = " << relativeDifference(fvc::d2dt2(D)(), fvmD2dt2)
        << endl;

    if (error > tolerance || rhoError > tolerance)
    {
        FatalErrorInFunction
            << "The explicit and implicit d2dt2 operators do not agree"
            << abort(FatalError);
    }

    Info<< "Test-fvcD2dt2: PASSED" << endl;

    return 0;
}


// ************************************************************************* //
