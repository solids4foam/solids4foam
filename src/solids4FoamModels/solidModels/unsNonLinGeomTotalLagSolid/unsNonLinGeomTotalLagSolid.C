/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "unsNonLinGeomTotalLagSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, unsNonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


scalar unsNonLinGeomTotalLagSolid::residual(const volVectorField& D) const
{
    return
        gMax
        (
#ifdef OPENFOAMESIORFOUNDATION
            DimensionedField<double, volMesh>
#endif
            (
                mag(D.internalField() - D.prevIter().internalField())
               /max
                (
                    gMax
                    (
#ifdef OPENFOAMESIORFOUNDATION
                        DimensionedField<double, volMesh>
#endif
                        (
                            mag
                            (
                                D.internalField() - D.oldTime().internalField()
                            )
                        )
                    ),
                    SMALL
                )
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsNonLinGeomTotalLagSolid::unsNonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Ff_
    (
        IOobject
        (
            "Ff",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    Finvf_
    (
        IOobject
        (
            "Finvf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(Ff_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    Jf_
    (
        IOobject
        (
            "Jf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(Ff_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    nonLinear_(solidModelDict().lookupOrDefault<Switch>("nonLinear", true)),
    debug_(solidModelDict().lookupOrDefault<Switch>("debug", false)),
    K_
    (
        solidModelDict().lookupOrDefault<dimensionedScalar>
        (
            "K",
            dimensionedScalar("K", dimless/dimTime, 0)
        )
    ),
    relativeTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "solutionTolerance",
            solutionTol()
        )
    )
{
    DisRequired();

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        mechanical().interpolate(D(), pointD(), false);
        mechanical().grad(D(), pointD(), gradD());
        mechanical().grad(D(), pointD(), gradDf_);
        gradDD() = gradD() - gradD().oldTime();
        Ff_ = I + gradDf_.T();
        Finvf_ = inv(Ff_);
        Jf_ = det(Ff_);

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool unsNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    scalar initialResidual = 0;
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<vector>::debug = 0;
#else
    lduSolverPerformance solverPerfD;
    blockLduMatrix::debug = 0;
#endif
    scalar res = 1.0;
    scalar maxRes = 0;
    scalar curConvergenceTolerance = solutionTol();

    // Reset enforceLinear switch
    enforceLinear() = false;

    do
    {
#ifdef OPENFOAMESIORFOUNDATION
        if (SolverPerformance<vector>::debug)
#else
        if (blockLduMatrix::debug)
#endif
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        // Store previous iteration to allow under-relaxation and residual
        // calculation
        D().storePrevIter();

        // Construct momentum equation in total Lagrangian form where gradients
        // are calculated directly at the faces
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
          + rho()*g()
        );

        // Add damping
        if (K_.value() > SMALL)
        {
            DEqn += K_*rho()*fvm::ddt(D());
        }

        // Enforce linear to improve convergence
        if (enforceLinear())
        {
            // Replace nonlinear terms with linear
            // Note: the mechanical law could still be nonlinear
            DEqn +=
                fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
              - fvc::div(mesh().Sf() & sigmaf_);
        }

        // Under-relax the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
        const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif

        // Solve the system
        solverPerfD = DEqn.solve();

        // Under-relax displacement field
        relaxField(D(), iCorr);

        if (iCorr == 0)
        {
            initialResidual = mag(solverPerfD.initialResidual());
        }

        // Interpolate D to pointD
        mechanical().interpolate(D(), pointD(), false);

        // Update gradient of displacement
        mechanical().grad(D(), pointD(), gradD());
        mechanical().grad(D(), pointD(), gradDf_);

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        Ff_ = I + gradDf_.T();

        // Inverse of the deformation gradient
        Finvf_ = inv(Ff_);

        // Jacobian of the deformation gradient
        Jf_ = det(Ff_);

        // Check if outer loops are diverging
        if (nonLinear_ && !enforceLinear())
        {
            checkEnforceLinear(Jf_);
        }

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);

        // Calculate relative momentum residual
        res = residual(D());

        if (res > maxRes)
        {
            maxRes = res;
        }

        curConvergenceTolerance = maxRes*relativeTol_;

        if (curConvergenceTolerance < solutionTol())
        {
            curConvergenceTolerance = solutionTol();
        }

        if
        (
#ifdef OPENFOAMESIORFOUNDATION
            SolverPerformance<vector>::debug
#else
            blockLduMatrix::debug
#endif
         || (iCorr % infoFrequency()) == 0
         || res < curConvergenceTolerance
         || maxIterReached() == nCorr()
        )
        {
            Info<< "Corr " << iCorr << ", relative residual = " << res << endl;
        }

        if (maxIterReached() == nCorr())
        {
            maxIterReached()++;
        }
    }
    while (res > curConvergenceTolerance && ++iCorr < nCorr());

    // Velocity
    U() = fvc::ddt(D());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Print summary of residuals
    Info<< solverPerfD.solverName() << ": Solving for " << D().name()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerfD.initialResidual()
        << ", No outer iterations = " << iCorr << nl
        << " Max relative residual = " << maxRes
        << ", Relative residual = " << res
        << ", enforceLinear = " << enforceLinear() << endl;

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    if (nonLinear_ && enforceLinear())
    {
        return false;
    }

    return true;
}


tmp<vectorField> unsNonLinGeomTotalLagSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impKf_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& gradD = gradDf_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigmaf_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    if (enforceLinear())
    {
        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - n*pressure)
                  - (n & sigma)
                  + (n & (impK*gradD))
                )*rImpK
            )
        );
    }
    else
    {
        // Patch total deformation gradient inverse
        const tensorField& Finv = Finvf_.boundaryField()[patchID];

        // Patch total Jacobian
        const scalarField& J = Jf_.boundaryField()[patchID];

        // Patch unit normals (deformed configuration)
        const vectorField nCurrent(J*Finv.T() & n);

        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - nCurrent*pressure)
                  - (nCurrent & sigma)
                  + (n & (impK*gradD))
                )*rImpK
            )
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
