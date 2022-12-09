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

#include "nonLinGeomTotalLagSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagSolid::nonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
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
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{
    DDisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(DD());

    // For consistent restarts, we will update the relative kinematic fields
    DD().correctBoundaryConditions();
    if (restart())
    {
        mechanical().grad(DD(), gradDD());
        gradD() = gradD().oldTime() + gradDD();
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector> solverPerfDD;
    SolverPerformance<vector>::debug = 0;
#else
    lduSolverPerformance solverPerfDD;
    blockLduMatrix::debug = 0;
#endif


    Info<< "Solving the total Lagrangian form of the momentum equation for DD"
        << endl;

    // Reset enforceLinear switch
    enforceLinear() = false;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Momentum equation incremental displacement total Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho(), DD())
          + fvc::d2dt2(rho().oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
          + rho()*g()
          + stabilisation().stabilisation(DD(), gradDD(), impK_)
        );

        // Enforce linear to improve convergence
        if (enforceLinear())
        {
            // Replace nonlinear terms with linear
            // Note: the mechanical law could still be nonlinear
            DDEqn +=
                fvc::div(J_*Finv_ & sigma(), "div(sigma)")
              - fvc::div(sigma());
        }

        // Under-relax the linear system
        DDEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DDEqn);

        // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
        const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field using fixed or adaptive under-relaxation
        relaxField(DD(), iCorr);

        // Update the total displacement
        D() = D().oldTime() + DD();

        // Update gradient of displacement increment
        mechanical().grad(DD(), gradDD());

        // Update the gradient of total displacement
        gradD() = gradD().oldTime() + gradDD();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Check if outer loops are diverging
        if (!enforceLinear())
        {
            checkEnforceLinear(J_);
        }

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DDEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Update impKf to improve convergence
        // Note: impK and rImpK are not updated as they are used for traction
        // boundaries
        // if (iCorr % 10 == 0)
        // {
        //     impKf_ = mechanical().impKf();
        // }
    }
    while
    (
       !converged
        (
            iCorr,
#ifdef OPENFOAMESIORFOUNDATION
            mag(solverPerfDD.initialResidual()),
            cmptMax(solverPerfDD.nIterations()),
#else
            solverPerfDD.initialResidual(),
            solverPerfDD.nIterations(),
#endif
            DD()
        ) && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD(), pointDD());

    // Increment of displacement
    D() = D().oldTime() + DD();

    // Increment of point displacement
    pointD() = pointD().oldTime() + pointDD();

    // Velocity
    U() = fvc::ddt(D());

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> nonLinGeomTotalLagSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

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
                  - (n & pSigma)
                  + impK*(n & pGradDD)
                )*rImpK
            )
        );
    }
    else
    {
        // Patch total deformation gradient inverse
        const tensorField& Finv = Finv_.boundaryField()[patchID];

        // Patch unit normals (deformed configuration)
        vectorField nCurrent(Finv.T() & n);
        nCurrent /= mag(nCurrent);

        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - nCurrent*pressure)
                  - (nCurrent & pSigma)
                  + impK*(n & pGradDD)
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
