/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "coupledLinGeomPressureDisplacementSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "linearElastic.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#ifdef FOAMEXTEND
    #include "fvBlockMatrix.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledLinGeomPressureDisplacementSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, coupledLinGeomPressureDisplacementSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, coupledLinGeomPressureDisplacementSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledLinGeomPressureDisplacementSolid::coupledLinGeomPressureDisplacementSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    k_(mechanical().bulkModulus()),
    rK_(1.0/k_),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    gradp_(fvc::grad(p_)),
#ifdef FOAMEXTEND
    Dp_
    (
        IOobject
        (
            "Dp",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector4("zero", dimless, vector4::zero)
    ),
#endif
    pressureRhieChowScaleFac_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "pressureRhieChowScaleFactor", 0.5
        )
    )
{
    // Force p oldTime to be stored
    p_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool coupledLinGeomPressureDisplacementSolid::evolve()
{
#ifdef FOAMEXTEND
    Info<< "Evolving solid solver" << endl;

    // Disable default writing of linear solver residuals
    blockLduMatrix::debug = 0;

    // Mesh update loop
    do
    {
        int iCorr = 0;
        BlockSolverPerformance< VectorN<double, 4> > solverPerfDp;

        Info<< "Solving the momentum equation for D and p" << endl;

        // Loop around displacement and pressure equations
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Initialize the Dp block system (matrix, source and reference to
            // Dp)
            fvBlockMatrix<vector4> DpEqn(Dp_);

            // Momentum equation in terms of displacement
            // Note: we explicitly pressure remove it so that we add implicitly
            // add it below as "pInD"
            fvVectorMatrix DEqn
            (
                fvm::d2dt2(rho(), D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(sigma(), "div(sigma)")
              + gradp_
              + rho()*g()
              + stabilisation().stabilisation(D(), gradD(), impK_)
            );

            // Store reciprocal of diagonal
            const surfaceScalarField rAUf = fvc::interpolate(1.0/DEqn.A());

            // Under-relaxation the linear system
            DEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DEqn);

            // Insert displacement equation into block system
            DpEqn.insertEquation(0, DEqn);

            // Store fields for under-relaxation and residual
            // calculation
            p_.storePrevIter();

            // Pressure equation
            // Note: div(D) is add implicitly below as "Dinp"
            fvScalarMatrix pEqn
            (
                fvm::Sp(rK_, p_)
              - fvm::laplacian(rAUf, p_, "laplacian(Dp,p)")
              + fvc::laplacian(rAUf, p_, "laplacian(Dp,p)")
             ==
                //- fvc::div(D())
              + pressureRhieChowScaleFac_
               *(
                    fvc::laplacian(rAUf, p_, "laplacian(Dp,p)")
                  - fvc::div(rAUf*mesh().Sf() & fvc::interpolate(gradp_))
                )
            );

            // Under-relaxation the linear system
            pEqn.relax();

            // Insert pressure equation into block system
            DpEqn.insertEquation(3, pEqn);

            // Insert coupling terms into block system
            {
                // Calculate grad p coupling matrix
                BlockLduSystem<vector, vector> pInD(fvm::grad(p_));

                // Calculate div D coupling
                BlockLduSystem<vector, scalar> DInp(fvm::UDiv(D()));

                DpEqn.insertBlockCoupling(0, 3, pInD, true);
                DpEqn.insertBlockCoupling(3, 0, DInp, false);

                // Not used. We could use this instead of explicitly leaving
                // out/removing these terms
                // Update source coupling: coupling terms eliminated from source
                //blockM.updateSourceCoupling();
            }

            // Solve the block matrix
            solverPerfDp = DpEqn.solve();

            // Retrieve solution
            DpEqn.retrieveSolution(0, D().internalField());
            DpEqn.retrieveSolution(3, p_.internalField());

            D().correctBoundaryConditions();
            p_.correctBoundaryConditions();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);
            p_.relax();

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Update the gradient of pressure
            gradp_ = fvc::grad(p_);

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());
        }
        while
        (
            !converged
            (
                iCorr,
                cmptMax(solverPerfDp.initialResidual()),
                solverPerfDp.nIterations(),
                D()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    blockLduMatrix::debug = 1;
#else
    notImplemented("Not implemented (yet) for this version of OpenFOAM/FOAM");
#endif

    return true;
}


tmp<vectorField> coupledLinGeomPressureDisplacementSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - n*pressure)
                  - (n & (pSigma - impK*pGradD))
                )*rImpK
            )
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
