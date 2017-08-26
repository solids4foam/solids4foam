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

#include "coupledHybridLinGeomSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"
#include "fvBlockMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledHybridLinGeomSolid, 0);
addToRunTimeSelectionTable(physicsModel, coupledHybridLinGeomSolid, solid);
addToRunTimeSelectionTable(solidModel, coupledHybridLinGeomSolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledHybridLinGeomSolid::coupledHybridLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    pD_
    (
        IOobject
        (
            "pD",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector4("zero", dimless, vector4::zero)
    ),
    gradp_
    (
        IOobject
        (
            "grad(" + p_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimPressure/dimLength, vector::zero)
    ),
    rK_(1.0/mechanical().K()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool coupledHybridLinGeomSolid::evolve()
{
    notImplemented("bool coupledHybridLinGeomSolid::evolve()");

//     Info<< "Evolving solid solver" << endl;

//     int iCorr = 0;
//     lduSolverPerformance solverPerfD;
//     lduSolverPerformance solverPerfp;
//     blockLduMatrix::debug = 0;
//     BlockLduMatrix<vector4>::debug = 0;

//     Info<< "Solving the momentum equation for D and p using a coupled "
//         << "approach" << endl;

//     // Momentum equation loop
//     do
//     {
//         // Store fields for under-relaxation and residual calculation
//         D_.storePrevIter();
//         p_.storePrevIter();

//         // Initialize the pD block system
//         fvBlockMatrix<vector4> pDEqn(pD_);

//         // Assemble and insert momentum equation with pressure
//         // term removed i.e. we add it back on
//         tmp<fvVectorMatrix> DEqn
//         (
//             rho_*fvm::d2dt2(D_)
//          == fvm::laplacian(impKf_, D_, "laplacian(DD,D)")
//           - fvc::laplacian(impKf_, D_, "laplacian(DD,D)")
//           + fvc::div(sigma_, "div(sigma)")
//           + gradp_
//           + rho_*g_
//           + mechanical().RhieChowCorrection(D_, gradD_)
//         );

//         // Under-relax the equation
//         DEqn().relax();

//         // Inser the equation into the block system
//         pDEqn.insertEquation(0, DEqn());

//         // Assemble and insert pressure equation

//         // Pressure parts of the continuity equation
//         const surfaceScalarField rUAf
//         (
//             "rUAf",
//             fvc::interpolate(1.0/(DEqn()).A())
//         );

//         // Clear the DEqn to reduce memory overhead
//         DEqn.clear();

//         // Pressure equation where we allow for explicit skewness correction
//         // and Rhie-Chow corrections are also employed
//         fvScalarMatrix pEqn
//         (
//             fvm::Sp(rK_, p_)
//           - fvm::laplacian(rUAf, p_)
//           + fvc::div(rUAf*(fvc::interpolate(fvc::grad(p_)) & mesh().Sf()))
//           + fvc::div(D_, "skewCorrectedDiv(D)")
//           - fvc::div(D_)
//         );

//         // Under-relax the equation
//         pEqn.relax();

//         // Insert the equation into the block system
//         pDEqn.insertEquation(3, pEqn);

//         // Assemble and insert coupling terms

//         // Calculate grad p coupling matrix
//         BlockLduSystem<vector, vector> pInD(fvm::grad(p_));

//         // Calculate div(D) coupling
//         BlockLduSystem<vector, scalar> DInp(fvm::UDiv(D_));

//         // Last argument in insertBlockCoupling says if the column direction
//         // should be incremented. This is needed for arbitrary positioning
//         // of D and p in the system
//         pDEqn.insertBlockCoupling(0, 3, pInD, true);
//         pDEqn.insertBlockCoupling(3, 0, DInp, false);

//         // Solve the block matrix
//         BlockSolverPerformance<vector4> blockSolverPerf = pDEqn.solve();

//         // Store the residuals and iterations
//         solverPerfD.initialResidual() =
//             cmptMax(blockSolverPerf.initialResidual());
//         solverPerfD.nIterations() = blockSolverPerf.nIterations();

//         // Retrieve solution
//         pDEqn.retrieveSolution(0, D_.internalField());
//         pDEqn.retrieveSolution(3, p_.internalField());

//         D_.correctBoundaryConditions();
//         p_.correctBoundaryConditions();

//         // Under-relax the field
//         D_.relax();

//         // Update gradient of displacement
//         mechanical().grad(D_, gradD_);

//         // Explicitly relax pressure for momentum corrector
//         p_.relax();

//         // Update gradient of pressure
//         gradp_ = fvc::grad(p_);

//         // Calculate the stress using run-time selectable mechanical law
//         mechanical().correct(sigma_);
//     }
//     while
//     (
//         !converged(iCorr, solverPerfD.initialResidual(), D())
//      && ++iCorr < nCorr_
//     );

//     // Interpolate cell displacements to vertices
//     mechanical().interpolate(D_, pointD_);

//     // Velocity
//     U_ = fvc::ddt(D_);

    return true;
}


tmp<vectorField> coupledHybridLinGeomSolid::tractionBoundarySnGrad
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
    const vectorField pN = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - pN*pressure)
              - (pN & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
