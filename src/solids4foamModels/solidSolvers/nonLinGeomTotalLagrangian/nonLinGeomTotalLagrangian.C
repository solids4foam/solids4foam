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

#include "nonLinGeomTotalLagrangian.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
//#include "solidTractionFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagrangian, 0);
addToRunTimeSelectionTable
(
    solidSolver, nonLinGeomTotalLagrangian, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool nonLinGeomTotalLagrangian::converged
(
    const int iCorr,
    const lduMatrix::solverPerformance& solverPerfD
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate displacement residual
    const scalar residualD =
        gMax
        (
            mag(D_.internalField() - D_.prevIter().internalField())
           /max
            (
                gMax(mag(D_.internalField() - D_.oldTime().internalField())),
                SMALL
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol_)
    {
        if
        (
            solverPerfD.initialResidual() < solutionTol_
         && residualD < solutionTol_
        )
        {
            Info<< "    Both residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            residualD < alternativeTol_
        )
        {
            Info<< "    The relative residual has converged" << endl;
            converged = true;
        }
        else if
        (
            solverPerfD.initialResidual() < alternativeTol_
        )
        {
            Info<< "    The solver residual has converged" << endl;
            converged = true;
        }
        else
        {
            converged = false;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res, relRes, matRes, iters" << endl;
    }
    else if (iCorr % infoFrequency_ == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfD.initialResidual()
            << ", " << residualD
            << ", " << materialResidual
            << ", " << solverPerfD.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr_ - 1)
    {
        maxIterReached_++;
        Warning
            << "Max iterations reached within momentum loop" << endl;
    }

    return converged;
}


// void nonLinGeomTotalLagrangian::checkJacobian(const volScalarField& J)
// {
//     const scalarField& JI = J.internalField();

//     if (gMax(JI) < 0.01)
//     {
//         forAll(JI, cellI)
//         {
//             if (JI[cellI] < SMALL)
//             {
//                 Pout<< "Cell " << cellI
//                     << " with centre " << mesh.C()[cellI]
//                     << " has a become inverted!" << endl;
//             }
//         }

//         FatalErrorIn(type() + "::evolve()")
//             << "Cells have become inverted! see details above."
//             << abort(FatalError);
//     }
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagrangian::nonLinGeomTotalLagrangian(fvMesh& mesh)
:
    solidSolver(typeName, mesh),
    D_
    (
        IOobject
        (
            "D",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    sigma_
    (
        IOobject
        (
            "sigmaCauchy",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        hinv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    rho_(mechanical().rho()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    DEqnRelaxFactor_
    (
        mesh.solutionDict().relax("DEqn")
      ? mesh.solutionDict().relaxationFactor("DEqn")
      : 1.0
    ),
    solutionTol_(lookupOrDefault<scalar>("solutionTolerance", 1e-06)),
    alternativeTol_(lookupOrDefault<scalar>("alternativeTolerance", 1e-07)),
    materialTol_(lookupOrDefault<scalar>("materialTolerance", 1e-05)),
    infoFrequency_(lookupOrDefault<int>("infoFrequency", 100)),
    nCorr_(lookupOrDefault<int>("nCorrectors", 1000)),
    maxIterReached_(0)
{
    D_.oldTime().oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// void nonLinGeomTotalLagrangian::setTraction
// (
//     const label patchID,
//     const vectorField& traction
// )
// {
//     if
//     (
//         D_.boundaryField()[patchID].type()
//      != solidTractionFvPatchVectorField::typeName
//     )
//     {
//         FatalErrorIn("void nonLinGeomTotalLagrangian::setTraction(...)")
//             << "Bounary condition on " << D_.name()
//             <<  " is "
//             << D_.boundaryField()[patchID].type()
//             << "for patch" << mesh().boundary()[patchID].name()
//             << ", instead "
//             << solidTractionFvPatchVectorField::typeName
//             << abort(FatalError);
//     }

//     solidTractionFvPatchVectorField& patchU =
//         refCast<solidTractionFvPatchVectorField>
//         (
//             D_.boundaryField()[patchID]
//         );

//     patchU.traction() = traction;
// }

// void nonLinGeomTotalLagrangian::setPressure
// (
//     const label patchID,
//     const scalarField& pressure
// )
// {
//     if
//     (
//         D_.boundaryField()[patchID].type()
//      != solidTractionFvPatchVectorField::typeName
//     )
//     {
//         FatalErrorIn("void nonLinGeomTotalLagrangian::setTraction(...)")
//             << "Bounary condition on " << D_.name()
//             <<  " is "
//             << D_.boundaryField()[patchID].type()
//             << "for patch" << mesh().boundary()[patchID].name()
//             << ", instead "
//             << solidTractionFvPatchVectorField::typeName
//             << abort(FatalError);
//     }

//     solidTractionFvPatchVectorField& patchU =
//         refCast<solidTractionFvPatchVectorField>
//         (
//             D_.boundaryField()[patchID]
//         );

//     patchU.pressure() = pressure;
// }

// void nonLinGeomTotalLagrangian::setTraction
// (
//     const label patchID,
//     const label zoneID,
//     const vectorField& faceZoneTraction
// )
// {
//   vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);

//     const label patchStart =
//         mesh().boundaryMesh()[patchID].start();

//     forAll(patchTraction, i)
//     {
//         patchTraction[i] =
//             faceZoneTraction
//             [
//                 mesh().faceZones()[zoneID].whichFace(patchStart + i)
//             ];
//     }

//     setTraction(patchID, patchTraction);
// }

// void nonLinGeomTotalLagrangian::setPressure
// (
//     const label patchID,
//     const label zoneID,
//     const scalarField& faceZonePressure
// )
// {
//     scalarField patchPressure(mesh().boundary()[patchID].size(), 0.0);

//     const label patchStart =
//         mesh().boundaryMesh()[patchID].start();

//     forAll(patchPressure, i)
//     {
//         patchPressure[i] =
//             faceZonePressure
//             [
//                 mesh().faceZones()[zoneID].whichFace(patchStart + i)
//             ];
//     }

//     setPressure(patchID, patchPressure);
// }


bool nonLinGeomTotalLagrangian::evolve()
{
    Info << "Evolving solid solver" << endl;

    int iCorr = 0;
    lduMatrix::solverPerformance solverPerfD;
    lduMatrix::debug = 0;

    Info<< "Solving the momentum equation for D" << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D_.storePrevIter();

        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            rho_*fvm::d2dt2(D_)
         == fvm::laplacian(impKf_, D_, "laplacian(DD,D)")
          + fvc::div
            (
                (J_*sigma_ & Finv_.T()) - impK_*gradD_,
                "div(sigma)"
            )
        );

        // Under-relax the linear system
        DEqn.relax(DEqnRelaxFactor_);

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the D field
        D_.relax();

        // Update gradient of displacement increment
        gradD_ = fvc::grad(D_);

        // Total deformation gradient
        F_ = I + gradD_.T();

        // Inverse of the deformation gradient
        Finv_ = hinv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma_);
    }
    while (!converged(iCorr, solverPerfD) && ++iCorr < nCorr_);

    // PC: rename this function or maybe even remove it
    // Update yield stress and plasticity total field e.g. epsilonP
    // Or updateTotalFields: actually, this should be called inside
    // updateTotalFields() that gets called in solidFoam
    mechanical().updateYieldStress();

    return true;
}


tmp<vectorField> nonLinGeomTotalLagrangian::tractionBoundarySnGrad
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
    const tensorField& gradD = gradD_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigma_.boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch total Jacobian
    const scalarField& J = J_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    const vectorField nCurrent = J*Finv.T() & n;

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (nCurrent & sigma)
              + (n & (impK*gradD))
            )*rImpK
        )
    );
}


void nonLinGeomTotalLagrangian::writeFields(const Time& runTime)
{
    // Update equivalent strain
    // volScalarField epsilonEq
    // (
    //     IOobject
    //     (
    //         "epsilonEq",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     sqrt((2.0/3.0)*magSqr(dev(epsilon_)))
    // );

    // Info<< "Max epsilonEq = " << max(epsilonEq).value()
    //     << endl;

    // Update equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaCauchyEq",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma_)))
    );

    Info<< "Max sigmaCauchyEq = " << gMax(sigmaEq) << endl;

    solidSolver::writeFields(runTime);
}


void nonLinGeomTotalLagrangian::end()
{
    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum momentum correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The momentum equation converged in all time-steps"
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
