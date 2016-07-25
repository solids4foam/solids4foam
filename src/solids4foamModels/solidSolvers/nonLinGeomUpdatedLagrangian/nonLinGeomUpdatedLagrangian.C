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

#include "nonLinGeomUpdatedLagrangian.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
//#include "solidTractionFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomUpdatedLagrangian, 0);
addToRunTimeSelectionTable
(
    solidSolver, nonLinGeomUpdatedLagrangian, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool nonLinGeomUpdatedLagrangian::converged
(
    const int iCorr,
    const lduMatrix::solverPerformance& solverPerfDD
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate displacement residual
    const scalar residualDD =
        gMax
        (
            mag(DD_.internalField() - DD_.prevIter().internalField())
           /max
            (
                gMax(mag(DD_.internalField())),
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
            solverPerfDD.initialResidual() < solutionTol_
         && residualDD < solutionTol_
        )
        {
            Info<< "    Both residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            residualDD < alternativeTol_
        )
        {
            Info<< "    The relative residual has converged" << endl;
            converged = true;
        }
        else if
        (
            solverPerfDD.initialResidual() < alternativeTol_
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
            << ", " << solverPerfDD.initialResidual()
            << ", " << residualDD
            << ", " << materialResidual
            << ", " << solverPerfDD.nIterations() << endl;

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


// void nonLinGeomUpdatedLagrangian::checkJacobian(const volScalarField& J)
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


void nonLinGeomUpdatedLagrangian::moveMesh(const pointField& oldPoints)
{
    //- Move mesh by interpolating displacement field to vertices
    // TO be checked: sync boundary and global points across procs to make sure
    // numiercal error does not build up and when end up with the error
    // "face area does not match neighbour..."
    // We could sync points as a pointVectorField just as we sync pointDD

    // Interppolate cell displacements to vertices
    volToPoint_.interpolate(DD_, pointDD_);

    // Ensure continuous displacement across processor boundary
    // Something strange is happening here
    pointDD_.correctBoundaryConditions();

    vectorField& pointDDI = pointDD_.internalField();

    vectorField newPoints = oldPoints;

    // Correct symmetryPlane points

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());

            const vector i(1, 0, 0);
            const vector j(0, 1, 0);
            const vector k(0, 0, 1);

            if (mag(avgN & i) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].x() = 0;
                }
            }
            else if (mag(avgN & j) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].y() = 0;
                }
            }
            else if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
        else if (isA<emptyPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());
            const vector k(0, 0, 1);

            if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
    }

    // Note: allPoints will have more points than pointDD if there are
    // globalFaceZones
    forAll (pointDDI, pointI)
    {
        newPoints[pointI] += pointDDI[pointI];
    }

    // Move unused globalFaceZone points
    updateGlobalFaceZoneNewPoints(pointDDI, newPoints);

    twoDPointCorrector twoDCorrector(mesh());
    twoDCorrector.correctPoints(newPoints);
    twoDCorrector.correctPoints(pointDD_.internalField());
    mesh().movePoints(newPoints);
    mesh().V00();
    mesh().moving(false);
    mesh().changing(false);

    // meshPhi does not need to be written
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomUpdatedLagrangian::nonLinGeomUpdatedLagrangian(fvMesh& mesh)
:
    solidSolver(typeName, mesh),
    DD_
    (
        IOobject
        (
            "DD",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    D_
    (
        IOobject
        (
            "D",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pMesh_(mesh),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
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
    volToPoint_(mesh),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
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
    relF_
    (
        IOobject
        (
            "relF",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + gradDD_.T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        hinv(relF_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    rho_(mechanical().rho()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    DDEqnRelaxFactor_
    (
        mesh.solutionDict().relax("DDEqn")
      ? mesh.solutionDict().relaxationFactor("DDEqn")
      : 1.0
    ),
    solutionTol_(lookupOrDefault<scalar>("solutionTolerance", 1e-06)),
    alternativeTol_(lookupOrDefault<scalar>("alternativeTolerance", 1e-07)),
    materialTol_(lookupOrDefault<scalar>("materialTolerance", 1e-05)),
    infoFrequency_(lookupOrDefault<int>("infoFrequency", 100)),
    nCorr_(lookupOrDefault<int>("nCorrectors", 1000)),
    maxIterReached_(0),
    stabilisePressure_(lookupOrDefault<Switch>("stabilisePressure", false))
{
    DD_.oldTime().oldTime();
    D_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// void nonLinGeomUpdatedLagrangian::setTraction
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
//         FatalErrorIn("void nonLinGeomUpdatedLagrangian::setTraction(...)")
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

// void nonLinGeomUpdatedLagrangian::setPressure
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
//         FatalErrorIn("void nonLinGeomUpdatedLagrangian::setTraction(...)")
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

// void nonLinGeomUpdatedLagrangian::setTraction
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

// void nonLinGeomUpdatedLagrangian::setPressure
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


bool nonLinGeomUpdatedLagrangian::evolve()
{
    Info << "Evolving solid solver" << endl;

    int iCorr = 0;
    lduMatrix::solverPerformance solverPerfDD;
    lduMatrix::debug = 0;

    // Store old points for moving the mesh
    const vectorField oldPoints = mesh().allPoints();

    Info<< "Solving the momentum equation for DD" << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD_.storePrevIter();

        // Momentum equation incremental updated Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho_, DD_)
          + fvc::d2dt2(rho_.oldTime(), D_.oldTime())
         == fvm::laplacian(impKf_, DD_, "laplacian(DDD,DD)")
          + fvc::div
            (
                (relJ_*sigma_ & relFinv_.T()) - impK_*gradDD_,
                "div(sigma)"
            )
        );

        // Under-relax the linear system
        DDEqn.relax(DDEqnRelaxFactor_);

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field
        DD_.relax();

        // Update gradient of displacement increment
        gradDD_ = fvc::grad(DD_);

        // Relative deformation gradient
        relF_ = I + gradDD_.T();

        // Inverse relative deformation gradient
        relFinv_ = hinv(relF_);

        // Total deformation gradient
        F_ = relF_ & F_.oldTime();

        // Relative Jacobian (Jacobian of relative deformation gradient)
        if (stabilisePressure_)
        {
            // Reconstruct face gradients to relJ: this will avoid
            // pressure oscillations
            // This is just det(I + gradDD.T()) where the gradDD has
            // been averaged from the face gradDD field
            relJ_ =
                det
                (
                    I
                    + fvc::reconstruct
                    (
                        mesh().magSf()*fvc::snGrad(DD_)
                    )().T()
                );
        }
        else
        {
            relJ_ = det(relF_);
        }

        // Bound relative Jacobian to improve robustness
        boundMinMax
        (
            relJ_,
            dimensionedScalar("smallJ", dimless, 0.1),
            dimensionedScalar("smallJ", dimless, 2.0)
        );

        // Jacobian of deformation gradient
        J_ = relJ_*J_.oldTime();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma_);
    }
    while (!converged(iCorr, solverPerfDD) && ++iCorr < nCorr_);

    // PC: rename this function or maybe even remove it
    // Update yield stress and plasticity total field e.g. epsilonP
    // Or updateTotalFields: actually, this should be called inside
    // updateTotalFields() that gets called in solidFoam
    mechanical().updateYieldStress();

    // Update gradient of total displacement
    gradD_ = fvc::grad(D_.oldTime() + DD_);

    // Total displacement
    D_ = D_.oldTime() + DD_;

    // Density
    rho_ = rho_.oldTime()/relJ_;

    // Update to current configuration after TEqn
    // Note: that the energy equation moves the mesh to the mid-step of
    // the time-step when using convective acceleration
    moveMesh(oldPoints);

    // Total displacement at points
    pointD_ = pointD_.oldTime() + pointDD_;

    return true;
}


tmp<vectorField> nonLinGeomUpdatedLagrangian::tractionBoundarySnGrad
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
    const tensorField& gradDD = gradDD_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigma_.boundaryField()[patchID];

    // Patch relative deformation gradient inverse
    const tensorField& relFinv = relFinv_.boundaryField()[patchID];

    // Patch relative Jacobian
    const scalarField& relJ = relJ_.boundaryField()[patchID];

    // Patch unit normals (updated configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    const vectorField nCurrent = relJ*relFinv.T() & n;

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (nCurrent & sigma)
              + (n & (impK*gradDD))
            )*rImpK
        )
    );
}


void nonLinGeomUpdatedLagrangian::writeFields(const Time& runTime)
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


void nonLinGeomUpdatedLagrangian::end()
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
