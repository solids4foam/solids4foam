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


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool coupledHybridLinGeomSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfD,
    const lduSolverPerformance& solverPerfp
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

    // Calculate pressure residual
    const scalar residualp =
        gMax
        (
            mag(p_.internalField() - p_.prevIter().internalField())
           /max
            (
                gMax(mag(p_.internalField() - p_.oldTime().internalField())),
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
         && solverPerfp.initialResidual() < solutionTolp_
         && residualp < solutionTolp_
        )
        {
            Info<< "    All residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            residualD < alternativeTol_
         && residualp < alternativeTolp_
        )
        {
            Info<< "    The relative residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            solverPerfD.initialResidual() < alternativeTol_
         && solverPerfp.initialResidual() < alternativeTolp_
        )
        {
            Info<< "    The solver residuals have converged" << endl;
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
        Info<< "    Corr, resPD, relResD, relResP, matRes, itersPD" << endl;
    }
    else if (iCorr % infoFrequency_ == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfD.initialResidual()
            << ", " << residualD
            << ", " << residualp
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledHybridLinGeomSolid::coupledHybridLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    D_
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
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
    U_
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength/dimTime, vector::zero)
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
    pMesh_(mesh()),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime.timeName(),
            mesh(),
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
            "sigma",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
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
    rho_(mechanical().rho()),
    rK_(1.0/mechanical().K()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    solutionTol_
    (
        solidProperties().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
    ),
    solutionTolp_
    (
        solidProperties().lookupOrDefault<scalar>
        (
            "solutionToleranceP", solutionTol_
        )
    ),
    alternativeTol_
    (
        solidProperties().lookupOrDefault<scalar>("alternativeTolerance", 1e-07)
    ),
    alternativeTolp_
    (
        solidProperties().lookupOrDefault<scalar>
        (
            "alternativeToleranceP", alternativeTol_
        )
    ),
    materialTol_
    (
        solidProperties().lookupOrDefault<scalar>("materialTolerance", 1e-05)
    ),
    infoFrequency_
    (
        solidProperties().lookupOrDefault<int>("infoFrequency", 100)
    ),
    nCorr_(solidProperties().lookupOrDefault<int>("nCorrectors", 10000)),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    maxIterReached_(0)
{
    D_.oldTime().oldTime();
    pointD_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> coupledHybridLinGeomSolid::faceZonePointDisplacementIncrement
(
    const label zoneID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(),
            vector::zero
        )
    );
    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& pointDI = pointD_.internalField();
    const vectorField& oldPointDI = pointD_.oldTime().internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    pointDI[procPoint] - oldPointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] =
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        tPointDisplacement() =
            vectorField
            (
                pointDI - oldPointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


tmp<tensorField> coupledHybridLinGeomSolid::faceZoneSurfaceGradientOfVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<tensorField> tVelocityGradient
    (
        new tensorField
        (
            mesh().faceZones()[zoneID]().size(),
            tensor::zero
        )
    );
    tensorField& velocityGradient = tVelocityGradient();

    vectorField pPointU =
        mechanical().volToPoint().interpolate
        (
            mesh().boundaryMesh()[patchID],
            U_
        );

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints += pointD_.boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    tensorField patchGradU = fvc::fGrad(patch, pPointU);

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart =
            mesh().boundaryMesh()[patchID].start();

        forAll(patchGradU, i)
        {
            velocityGradient
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] =
                patchGradU[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(velocityGradient, sumOp<tensorField>());
    }
    else
    {
        velocityGradient = patchGradU;
    }

    return tVelocityGradient;
}


tmp<vectorField> coupledHybridLinGeomSolid::currentFaceZonePoints
(
    const label zoneID
) const
{
    vectorField pointDisplacement
    (
        mesh().faceZones()[zoneID]().localPoints().size(),
        vector::zero
    );

    const vectorField& pointDI = pointD_.internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone
        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    pointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] =
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        pointDisplacement =
            vectorField
            (
                pointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    tmp<vectorField> tCurrentPoints
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints()
          + pointDisplacement
        )
    );

    return tCurrentPoints;
}


tmp<vectorField> coupledHybridLinGeomSolid::faceZoneNormal
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tNormals
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& normals = tNormals();

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints += pointD_.boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    vectorField patchNormals(patch.size(), vector::zero);

    forAll(patchNormals, faceI)
    {
        patchNormals[faceI] =
            localFaces[faceI].normal(localPoints);
    }

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart =
            mesh().boundaryMesh()[patchID].start();

        forAll(patchNormals, i)
        {
            normals
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] =
                patchNormals[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(normals, sumOp<vectorField>());
    }
    else
    {
        normals = patchNormals;
    }

    return tNormals;
}


void coupledHybridLinGeomSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != solidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void coupledHybridLinGeomSolid::setTraction(...)")
            << "Bounary condition on " << D_.name()
            <<  " is "
            << D_.boundaryField()[patchID].type()
            << "for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    solidTractionFvPatchVectorField& patchU =
        refCast<solidTractionFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.traction() = traction;
}


void coupledHybridLinGeomSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != solidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void coupledHybridLinGeomSolid::setTraction(...)")
            << "Bounary condition on " << D_.name()
            <<  " is "
            << D_.boundaryField()[patchID].type()
            << "for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    solidTractionFvPatchVectorField& patchU =
        refCast<solidTractionFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.pressure() = pressure;
}


bool coupledHybridLinGeomSolid::evolve()
{
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
//     while (!converged(iCorr, solverPerfD, solverPerfp) && ++iCorr < nCorr_);

//     // Interpolate cell displacements to vertices
//     mechanical().interpolate(D_, pointD_);

//     // Velocity
//     U_ = fvc::ddt(D_);

    return true;
}


// void coupledHybridLinGeomSolid::predict()
// {
//     Info << "Predicting solid model" << endl;

//     D_ = D_ + U_*runTime().deltaT();

//     if (interface().valid())
//     {
//         interface()->updateDisplacement(pointD_);
//         interface()->updateDisplacementGradient(gradD_, gradDf_);
//     }
//     else
//     {
//         volToPoint_.interpolate(D_, pointD_);
//         gradD_ = fvc::grad(D_, pointD_);
//         gradDf_ = fvc::fGrad(D_, pointD_);
//     }

//     D_.correctBoundaryConditions();
// }


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
    const tensorField& gradD = gradD_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigma_.boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (sigma - impK*gradD))
            )*rImpK
        )
    );
}


void coupledHybridLinGeomSolid::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void coupledHybridLinGeomSolid::writeFields(const Time& runTime)
{
    // Calculate strain
    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm(gradD_)
    );

    // Calculaute equivalent strain
    volScalarField epsilonEq
    (
        IOobject
        (
            "epsilonEq",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    );

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
        << endl;

    // Calculaute equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma_)))
    );

    Info<< "Max sigmaEq = " << gMax(sigmaEq) << endl;

    solidModel::writeFields(runTime);
}


void coupledHybridLinGeomSolid::end()
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

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
