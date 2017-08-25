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

#include "unsNonLinGeomUpdatedLagSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsNonLinGeomUpdatedLagSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, unsNonLinGeomUpdatedLagSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, unsNonLinGeomUpdatedLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool unsNonLinGeomUpdatedLagSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfDD
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


void unsNonLinGeomUpdatedLagSolid::moveMesh(const pointField& oldPoints)
{
    Info<< "Moving the mesh to the deformed configuration" << nl << endl;

    //- Move mesh by interpolating displacement field to vertices
    // To be checked: sync boundary and global points across procs to make sure
    // numiercal error does not build up and when end up with the error
    // "face area does not match neighbour..."
    // We could sync points as a pointVectorField just as we sync pointDD

    // Interpolate cell displacements to vertices
    mechanical().volToPoint().interpolate(DD_, pointDD_);

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

unsNonLinGeomUpdatedLagSolid::unsNonLinGeomUpdatedLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    DD_
    (
        IOobject
        (
            "DD",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    D_
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
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
    pMesh_(mesh()),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            runTime.timeName(),
            mesh(),
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
            "sigmaCauchy",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    sigmaf_
    (
        IOobject
        (
            "sigmaCauchyf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDDf_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
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
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
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
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(Ff_)
    ),
    relF_
    (
        IOobject
        (
            "relF",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + gradDD_.T()
    ),
    relFf_
    (
        IOobject
        (
            "relFf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + gradDDf_.T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relF_)
    ),
    relFinvf_
    (
        IOobject
        (
            "relFinvf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relFf_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    relJf_
    (
        IOobject
        (
            "relJf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relFf_)
    ),
    rho_(mechanical().rho()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    solutionTol_
    (
        solidProperties().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
    ),
    alternativeTol_
    (
        solidProperties().lookupOrDefault<scalar>("alternativeTolerance", 1e-07)
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
    maxIterReached_(0),
    stabilisePressure_(lookupOrDefault<Switch>("stabilisePressure", false))
{
    DD_.oldTime().oldTime();
    D_.oldTime();
    pointDD_.oldTime();
    pointD_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField>
unsNonLinGeomUpdatedLagSolid::faceZonePointDisplacementIncrement
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

    const vectorField& pointDDI = pointDD_.internalField();

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

            if (zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    pointDDI[procPoint];

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
                pointDDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


tmp<tensorField> unsNonLinGeomUpdatedLagSolid::faceZoneSurfaceGradientOfVelocity
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
        mechanical().volToPoint().interpolate(mesh().boundaryMesh()[patchID], U_);

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
            ] = patchGradU[i];
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


tmp<vectorField> unsNonLinGeomUpdatedLagSolid::currentFaceZonePoints
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

            if (zoneMeshPoints[localPoint] < mesh().nPoints())
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


tmp<vectorField> unsNonLinGeomUpdatedLagSolid::faceZoneNormal
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
            ] = patchNormals[i];
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


void unsNonLinGeomUpdatedLagSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        DD_.boundaryField()[patchID].type()
     != solidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsNonLinGeomUpdatedLagSolid::setTraction(...)")
            << "Boundary condition on " << DD_.name()
            <<  " is "
            << DD_.boundaryField()[patchID].type()
            << " for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    solidTractionFvPatchVectorField& patchDD =
        refCast<solidTractionFvPatchVectorField>
        (
            DD_.boundaryField()[patchID]
        );

    patchDD.traction() = traction;
}


void unsNonLinGeomUpdatedLagSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        DD_.boundaryField()[patchID].type()
     != solidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsNonLinGeomUpdatedLagSolid::setTraction(...)")
            << "Boundary condition on " << DD_.name()
            <<  " is "
            << DD_.boundaryField()[patchID].type()
            << " for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    solidTractionFvPatchVectorField& patchDD =
        refCast<solidTractionFvPatchVectorField>
        (
            DD_.boundaryField()[patchID]
        );

    patchDD.pressure() = pressure;
}


bool unsNonLinGeomUpdatedLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfDD;
    blockLduMatrix::debug = 0;

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
          - fvc::laplacian(impKf_, DD_, "laplacian(DDD,DD)")
          + fvc::div((relJf_*relFinvf_.T() & mesh().Sf()) & sigmaf_)
          + rho_*g_
        );

        // Under-relax the linear system
        DDEqn.relax();

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field
        DD_.relax();

        // Interpolate DD to pointDD
        mechanical().interpolate(DD_, pointDD_, false);

        // Update gradient of displacement
        mechanical().grad(DD_, pointDD_, gradDD_);
        mechanical().grad(DD_, pointDD_, gradDDf_);

        // Relative deformation gradient
        relFf_ = I + gradDDf_.T();

        // Inverse relative deformation gradient
        relFinvf_ = inv(relFf_);

        // Total deformation gradient
        Ff_ = relFf_ & Ff_.oldTime();

        // Relative Jacobian
        relJf_ = det(relFf_);

        // Jacobian of deformation gradient
        Jf_ = relJf_*Jf_.oldTime();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);
    }
    while (!converged(iCorr, solverPerfDD) && ++iCorr < nCorr_);

    // Relative deformation gradient
    relF_ = I + gradDD_.T();

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma_);

    // Update gradient of total displacement
    gradD_ = fvc::grad(D_.oldTime() + DD_);

    // Total displacement
    D_ = D_.oldTime() + DD_;

    // Total displacement at points
    pointD_ = pointD_.oldTime() + pointDD_;

    // Velocity
    U_ = fvc::ddt(D_);

    return true;
}


tmp<vectorField> unsNonLinGeomUpdatedLagSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impKf_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& gradDD = gradDDf_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigmaf_.boundaryField()[patchID];

    // Patch relative deformation gradient inverse
    const tensorField& relFinv = relFinvf_.boundaryField()[patchID];

    // Patch relative Jacobian
    const scalarField& relJ = relJf_.boundaryField()[patchID];

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


void unsNonLinGeomUpdatedLagSolid::updateTotalFields()
{
    // Density
    rho_ = rho_.oldTime()/relJ_;

    // Move the mesh to the deformed configuration
    const vectorField oldPoints = mesh().allPoints();
    moveMesh(oldPoints);

    mechanical().updateTotalFields();
}


void unsNonLinGeomUpdatedLagSolid::writeFields(const Time& runTime)
{
    // Calculate equivalent (von Mises) stress
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

    solidModel::writeFields(runTime);
}


void unsNonLinGeomUpdatedLagSolid::end()
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
