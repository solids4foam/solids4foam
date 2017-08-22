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

#include "coupledUnsLinGeomLinearElasticSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "blockSolidTractionFvPatchVectorField.H"
#include "fvcGradf.H"
#include "BlockFvmDivSigma.H"
#include "linearElastic.H"

#include "blockFixedDisplacementZeroShearFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledUnsLinGeomLinearElasticSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, coupledUnsLinGeomLinearElasticSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledUnsLinGeomLinearElasticSolid::coupledUnsLinGeomLinearElasticSolid
(
    dynamicFvMesh& mesh
)
:
    solidModel(typeName, mesh),
    extendedMesh_(mesh),
    solutionVec_
    (
        IOobject
        (
            "solutionVec",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField(extendedMesh_.nVariables(), vector::zero)
    ),
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
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pMesh_(mesh),
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    volToPoint_(mesh),
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
    rho_(mechanical().rho()),
    muf_
    (
        IOobject
        (
            "interpolate(mu)",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    lambdaf_
    (
        IOobject
        (
            "interpolate(lambda)",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    g_
    (
        IOobject
        (
            "g",
            runTime().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    DEqnRelaxFactor_
    (
        mesh.solutionDict().relax("DEqn")
      ? mesh.solutionDict().relaxationFactor("DEqn")
      : 1.0
    )
{
    D_.oldTime().oldTime();
    pointD_.oldTime();

    // We will directly read the linearElastic mechanicalLaw
    const PtrList<mechanicalLaw>& mechLaws = mechanical();
    if (mechLaws.size() != 1)
    {
        FatalErrorIn
        (
            "coupledUnsLinGeomLinearElasticSolid::"
            "coupledUnsLinGeomLinearElasticSolid"
        )   << type() << " can currently only be used with a single material"
            << "\nConsider using one of the other solidModels."
            << abort(FatalError);
    }
    else if (!isA<linearElastic>(mechLaws[0]))
    {
        FatalErrorIn
        (
            "coupledUnsLinGeomLinearElasticSolid::"
            "coupledUnsLinGeomLinearElasticSolid"
        )   << type() << " can only be used with the linearElastic "
            << "mechanicalLaw" << nl
            << "Consider using one of the other linearGeometry solidModels."
            << abort(FatalError);
    }

    // Cast the mechanical law to a linearElastic mechanicalLaw
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set mu and lambda fields
    muf_ = mech.mu();
    lambdaf_ = mech.lambda();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField>
coupledUnsLinGeomLinearElasticSolid::faceZonePointDisplacementIncrement
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

            if (zoneMeshPoints[localPoint] < mesh().nPoints())
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


tmp<tensorField>
coupledUnsLinGeomLinearElasticSolid::faceZoneSurfaceGradientOfVelocity
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
        volToPoint_.interpolate(mesh().boundaryMesh()[patchID], U_);

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


tmp<vectorField> coupledUnsLinGeomLinearElasticSolid::currentFaceZonePoints
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


tmp<vectorField> coupledUnsLinGeomLinearElasticSolid::faceZoneNormal
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


void coupledUnsLinGeomLinearElasticSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != blockSolidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn
        (
            "void coupledUnsLinGeomLinearElasticSolid::setTraction(...)"
        )   << "Bounary condition on " << D_.name()
            <<  " is "
            << D_.boundaryField()[patchID].type()
            << " for patch" << mesh().boundary()[patchID].name()
            << ", instead it should be "
            << blockSolidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    blockSolidTractionFvPatchVectorField& patchU =
        refCast<blockSolidTractionFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.traction() = traction;
}


void coupledUnsLinGeomLinearElasticSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != blockSolidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn
        (
            "void coupledUnsLinGeomLinearElasticSolid::setTraction(...)"
        )   << "Bounary condition on " << D_.name()
            <<  " is "
            << D_.boundaryField()[patchID].type()
            << "for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << blockSolidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    blockSolidTractionFvPatchVectorField& patchU =
        refCast<blockSolidTractionFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.pressure() = pressure;
}


void coupledUnsLinGeomLinearElasticSolid::setTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setTraction(patchID, patchTraction);
}


void coupledUnsLinGeomLinearElasticSolid::setPressure
(
    const label patchID,
    const label zoneID,
    const scalarField& faceZonePressure
)
{
    scalarField patchPressure(mesh().boundary()[patchID].size(), 0.0);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchPressure, i)
    {
        patchPressure[i] =
            faceZonePressure
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setPressure(patchID, patchPressure);
}


bool coupledUnsLinGeomLinearElasticSolid::evolve()
{
    Info << "Evolving solid solver" << endl;

    BlockSolverPerformance<vector> solverPerfD("undefined", "blockD");
    BlockLduMatrix<vector>::debug = 5;

    Info<< "Solving the momentum equation for D" << endl;

    // Note: there is no momentum loop as everything is linear and we have
    // hard-coded the linear definition of stress

    // Global coefficients are currently stored in the extendedMesh so we
    // must clear them out each time before constructing a new equation
    // Using a global patch may be a nicer solution
    //extendedMesh_.clearOut();
    extendedMesh_.clearOutGlobalCoeffs();

    // Store fields for under-relaxation and residual calculation
    D_.storePrevIter();

    // Create source vector for block matrix
    vectorField blockB(solutionVec_.size(), vector::zero);

    // Create block system
    BlockLduMatrix<vector> blockM(extendedMesh_);

    // Grab block diagonal and set it to zero
    Field<tensor>& d = blockM.diag().asSquare();
    d = tensor::zero;

    // Grab linear off-diagonal and set it to zero
    Field<tensor>& l = blockM.lower().asSquare();
    Field<tensor>& u = blockM.upper().asSquare();
    u = tensor::zero;
    l = tensor::zero;

    // Insert coefficients

    // For now we create separate matrices for each term of the three
    // diffusion terms in the momentum equation and add the contributions
    // together manually; a BlockFvMatrix class would help make this look
    // nicer

    // Laplacian
    // non-orthogonal correction is treated implicitly
    BlockLduMatrix<vector> blockMatLap =
        BlockFvm::laplacian(extendedMesh_, muf_, D_, blockB);

    // Laplacian transpose == div(mu*gradU.T())
    BlockLduMatrix<vector> blockMatLapTran =
        BlockFvm::laplacianTranspose(extendedMesh_, muf_, D_, blockB);

    // Laplacian trace == div(lambda*I*tr(gradU))
    BlockLduMatrix<vector> blockMatLapTrac =
        BlockFvm::laplacianTrace(extendedMesh_, lambdaf_, D_, blockB);

    // Add diagonal contributions
    d += blockMatLap.diag().asSquare();
    d += blockMatLapTran.diag().asSquare();
    d += blockMatLapTrac.diag().asSquare();

    // Add off-diagonal contributions
    u += blockMatLap.upper().asSquare();
    u += blockMatLapTran.upper().asSquare();
    u += blockMatLapTrac.upper().asSquare();
    l += blockMatLap.lower().asSquare();
    l += blockMatLapTran.lower().asSquare();
    l += blockMatLapTrac.lower().asSquare();

    // Add contribution for processor boundaries
    blockM.interfaces() = blockMatLap.interfaces();
    forAll(mesh().boundaryMesh(), patchI)
    {
        const word& patchType = mesh().boundaryMesh()[patchI].type();
        if (patchType == processorPolyPatch::typeName)
        {
            Field<tensor>& coupleUpper =
                blockM.coupleUpper()[patchI].asSquare();

            coupleUpper = blockMatLap.coupleUpper()[patchI].asSquare();
            coupleUpper +=
                blockMatLapTran.coupleUpper()[patchI].asSquare();
            coupleUpper +=
                blockMatLapTrac.coupleUpper()[patchI].asSquare();
        }
    }

    // We manually add the boundary conditions equations
    // More thinking is required to get it to fit cleanly with the rest of
    // the block coupled machinery
    extendedMesh_.insertBoundaryConditions
    (
        blockM, blockB, muf_, lambdaf_, D_
    );

    // Add terms temporal and gravity terms to the block matrix and source
    extendedMesh_.addFvMatrix
    (
        blockM,
        blockB,
        rho_*fvm::d2dt2(D_) - rho_*g_,
        true
    );

    // Under-relax the linear system
    if (mesh().solutionDict().relax("DEqn"))
    {
        FatalError
            << "Equation under-relaxation disabled!"
            << abort(FatalError);
        //Info<< "Under-relaxing the equation" << endl;
        //blockM.relax
        //(
        //    solutionVec_,
        //    blockB,
        //    mesh().solutionDict().relaxationFactor("DEqn")
        //);
    }

    // Block coupled solver call
    solverPerfD =
        BlockLduSolver<vector>::New
        (
            D_.name(),
            blockM,
            mesh().solutionDict().solver("blockD")
        )->solve(solutionVec_, blockB);

    // Transfer solution vector to D field
    extendedMesh_.copySolutionVector(solutionVec_, D_);

    // Under-relax the field
    D_.relax();

    // Update gradient of displacement
    volToPoint_.interpolate(D_, pointD_);

    // Enforce zero normal displacement on symmetry
    // Todo: we should create a blockSymmetry boundary condition
    // pointField& pointDI = pointD_.internalField();
    // forAll(mesh().boundaryMesh(), patchI)
    // {
    //     if
    //     (
    //         D_.boundaryField()[patchI].type()
    //      == blockFixedDisplacementZeroShearFvPatchVectorField::typeName
    //     )
    //     {
    //         WarningIn("coupledUnsLinGeomLinearElasticSolid::evolve()")
    //             << "Setting the normal displacement on patch "
    //             << mesh().boundaryMesh()[patchI].name() << " to zero"
    //             << endl;

    //         const labelList& meshPoints =
    //             mesh().boundaryMesh()[patchI].meshPoints();
    //         const vectorField& pointNormals =
    //             mesh().boundaryMesh()[patchI].pointNormals();

    //         // Remove normal component
    //         forAll(meshPoints, pI)
    //         {
    //             const label pointID = meshPoints[pI];
    //             pointDI[pointID] =
    //                 ((I - sqr(pointNormals[pI])) & pointDI[pointID]);
    //         }
    //     }
    // }
    gradD_ = fvc::grad(D_, pointD_);

    // We will call fvc::grad a second time as fixed displacement boundaries
    // need the patch internal field values to calculate snGrad, which
    // is then used to set snGrad on the gradD boundary
    D_.correctBoundaryConditions();
    gradD_ = fvc::grad(D_, pointD_);

    // Update the strain
    epsilon_ = symm(gradD_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma_);

    // Velocity
    U_ = fvc::ddt(D_);

    return true;
}


// void coupledUnsLinGeomLinearElasticSolid::predict()
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


tmp<vectorField> coupledUnsLinGeomLinearElasticSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& mu = muf_.boundaryField()[patchID];
    const scalarField& lambda = lambdaf_.boundaryField()[patchID];

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
              - (n & (sigma - (2.0*mu + lambda)*gradD))
            )/(2.0*mu + lambda)
        )
    );
}


void coupledUnsLinGeomLinearElasticSolid::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void coupledUnsLinGeomLinearElasticSolid::writeFields(const Time& runTime)
{
    // Update equivalent strain
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
        sqrt((2.0/3.0)*magSqr(dev(epsilon_)))
    );

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
        << endl;

    // Update equivalent (von Mises) stress
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


void coupledUnsLinGeomLinearElasticSolid::end()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
