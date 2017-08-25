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

#include "coupledUnsLinGeomSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"
#include "BlockFvmDivSigma.H"

#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledUnsLinGeomSolid, 0);
addToRunTimeSelectionTable(physicsModel, coupledUnsLinGeomSolid, solid);
addToRunTimeSelectionTable(solidModel, coupledUnsLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool coupledUnsLinGeomSolid::converged
(
    const int iCorr,
    const BlockSolverPerformance<vector>& solverPerfD
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
            mag(solverPerfD.initialResidual()) < solutionTol_
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
            mag(solverPerfD.initialResidual()) < alternativeTol_
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
            << ", " << mag(solverPerfD.finalResidual())
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledUnsLinGeomSolid::coupledUnsLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    extendedMesh_(mesh()),
    solutionVec_
    (
        IOobject
        (
            "solutionVec",
            runTime.timeName(),
            mesh(),
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
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonf_
    (
        IOobject
        (
            "epsilonf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    volToPoint_(mesh()),
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
    gradDf_
    (
        IOobject
        (
            "grad(" + D_.name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    rho_(mechanical().rho()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    muf_("muf", impKf_/3.5), // assuming a Poisson's ratio of 0.3
    lambdaf_("lambdaf", 1.5*muf_),
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
        solidProperties().lookupOrDefault<int>("infoFrequency", 1)
    ),
    nCorr_(solidProperties().lookupOrDefault<int>("nCorrectors", 100)),
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


tmp<vectorField> coupledUnsLinGeomSolid::faceZonePointDisplacementIncrement
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


tmp<tensorField> coupledUnsLinGeomSolid::faceZoneSurfaceGradientOfVelocity
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


tmp<vectorField> coupledUnsLinGeomSolid::currentFaceZonePoints
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


tmp<vectorField> coupledUnsLinGeomSolid::faceZoneNormal
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


void coupledUnsLinGeomSolid::setTraction
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
        FatalErrorIn("void coupledUnsLinGeomSolid::setTraction(...)")
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


void coupledUnsLinGeomSolid::setPressure
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
        FatalErrorIn("void coupledUnsLinGeomSolid::setTraction(...)")
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


void coupledUnsLinGeomSolid::setTraction
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


void coupledUnsLinGeomSolid::setPressure
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


bool coupledUnsLinGeomSolid::evolve()
{
    Info << "Evolving solid solver" << endl;

    int iCorr = 0;
    BlockSolverPerformance<vector> solverPerfD("undefined", "blockD");
    BlockLduMatrix<vector>::debug = 5;

    Info<< "Solving the momentum equation for D" << endl;

    // Momentum equation loop
    // In this case, a Hookean elastic constutive equation is assumed, and outer
    // corrections are performed is a different (possibily nonlinear )definition
    // of stress is used
    do
    {
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


        // TESTING start
        // The implicit component assumes small strain Hookean elastic, so we
        // will explicitly remove this and then explicitly add on the divergence
        // of the actual stress
        /*
        // Store the div(sigma) coeffs, cells and boundary conditions
        vectorField AxDivSigma = vectorField(solutionVec_.size(), vector::zero);
        blockM.Amul(AxDivSigma, solutionVec_);

        // Store the div(sigma) coeffs, just cells
        vectorField AxDivSigmaCells =
            SubField<vector>(AxDivSigma, mesh().nCells(), 0);
            */
        // We will leave the boundary conditions untouched to be will set the
        // contribution to zero
        // for (int varI = mesh().nCells(); varI < Ax.size(); varI++)
        // {
        //     Ax[varI] = vector::zero;
        // }

        // Remove contribution to internal cells
        //blockB -= Ax;
        // TESTING end



        // Add terms temporal and gravity terms to the block matrix and source
        // Alos, we explicitly remove the implicit Hooke's law definition of
        // stress and explicitly add the definition of stress from the
        // mechanicalLaw
        extendedMesh_.addFvMatrix
        (
            blockM,
            blockB,
            rho_*fvm::d2dt2(D_) - rho_*g_,
            true
        );

        /*
        // Store the div(sigma) coeffs, cells and boundary conditions
        vectorField AxD2dt2RhoG =
            vectorField(solutionVec_.size(), vector::zero);
        blockM.Amul(AxD2dt2RhoG, solutionVec_);
        AxD2dt2RhoG -= AxDivSigma;

        // Store the just cell equations
        vectorField AxD2dt2RhoGCells =
            SubField<vector>(AxD2dt2RhoG, mesh().nCells(), 0);


        // TESTING start
        // Explicitly add the real definition of stress
        volVectorField newDivSigma =
            volVectorField("newDivSigma", fvc::div(mesh().Sf() & sigmaf_));
        newDivSigma.internalField() *= mesh().V();
        newDivSigma.write();
        const surfaceVectorField force = mesh().Sf() & sigmaf_;
        Info<< "Force[0] is: " << max(mag(force.boundaryField()[0])) << nl
            << "Force[1] is: " << max(mag(force.boundaryField()[1])) << endl;

        vectorField Ax = vectorField(solutionVec_.size(), vector::zero);
        blockM.Amul(Ax, solutionVec_);
        volVectorField AxField("Ax", 0.0*newDivSigma);
        AxField.internalField() = SubField<vector>(Ax, mesh().nCells(), 0);
        AxField.write();

        volVectorField blockBField("blockB", 0.0*newDivSigma);
        blockBField.internalField() =
            SubField<vector>(blockB, mesh().nCells(), 0);
        blockBField.write();

        // volVectorField AxD2dt2RhoGCellsField("AxD2dt2RhoG", 0.0*newDivSigma);
        // AxD2dt2RhoGCellsField.internalField() = AxD2dt2RhoGCells;
        // AxD2dt2RhoGCellsField.write();

        volVectorField AxDivSigmaCellsField("AxDivSigmaCells", 0.0*newDivSigma);
        AxDivSigmaCellsField.internalField() = AxDivSigmaCells;
        AxDivSigmaCellsField.write();

        volVectorField RhoTerm("RhoTerm", 0.0*newDivSigma);
        RhoTerm.internalField() = rho_*fvc::d2dt2(D_) - rho_*g_;
        RhoTerm.internalField() *= mesh().V();
        RhoTerm.write();

        // Minus means we add to the right-hand side
        volVectorField contrib =
            volVectorField("contrib", AxDivSigmaCellsField - newDivSigma);
        const unallocLabelList& faceCells0 =
            mesh().boundaryMesh()[0].faceCells();
        forAll(contrib.boundaryField()[0], faceI)
        {
            contrib.internalField()[faceCells0[faceI]] = vector::zero;
        }
        const unallocLabelList& faceCells1 =
            mesh().boundaryMesh()[1].faceCells();
        forAll(contrib.boundaryField()[1], faceI)
        {
            contrib.internalField()[faceCells1[faceI]] = vector::zero;
        }
        contrib.write();
        // extendedMesh_.addFvSource
        // (
        //     blockB,
        //     contrib
        // );
*/
        // Block coupled solver call
        // solverPerfD =
        //     BlockLduSolver<vector>::New
        //     (
        //         D_.name(),
        //         blockM,
        //         mesh().solutionDict().solver("blockD")
        //     )->solve(solutionVec_, blockB);

        // Under-relax the linear system
        if (mesh().solutionDict().relax("DEqn"))
        {
            FatalError
                << "Equation under-relaxation disabled"
                << abort(FatalError);

            // Info<< "Under-relaxing the equation" << endl;
            // blockM.relax
            // (
            //     solutionVec_,
            //     blockB,
            //     mesh().solutionDict().relaxationFactor("DEqn")
            // );
        }

        // Create the linear solver
        autoPtr<BlockLduSolver<vector> > solver =
            BlockLduSolver<vector>::New
            (
                D_.name(),
                blockM,
                mesh().solutionDict().solver("blockD")
            );

        // Solve the linear system
        solver->solve(solutionVec_, blockB);

        // Transfer solution vector to D field
        extendedMesh_.copySolutionVector(solutionVec_, D_);

        // Under-relax the field
        D_.relax();

        // Update gradient of displacement
        volToPoint_.interpolate(D_, pointD_);
        gradD_ = fvc::grad(D_, pointD_);
        gradDf_ = fvc::fGrad(D_, pointD_);

        // We will call fvc::grad a second time as fixed displacement boundaries
        // need the patch internal field values to calculate snGrad, which
        // is then used to set snGrad on the gradD boundary
        D_.correctBoundaryConditions();
        gradD_ = fvc::grad(D_, pointD_);
        gradDf_ = fvc::fGrad(D_, pointD_);

        // Update the strain
        epsilonf_ = symm(gradDf_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);

        // TESTING
        mechanical().correct(sigma_);
    }
    while (!converged(iCorr, solverPerfD) && ++iCorr < nCorr_);

    // Calculate cell strain
    epsilon_ = symm(gradD_);

    // Calculate cell stress
    mechanical().correct(sigma_);

    // Velocity
    U_ = fvc::ddt(D_);

    return true;
}


// void coupledUnsLinGeomSolid::predict()
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


tmp<vectorField> coupledUnsLinGeomSolid::tractionBoundarySnGrad
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
    const tensorField& gradD = gradDf_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigmaf_.boundaryField()[patchID];

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


void coupledUnsLinGeomSolid::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void coupledUnsLinGeomSolid::writeFields(const Time& runTime)
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


void coupledUnsLinGeomSolid::end()
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
