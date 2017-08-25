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

#include "thermalLinGeomSolid.H"
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

defineTypeNameAndDebug(thermalLinGeomSolid, 0);
addToRunTimeSelectionTable(physicsModel, thermalLinGeomSolid, solid);
addToRunTimeSelectionTable(solidModel, thermalLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool thermalLinGeomSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfT,
    const lduSolverPerformance& solverPerfD
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar absResidualT =
        gMax(mag(T_.internalField() - T_.oldTime().internalField()));
    const scalar residualT =
        absResidualT
       /max
        (
            gMax(mag(T_.internalField() - T_.oldTime().internalField())),
            SMALL
        );

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
        bool convergedD = false;
        bool convergedT = false;

        if
        (
            (
                solverPerfD.initialResidual() < solutionTol_
             && residualD < solutionTol_
            )
         || solverPerfD.initialResidual() < alternativeTol_
         || residualD < alternativeTol_
        )
        {
            convergedD = true;
        }

        if
        (
            (
                solverPerfT.initialResidual() < solutionTol_
             && residualT < solutionTol_
            )
         || solverPerfT.initialResidual() < alternativeTol_
         || residualT < alternativeTol_
         || absResidualT < absTTol_
        )
        {
            convergedT = true;
        }


        if (convergedD && convergedT)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (T & D), relRes (T & D), matRes, iters (T & D)"
            << endl;
    }
    else if (iCorr % infoFrequency_ == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfT.initialResidual()
            << ", " << solverPerfD.initialResidual()
            << ", " << residualT
            << ", " << residualD
            << ", " << materialResidual
            << ", " << solverPerfT.nIterations()
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
            << "Max iterations reached within the enery-momentum loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalLinGeomSolid::thermalLinGeomSolid
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
    rho_(mechanical().rho()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    T_
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    T0_
    (
        IOobject
        (
            "T0",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("T0", dimTemperature, 0.0)
    ),
    gradT_
    (
        IOobject
        (
            "grad(T)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimTemperature/dimLength, vector::zero)
    ),
    rhoC_
    (
        IOobject
        (
            "rhoC",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    solutionTol_
    (
        solidProperties().lookupOrDefault<scalar>("solutionToleranceT", 1e-06)
    ),
    alternativeTol_
    (
        solidProperties().lookupOrDefault<scalar>
        (
            "alternativeToleranceT", 1e-06
        )
    ),
    materialTol_
    (
        solidProperties().lookupOrDefault<scalar>("materialTolerance", 1e-05)
    ),
    absTTol_
    (
        solidProperties().lookupOrDefault<scalar>
        (
            "absoluteTemperatureTolerance",
            1e-06
        )
    ),
    infoFrequency_
    (
        solidProperties().lookupOrDefault<int>("infoFrequencyT", 100)
    ),
    nCorr_(solidProperties().lookupOrDefault<int>("nCorrectorsT", 10000)),
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
    gradD_.oldTime();
    sigma_.oldTime();
    T_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> thermalLinGeomSolid::faceZonePointDisplacementIncrement
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


tmp<tensorField> thermalLinGeomSolid::faceZoneSurfaceGradientOfVelocity
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


tmp<vectorField> thermalLinGeomSolid::currentFaceZonePoints
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


tmp<vectorField> thermalLinGeomSolid::faceZoneNormal
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


void thermalLinGeomSolid::setTraction
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
        FatalErrorIn("void thermalLinGeomSolid::setTraction(...)")
            << "Boundary condition on " << D_.name()
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


void thermalLinGeomSolid::setPressure
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
        FatalErrorIn("void thermalLinGeomSolid::setTraction(...)")
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


bool thermalLinGeomSolid::evolve()
{
    Info << "Evolving thermal solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfD;
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;

    Info<< "Solving coupled energy and displacements equation for T and D"
        << endl;

    // Energy equation non-orthogonality corrector loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        T_.storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            rhoC_*fvm::ddt(T_)
         == fvm::laplacian(k_, T_, "laplacian(k,T)")
          + (sigma_ && fvc::grad(U_))
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        solverPerfT = TEqn.solve();

        // Under-relax the field
        T_.relax();

        // Update gradient of temperature
        gradT_ = fvc::grad(T_);

        // Store fields for under-relaxation and residual calculation
        D_.storePrevIter();

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho_*fvm::d2dt2(D_)
         == fvm::laplacian(impKf_, D_, "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D_, "laplacian(DD,D)")
          + fvc::div(sigma_, "div(sigma)")
          + rho_*g_
          + mechanical().RhieChowCorrection(D_, gradD_)
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the field
        D_.relax();

        // Update velocity
        U_ = fvc::ddt(D_);

        // Update gradient of displacement
        mechanical().grad(D_, gradD_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma_);

        // Update impKf to improve convergence
        // Note: impK and rImpK are not updated as they are used for traction
        // boundaries
        if (iCorr % 10 == 0)
        {
            impKf_ = mechanical().impKf();
        }
    }
    while (!converged(iCorr, solverPerfT, solverPerfD) && ++iCorr < nCorr_);

    // Interpolate cell displacements to vertices
    Info<< "Interpolating D to pointD" << endl;
    mechanical().interpolate(D_, pointD_);

    return true;
}


tmp<vectorField> thermalLinGeomSolid::tractionBoundarySnGrad
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


void thermalLinGeomSolid::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void thermalLinGeomSolid::writeFields(const Time& runTime)
{
    Info<< "Max T = " << max(T_).value() << nl
        << "Min T = " << min(T_).value() << endl;

    // Heat flux
    volVectorField heatFlux
    (
        IOobject
        (
            "heatFlux",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       -k_*gradT_
    );

    Info<< "Max magnitude of heat flux = " << max(mag(heatFlux)).value()
        << endl;

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
        sqrt((2.0/3.0)*magSqr(dev(symm(gradD_))))
    );

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
        << endl;

    // Calculate equivalent (von Mises) stress
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


void thermalLinGeomSolid::end()
{
    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum energy correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The energy equation converged in all time-steps"
            << nl << endl;
    }

    solidModel::end();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
