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


void unsNonLinGeomUpdatedLagSolid::moveMesh(const pointField& oldPoints)
{
    Info<< "Moving the mesh to the deformed configuration" << nl << endl;

    //- Move mesh by interpolating displacement field to vertices
    // To be checked: sync boundary and global points across procs to make sure
    // numiercal error does not build up and when end up with the error
    // "face area does not match neighbour..."
    // We could sync points as a pointVectorField just as we sync pointDD

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD(), pointDD());

    // Ensure continuous displacement across processor boundary
    // Something strange is happening here
    pointDD().correctBoundaryConditions();

    vectorField& pointDDI = pointDD().internalField();

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
    twoDCorrector.correctPoints(pointDD().internalField());
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
    gradDDf_
    (
        IOobject
        (
            "grad(" + DD().name() + ")f",
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
        I + gradDD().T()
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
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mechanical().rho()
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    stabilisePressure_(lookupOrDefault<Switch>("stabilisePressure", false))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


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
        DD().storePrevIter();

        // Momentum equation incremental updated Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho_, DD())
          + fvc::d2dt2(rho_.oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div((relJf_*relFinvf_.T() & mesh().Sf()) & sigmaf_)
          + rho()*g()
        );

        // Under-relax the linear system
        DDEqn.relax();

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field
        relaxField(DD(), iCorr);

        // Interpolate DD to pointDD
        mechanical().interpolate(DD(), pointDD(), false);

        // Update gradient of displacement increment
        mechanical().grad(DD(), pointDD(), gradDD());
        mechanical().grad(DD(), pointDD(), gradDDf_);

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
    while
    (
       !converged
        (
            iCorr,
            solverPerfDD.initialResidual(),
            solverPerfDD.nIterations(),
            DD()
        ) && ++iCorr < nCorr()
    );

    // Relative deformation gradient
    relF_ = I + gradDD().T();

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Update gradient of total displacement
    gradD() = fvc::grad(D().oldTime() + DD());

    // Total displacement
    D() = D().oldTime() + DD();

    // Total displacement at points
    pointD() = pointD().oldTime() + pointDD();

    // Velocity
    U() = fvc::ddt(D());

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

    solidModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //