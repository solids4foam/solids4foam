/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_COM

#include "meshMotionTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(meshMotionTotalDispSolid, 0);
addToRunTimeSelectionTable(solidModel, meshMotionTotalDispSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void meshMotionTotalDispSolid::updateKinematics()
{
    // Gradient of the total displacement on the current mesh
    gradD() = fvc::grad(D());

    if (!finiteStrain_)
    {
        return;
    }

    // The current points are x = X + D, so the inverse deformation gradient
    // is dX/dx = I - gradD^T
    Finv_ = I - gradD().T();

    const scalar minDetFinv = gMin(det(Finv_.primitiveField())());

    if (minDetFinv < SMALL)
    {
        FatalErrorInFunction
            << "The mesh motion displacement inverts the mesh: the minimum "
            << "determinant of the inverse deformation gradient is "
            << minDetFinv << abort(FatalError);
    }

    F_ = inv(Finv_);
    J_ = det(F_);
}


void meshMotionTotalDispSolid::correctStress()
{
    if (finiteStrain_)
    {
        // Cauchy stress from the total deformation gradient
        mechanicalManager().updateStressFiniteStrain
        (
            F_,
            F_.oldTime(),
            Finv_,
            Finv_.oldTime(),
            J_,
            J_.oldTime(),
            mesh().time().deltaTValue(),
            sigma()
        );
    }
    else
    {
        mechanicalManager().updateStressSmallStrain
        (
            gradD(),
            gradD().oldTime(),
            mesh().time().deltaTValue(),
            sigma()
        );
    }
}


void meshMotionTotalDispSolid::updateImpK()
{
    // A tangent query, so it neither writes a stress nor disturbs history
    if (finiteStrain_)
    {
        mechanicalManager().updateScalarTangentFiniteStrain
        (
            F_,
            F_.oldTime(),
            Finv_,
            Finv_.oldTime(),
            J_,
            J_.oldTime(),
            mesh().time().deltaTValue(),
            impK_,
            tangentRequest::scalar
        );
    }
    else
    {
        mechanicalManager().updateScalarTangent
        (
            gradD(),
            gradD().oldTime(),
            mesh().time().deltaTValue(),
            impK_,
            tangentRequest::scalar
        );
    }

    impKf_ = fvc::interpolate(impK_);
    rImpK_ = 1.0/impK_;
}


void meshMotionTotalDispSolid::enforceTractionBoundaries
(
    surfaceVectorField& traction,
    const surfaceVectorField& n
) const
{
    forAll(D().boundaryField(), patchI)
    {
        const fvPatchVectorField& DPatch = D().boundaryField()[patchI];
        const vectorField& nPatch = n.boundaryField()[patchI];

        if (isA<solidTractionFvPatchVectorField>(DPatch))
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>(DPatch);

            boundaryFieldRef(traction)[patchI] =
                tracPatch.traction() - nPatch*tracPatch.pressure();
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>(DPatch)
         || isA<symmetryFvPatchVectorField>(DPatch)
         || isA<slipFvPatchVectorField>(DPatch)
        )
        {
            // Set the shear traction to zero
            boundaryFieldRef(traction)[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
        }
    }
}


tmp<surfaceVectorField> meshMotionTotalDispSolid::faceForce()
{
    // Unit normal vectors at the faces of the current mesh
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma()));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the
    // stabilisation is set to zero on traction boundaries
    momentumStabilisation().updateVector(D(), &gradD());
    traction += impKf_*momentumStabilisation().faceVector();

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, n);

    return tmp<surfaceVectorField>
    (
        new surfaceVectorField("faceForce", mesh().magSf()*traction)
    );
}


void meshMotionTotalDispSolid::updateDerivedFields()
{
    // Increments of displacement
    DD() = D() - D().oldTime();
    gradDD() = gradD() - gradD().oldTime();

    // Interpolate cell displacements to vertices
    interpolatePointDisplacement(D(), gradD(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();
}


bool meshMotionTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving the mesh motion solid model on the current mesh" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // The mesh has moved since the last solution, so the kinematics, the
    // stiffness and the stress are recalculated on the current mesh
    updateKinematics();
    updateImpK();
    correctStress();

    SolverPerformance<vector>::debug = 0;

    int iCorr = 0;
    scalar currentResidualNorm = 0;
    scalar initialResidualNorm = 0;
    scalar deltaXNorm = 0;
    scalar xNorm = 0;
    const convergenceParameters convParam =
        readConvergenceParameters(solidModelDict());

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Quasi-static momentum equation on the current mesh
        fvVectorMatrix DEqn
        (
            fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(faceForce())
        );

        // Under-relax the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system and store the residual
        currentResidualNorm = mag(DEqn.solve().initialResidual());

        // Norm of the solution correction
        deltaXNorm =
            sqrt
            (
                gSum
                (
                    magSqr
                    (
                        D().primitiveField() - D().prevIter().primitiveField()
                    )
                )
            );

        // Norm of the solution
        xNorm = sqrt(gSum(magSqr(D().primitiveField())));

        // Store the initial residual
        if (iCorr == 0)
        {
            initialResidualNorm = currentResidualNorm;
            Info<< "Initial Residual Norm = " << initialResidualNorm << nl
                << "Initial Solution Norm = " << xNorm << endl;
        }

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Update the kinematics and the stress
        updateKinematics();
        correctStress();
    }
    while
    (
        !checkConvergence
        (
            currentResidualNorm,
            initialResidualNorm,
            deltaXNorm,
            xNorm,
            ++iCorr,
            convParam
        )
    );

    updateDerivedFields();

    SolverPerformance<vector>::debug = 1;

    return true;
}


bool meshMotionTotalDispSolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Evolving the mesh motion solid model on the current mesh using "
        << "PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // The mesh has moved since the last solution, so the kinematics, the
    // stiffness and the stress are recalculated on the current mesh, as in
    // the segregated loop. This stops if the boundary motion inverts the
    // mesh. The implicit stiffness is then fixed during the solution, so the
    // residual is that of the segregated loop
    updateKinematics();
    updateImpK();
    correctStress();

    // Seed the SNES solution vector from the current D on every solution:
    // the initial guess is the last solution, with the new boundary values,
    // as in the segregated loop, and the owner may have set D externally,
    // e.g. on a time-step retry
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        primitiveFieldRef(D()),
        foamPetscSnesHelper::solution(),
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Map the PETSc solution to the D field
    vectorField& DI = D();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    D().correctBoundaryConditions();

    // Update the kinematics and the stress from the solution, which also
    // stops if the solution inverts the mesh
    updateKinematics();
    correctStress();

    updateDerivedFields();

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshMotionTotalDispSolid::meshMotionTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "D",
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solutionLocation::CELLS,
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
    ),
    finiteStrain_
    (
        solidModelDict().lookupOrDefault<Switch>("finiteStrain", false)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
        mesh(),
        dimensionedScalar("one", dimless, 1)
    ),
    impK_
    (
        IOobject
        (
            "impK",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimForce/dimArea, 0),
        calculatedFvPatchScalarField::typeName
    ),
    impKf_(fvc::interpolate(impK_)),
    rImpK_
    (
        IOobject
        (
            "rImpK",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimArea/dimForce, 0),
        calculatedFvPatchScalarField::typeName
    ),
    blockSize_(solidModel::twoD() ? 2 : 3)
{
    if (!sharedMesh())
    {
        FatalErrorInFunction
            << "The " << typeName << " solid model moves a mesh owned by "
            << "another model, e.g. with the solids4foamSolidModel mesh "
            << "motion solver, and is not for solid mechanics"
            << exit(FatalError);
    }

    DisRequired();

    // A multi-material run needs a material-aware gradient
    checkFrameworkGradScheme(D().name());

    // The traction is applied to the current face areas
    forAll(D().boundaryField(), patchI)
    {
        if (isA<solidTractionFvPatchVectorField>(D().boundaryField()[patchI]))
        {
            if
            (
                refCast<const solidTractionFvPatchVectorField>
                (
                    D().boundaryField()[patchI]
                ).useUndeformedArea()
            )
            {
                FatalErrorInFunction
                    << "useUndeformedArea is not supported by the "
                    << typeName << " solid model, on patch "
                    << mesh().boundary()[patchI].name() << exit(FatalError);
            }
        }
    }

    if
    (
        solutionAlg() != solutionAlgorithm::IMPLICIT_SEGREGATED
     && solutionAlg() != solutionAlgorithm::PETSC_SNES
    )
    {
        FatalErrorInFunction
            << "The " << typeName << " solid model supports the "
            << solutionAlgorithmNames_[solutionAlgorithm::IMPLICIT_SEGREGATED]
            << " and "
            << solutionAlgorithmNames_[solutionAlgorithm::PETSC_SNES]
            << " solution algorithms, not "
            << solutionAlgorithmNames_[solutionAlg()] << exit(FatalError);
    }

    if (highOrderResidual() || highOrderJacobian())
    {
        FatalErrorInFunction
            << "highOrderResidual and highOrderJacobian are not supported by "
            << "the " << typeName << " solid model" << exit(FatalError);
    }

    Info<< "finiteStrain = " << finiteStrain_ << nl
        << "solutionAlgorithm = "
        << solutionAlgorithmNames_[solutionAlg()] << endl;

    D().correctBoundaryConditions();
    updateKinematics();
    updateImpK();
    correctStress();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool meshMotionTotalDispSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }

    return evolveImplicitSegregated();
}


#ifdef USE_PETSC

label meshMotionTotalDispSolid::initialiseJacobian(Mat& jac)
{
    // Initialise based on the compact stencil of the mesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label meshMotionTotalDispSolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label meshMotionTotalDispSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    // Copy x into the D field
    vectorField& DI = D();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    D().correctBoundaryConditions();

    // Update the kinematics on the current mesh, which stops if D inverts
    // the mesh, and the stress, as in the segregated loop
    updateKinematics();
    correctStress();

    // Enforce the boundary conditions again, now that gradD and sigma have
    // been updated from x, so that the traction conditions use the stress of
    // this x
    D().correctBoundaryConditions();

    // The residual is the divergence of the face forces, i.e. the explicit
    // part of the segregated momentum equation, made extensive
    vectorField residual(fvc::div(faceForce())().primitiveField());
    residual *= mesh().V();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    return 0;
}


label meshMotionTotalDispSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Copy x into the D field
    vectorField& DI = D();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    D().correctBoundaryConditions();

    // The implicit stiffness Laplacian of the segregated loop, on the current
    // mesh. It is built here rather than taken from the momentum
    // stabilisation, whose cached matrix would keep the geometry and the
    // stiffness of the mesh it was first built on
    fvVectorMatrix approxJ(fvm::laplacian(impKf_, D(), "laplacian(DD,D)"));

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    // Convert the fvMatrix to a PETSc matrix
    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
    (
        approxJ, jac, 0, 0, blockSize_
    );

    return 0;
}

#endif // USE_PETSC


tmp<vectorField> meshMotionTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals of the current mesh
    const vectorField n(patch.nf());

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

#endif // OPENFOAM_COM

// ************************************************************************* //
