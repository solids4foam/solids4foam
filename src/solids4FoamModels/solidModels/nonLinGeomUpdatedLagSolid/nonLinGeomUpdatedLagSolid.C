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

#include "nonLinGeomUpdatedLagSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomUpdatedLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomUpdatedLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void nonLinGeomUpdatedLagSolid::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    DD() = U()*runTime().deltaT() + 0.5*sqr(runTime().deltaT())*A_;

    // Update gradient of displacement increment
    mechanical().grad(DD(), gradDD());

    // Relative deformation gradient
    relF_ = I + gradDD().T();

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian (Jacobian of relative deformation gradient)
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


void nonLinGeomUpdatedLagSolid::enforceTractionBoundaries
(
    surfaceVectorField& traction,
    const volVectorField& DD,
    const surfaceVectorField& n
) const
{
    // Enforce traction conditions
    forAll(DD.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                DD.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    DD.boundaryField()[patchI]
                );

            const vectorField& nPatch = n.boundaryField()[patchI];

            traction.boundaryFieldRef()[patchI] =
                tracPatch.traction() - nPatch*tracPatch.pressure();
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                DD.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                DD.boundaryField()[patchI]
            )
        )
        {
            // Unit normals
            const vectorField& nPatch = n.boundaryField()[patchI];

            // Set shear traction to zero
            traction.boundaryFieldRef()[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
        }
    }
}


bool nonLinGeomUpdatedLagSolid::evolveImplicitSegregated()
{
   Info<< "Evolving solid solver using an implicit segregated approach"
       << endl;

    if (predictor_ && newTimeStep())
    {
        predict();
    }

    int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solverPerfDD;
    SolverPerformance<vector>::debug = 0;
#else
    lduSolverPerformance solverPerfDD;
    blockLduMatrix::debug = 0;
#endif

    Info<< "Solving the updated Lagrangian form of the momentum equation for DD"
        << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Momentum equation incremental updated Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho_, DD())
          + fvc::d2dt2(rho_, D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div(relJ_*relFinv_ & sigma(), "div(sigma)")
          + rho_*g()
          + stabilisation().stabilisation(DD(), gradDD(), impK_)
        );

        // Under-relax the linear system
        DDEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DDEqn);

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field using fixed or adaptive under-relaxation
        relaxField(DD(), iCorr);

        // Update the total displacement
        D() = D().oldTime() + DD();

        // Update gradient of displacement increment
        mechanical().grad(DD(), gradDD());

        // Relative deformation gradient
        relF_ = I + gradDD().T();

        // Inverse relative deformation gradient
        relFinv_ = inv(relF_);

        // Total deformation gradient
        F_ = relF_ & F_.oldTime();

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJ_ = det(relF_);

        // Jacobian of deformation gradient
        J_ = relJ_*J_.oldTime();

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DDEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
       !converged
        (
            iCorr,
#ifdef OPENFOAM_NOT_EXTEND
            mag(solverPerfDD.initialResidual()),
            cmptMax(solverPerfDD.nIterations()),
#else
            solverPerfDD.initialResidual(),
            solverPerfDD.nIterations(),
#endif
            DD()
        ) && ++iCorr < nCorr()
    );

    // Update gradient of total displacement
    // Do we need this?
    gradD() = fvc::grad(D().oldTime() + DD());

    // Total displacement
    D() = D().oldTime() + DD();

    // Interpolate cell displacement increments to vertices
    mechanical().interpolate(DD(), gradDD(), pointDD());

    // Total displacement at points
    pointD() = pointD().oldTime() + pointDD();

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


bool nonLinGeomUpdatedLagSolid::evolveSnes()
{
    Info<< "Solving the momentum equation for DD using PETSc SNES" << endl;

    // Update D boundary conditions
    DD().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Use the segregated solver as a predictor
        //evolveImplicitSegregated();

        // Map the DD field to the SNES solution vector
        foamPetscSnesHelper::InsertFieldComponents<vector>
        (
            DD().primitiveFieldRef(),
            foamPetscSnesHelper::solution(),
            solidModel::twoD() ? 2 : 3, // Block size of x
            0                           // Location of first component
        );
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the DD field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        DD().primitiveFieldRef(),
        solidModel::twoD() ? 2 : 3, // Block size of x
        0                           // Location of first component
    );

    // Total displacement
    D() = D().oldTime() + DD();

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD(), gradDD(), pointDD());
    pointDD().correctBoundaryConditions();

    // Total point displacement
    pointD() = pointD().oldTime() + pointDD();;

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomUpdatedLagSolid::nonLinGeomUpdatedLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solidModel::twoD() ? 2 : 3,
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
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
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    relF_
    (
        IOobject
        (
            "relF",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        I + gradDD().T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(relF_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mechanical().rho()
    ),
    A_
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimVelocity/dimTime, vector::zero)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false))
{
    DDisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(rho_, DD());
    fvc::d2dt2(rho_, D().oldTime());

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + DD().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + DD().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << DD().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // For consistent restarts, we will update the relative kinematic fields
    DD().correctBoundaryConditions();
    if (restart())
    {
        mechanical().grad(DD(), gradDD());
        relF_ = I + gradDD().T();
        relFinv_ = inv(relF_);
        relJ_ = det(relF_);

        F_.storeOldTime();
        J_.storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }

    // Check the gradScheme
    const word gradDScheme
    (
        mesh().gradScheme("grad(" + DD().name() +')')
    );

    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        if (gradDScheme != "leastSquaresS4f")
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithm" << abort(FatalError);
        }

        // Set extrapolateValue to true for solidTraction boundaries
        forAll(DD().boundaryField(), patchI)
        {
            if
            (
                isA<solidTractionFvPatchVectorField>
                (
                    DD().boundaryField()[patchI]
                )
            )
            {
                Info<< "    Setting `extrapolateValue` to `true` on the "
                    << mesh().boundary()[patchI].name() << " patch of the D "
                    << "field" << endl;

                solidTractionFvPatchVectorField& tracPatch =
                    refCast<solidTractionFvPatchVectorField>
                    (
                        DD().boundaryFieldRef()[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomUpdatedLagSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    // else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    // {
    //     // Not yet implmented, although coupledUnsLinGeomLinearElasticSolid
    //     // could be combined with PETSc to achieve this.. todo!
    //     return evolveImplicitCoupled();
    // }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    // else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    // {
    //     return evolveExplicit();
    // }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            // << solutionAlgorithmNames_.names() << endl;
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
               ]
            // << solidModel::solutionAlgorithmNames_
            //    [
            //        solidModel::solutionAlgorithm::EXPLICIT
            //    ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


label nonLinGeomUpdatedLagSolid::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
)
{
    // Copy x into the DD field
    volVectorField& DD = const_cast<volVectorField&>(this->DD());
    vectorField& DDI = DD;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DDI,
        0,                          // Location of first DDI component
        solidModel::twoD() ? 2 : 3  // Number of components to extract
    );

    // Enforce the boundary conditions
    DD.correctBoundaryConditions();

    // Update total displacement
    D() = D().oldTime() + DD;

    // Update displacement increment gradient
    mechanical().grad(DD, gradDD());

    // Relative deformation gradient
    relF_ = I + gradDD().T();

    // Inverse relative deformation gradient
    relFinv_ = inv(relF_);

    // Total deformation gradient
    F_ = relF_ & F_.oldTime();

    // Relative Jacobian (Jacobian of relative deformation gradient)
    relJ_ = det(relF_);

    // Jacobian of deformation gradient
    J_ = relJ_*J_.oldTime();

    // Update the momentum equation inverse diagonal field
    // This may be used by the mechanical law when calculating the
    // hydrostatic pressure
    //const volScalarField DEqnA("DEqnA", DDEqn.A());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());
    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(relJ_*relFinv_.T()) & mesh().Sf()
    );
    const surfaceScalarField magSfCurrent(mag(SfCurrent));
    const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

    // Traction vectors at the faces
    //surfaceVectorField traction(n & fvc::interpolate(sigma()));
    surfaceVectorField traction(nCurrent & fvc::interpolate(sigma()));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDDf(fvc::interpolate(gradDD()));
    traction += scaleFactor*impKf_*(fvc::snGrad(DD) - (n & gradDDf));

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, DD, nCurrent);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual
    (
        fvc::div(magSfCurrent*traction)
      + rho()
       *(
           g() - fvc::d2dt2(D()) - dampingCoeff()*fvc::ddt(D())
        )
    );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh().V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    //residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0,                          // Location of first DI component
        solidModel::twoD() ? 2 : 3  // Number of components to extract
    );

    return 0;
}


label nonLinGeomUpdatedLagSolid::formJacobian
(
    Mat jac,
    const PetscScalar *x
)
{
    // Copy x into the DD field
    volVectorField& DD = const_cast<volVectorField&>(this->DD());
    vectorField& DDI = DD;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DDI,
        0,                          // Location of first DDI component
        solidModel::twoD() ? 2 : 3  // Number of components to extract
    );

    // Enforce the boundary conditions
    DD.correctBoundaryConditions();

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, DD, "laplacian(DDD,DD)")
      - rho()*fvm::d2dt2(DD)
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(DD);
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    // Convert fvMatrix matrix to PETSc matrix
    foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
    (
        approxJ, jac, 0, 0, solidModel::twoD() ? 2 : 3
    );

    return 0;
}


tmp<vectorField> nonLinGeomUpdatedLagSolid::tractionBoundarySnGrad
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
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch relative deformation gradient inverse
    const tensorField& relFinv = relFinv_.boundaryField()[patchID];

    // Patch unit normals (updated configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(relFinv.T() & n);
    nCurrent /= mag(nCurrent);

    // Testing: let us instead calculate the deformed normals by interpolating
    // displacements to the points and calculating the normals on the deformed
    // patch; as this is how we will actually move the mesh, it will be more
    // consistent.
    // This, however, begs the question: is the cell-centred deformation
    // gradient field 'F' consistent with our point displacement field?"
    // i.e. we can calculate the deformed cell volumes two ways (at least):
    //     1. V = J*Vold
    //     2. Move the mesh with pointD and then directly calculate V
    // The answers from 1. and 2. are only approximately equal: this causes a
    // slight inconsistency. The equalavent can be said for the deformed face
    // areas.
    // In Maneeratana, the mesh is never moved, instead method 1. is used for
    // the deformed volumes and areas.

    // standAlonePatch deformedPatch =
    //     standAlonePatch
    //     (
    //         mesh().boundaryMesh()[patchID].localFaces(),
    //         mesh().boundaryMesh()[patchID].localPoints()
    //     );

    // // Calculate the deformed points
    // const pointField deformedPoints =
    //     mechanical().volToPoint().interpolate
    //     (
    //         mesh().boundaryMesh()[patchID],
    //         DD_
    //     )
    //   + mesh().boundaryMesh()[patchID].localPoints();

    // // Move the standAlonePatch points
    // const_cast<pointField&>(deformedPatch.points()) = deformedPoints;

    // // Patch unit normals (deformed configuration)
    // const vectorField& nCurrent = deformedPatch.faceNormals();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradDD)
            )*rImpK
        )
    );
}


void nonLinGeomUpdatedLagSolid::updateTotalFields()
{
    // Density
    rho_ = rho_.oldTime()/relJ_;

    // Move the mesh to the deformed configuration
#ifdef OPENFOAM_NOT_EXTEND
    const vectorField oldPoints = mesh().points();
#else
    const vectorField oldPoints = mesh().allPoints();
#endif
    moveMesh(oldPoints, pointDD());

    solidModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
