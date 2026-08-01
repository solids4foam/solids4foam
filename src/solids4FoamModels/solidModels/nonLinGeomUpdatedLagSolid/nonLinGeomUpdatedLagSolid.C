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
#include "slipFvPatchFields.H"
#include "compatibilityFunctions.H"
#ifndef FOAMEXTEND
    #include "hofvm.H"
#endif

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


#ifndef FOAMEXTEND
void nonLinGeomUpdatedLagSolid::makeQuadratureKinematics() const
{
    if (!FQuadOldPtr_.empty() || !gradDTotalQuadPtr_.empty())
    {
        FatalErrorInFunction
            << "pointer already set!" << abort(FatalError);
    }

    const CompactListList<point>& faceQuadPts =
        displacementMLS().quadrature().faceQuadPoints();

    labelList rowSizes(faceQuadPts.size(), 0);
    forAll(faceQuadPts, faceI)
    {
        rowSizes[faceI] = faceQuadPts[faceI].size();
    }

    FQuadOldPtr_.set(new CompactListList<tensor>(rowSizes));
    gradDTotalQuadPtr_.set(new CompactListList<tensor>(rowSizes));

    CompactListList<tensor>& FQuadOld = FQuadOldPtr_();
    CompactListList<tensor>& gradDTotalQuad = gradDTotalQuadPtr_();

    forAll(FQuadOld, faceI)
    {
        forAll(FQuadOld[faceI], qpI)
        {
            FQuadOld[faceI][qpI] = I;
            gradDTotalQuad[faceI][qpI] = tensor::zero;
        }
    }
}


void nonLinGeomUpdatedLagSolid::updateQuadratureKinematics() const
{
    if (FQuadOldPtr_.empty())
    {
        makeQuadratureKinematics();
    }

    // Displacement increment gradient at the face quadrature points, i.e. the
    // gradient relative to the configuration at the start of the time step
    const CompactListList<tensor>& gradDDQuad = gradDQuad();

    const CompactListList<tensor>& FQuadOld = FQuadOldPtr_();
    CompactListList<tensor>& gradDTotalQuad = gradDTotalQuadPtr_();

    forAll(gradDTotalQuad, faceI)
    {
        const UList<tensor> faceGradDDQuad = gradDDQuad[faceI];
        const UList<tensor> faceFQuadOld = FQuadOld[faceI];
        UList<tensor> faceGradDTotalQuad = gradDTotalQuad[faceI];

        forAll(faceGradDTotalQuad, qpI)
        {
            // Relative deformation gradient
            const tensor relF(I + faceGradDDQuad[qpI].T());

            // Total deformation gradient
            const tensor F(relF & faceFQuadOld[qpI]);

            // The mechanical law reconstructs the total deformation gradient
            // as F = I + gradD.T(), so we pass the equivalent total
            // displacement gradient
            faceGradDTotalQuad[qpI] = (F - I).T();
        }
    }
}


void nonLinGeomUpdatedLagSolid::storeQuadratureKinematicsOldTime() const
{
    if (FQuadOldPtr_.empty())
    {
        makeQuadratureKinematics();
    }

    const CompactListList<tensor>& gradDTotalQuad = gradDTotalQuadPtr_();
    CompactListList<tensor>& FQuadOld = FQuadOldPtr_();

    forAll(FQuadOld, faceI)
    {
        const UList<tensor> faceGradDTotalQuad = gradDTotalQuad[faceI];
        UList<tensor> faceFQuadOld = FQuadOld[faceI];

        forAll(faceFQuadOld, qpI)
        {
            faceFQuadOld[qpI] = I + faceGradDTotalQuad[qpI].T();
        }
    }
}
#endif


void nonLinGeomUpdatedLagSolid::enforceTractionBoundaries
(
    surfaceVectorField& force,
    const volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
        vectorField& forceP = boundaryFieldRef(force)[patchI];

        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // Specified traction at the patch faces
            vectorField tracP(nPatch.size(), vector::zero);

            if (highOrderResidual())
            {
#ifndef FOAMEXTEND
                // Face quadrature points weights
                const CompactListList<scalar>& faceQuadWeights =
                    displacementMLS().quadrature().faceQuadWeights();

                const surfaceScalarField& magSf = mesh().magSf();

                // Get value at patch faces quadrature points
                autoPtr<CompactListList<vector>> patchQuadraturePointsValue =
                    tracPatch.evaluateQuadrature();

                const CompactListList<vector>& quadratureValues =
                    patchQuadraturePointsValue();

                forAll(mesh().boundaryMesh()[patchI], faceI)
                {
                    const label start = mesh().boundaryMesh()[patchI].start();

                    // Get global face index
                    const label faceID = faceI + start;

                    // Get the number of quadrature points for this face
                    const label nPoints = faceQuadWeights[faceID].size();

                    // Loop over quadrature points and add their contribution
                    for (label pointI = 0; pointI < nPoints; ++pointI)
                    {
                        tracP[faceI] +=
                            quadratureValues[faceI][pointI]
                           *faceQuadWeights[faceID][pointI];
                    }
                    // Divide with area because we use physical weights
                    tracP[faceI] *=
                        (1.0/(magSf.boundaryField()[patchI][faceI]));
                }
#endif
            }
            else
            {
                tracP = tracPatch.traction() - nPatch*tracPatch.pressure();
            }

            if (tracPatch.useUndeformedArea())
            {
                notImplemented("Not implemented for updated Lagrangian");

                // const scalarField& magSfPatch =
                //     D.mesh().boundary()[patchI].magSf();

                // forceP = tracP*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                forceP = tracP*magSfCurrentPatch;
            }
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<slipFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // Set shear traction to zero
            // traction.boundaryFieldRef()[patchI] =
                // sqr(nPatch) & traction.boundaryField()[patchI];
            forceP = sqr(nPatch) & force.boundaryField()[patchI];
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

    // Updated (end of last time step) unit normal vectors at the faces
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Calculate deformed area vectors and normals
        const surfaceVectorField SfCurrent
        (
            fvc::interpolate(relJ_*relFinv_.T()) & mesh().Sf()
        );
        const surfaceScalarField magSfCurrent(mag(SfCurrent));
        const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

        // Traction vectors at the faces
        surfaceVectorField traction(nCurrent & fvc::interpolate(sigma()));

        // Add stabilisation to the traction
        // We add this before enforcing the traction condition as the stabilisation
        // is set to zero on traction boundaries
        momentumStabilisation().updateVector(DD(), &gradDD());
        traction += impKf_*momentumStabilisation().faceVector();

        // Calculate the force at the faces
        surfaceVectorField force(magSfCurrent*traction);

        // Enforce traction boundary conditions
        enforceTractionBoundaries(force, DD(), nCurrent, magSfCurrent);

        // Momentum equation incremental updated Lagrangian form
        // Assemble the RHS in stages.
        tmp<fvVectorMatrix> tRhsEqn
        (
            fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
        );
        tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)");
        tmpRef(tRhsEqn) += fvc::div(force);
        tmpRef(tRhsEqn) += rho_*g();

        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho_, DD())
          + fvc::d2dt2(rho_, D().oldTime())
         == tRhsEqn
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

    // Total displacement
    D() = D().oldTime() + DD();

    // Update gradient of total displacement
    // Do we need this?
    gradD() = fvc::grad(D());

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
#ifdef USE_PETSC
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
        // Map the D field to the SNES solution vector
        foamPetscSnesHelper::InsertFieldComponents<vector>
        (
            primitiveFieldRef(DD()),
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD()
          ? makeList<label>({0,1})
          : makeList<label>({0,1,2})
        );
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the DD field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        primitiveFieldRef(DD()),
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    DD().correctBoundaryConditions();

    // Total displacement
    D() = D().oldTime() + DD();

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Update the kinematic fields using the high-order gradient
        gradDD() = displacementMLS().grad(DD());
        relF_ = I + gradDD().T();
        relFinv_ = inv(relF_);
        relJ_ = det(relF_);
        F_ = relF_ & F_.oldTime();
        J_ = relJ_*J_.oldTime();

        // Update the total deformation gradient at the face quadrature points
        // so that it is consistent with the converged solution
        mechanical().grad(DD(), gradDQuad());
        updateQuadratureKinematics();

        // Calculate the cell centre stress using run-time selectable
        // mechanical law
        mechanical().correct(sigma());
#endif
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD(), gradDD(), pointDD());
    pointDD().correctBoundaryConditions();

    // Total point displacement
    pointD() = pointD().oldTime() + pointDD();;

    // Velocity
    U() = fvc::ddt(D());

    // Acceleration
    A_ = fvc::d2dt2(D());

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

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
        "DD",
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
    rKappaPtr_(),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    )
{
    DDisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(rho_, DD());
    fvc::d2dt2(rho_, D().oldTime());

    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        // It is important to call the stress calculation procedure during the
        // constructor to allow it to correctly initialise fields
        mechanical().correct(sigma());
    }

    if (solvePressure())
    {
        if (solutionAlg() != solutionAlgorithm::PETSC_SNES)
        {
            FatalErrorInFunction
                << "The solution algorithm must be "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " when solvePressure is enabled" << abort(FatalError);
        }

        // Ensure p is created and the oldTime is stored
        p().oldTime();
    }

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
#ifdef OPENFOAM_NOT_EXTEND
        mesh().gradScheme("grad(" + DD().name() +')')
#else
        mesh().schemesDict().gradScheme("grad(" + DD().name() +')')
#endif
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
                        boundaryFieldRef(DD())[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
}


void nonLinGeomUpdatedLagSolid::makeRKappa() const
{
    if (rKappaPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rKappaPtr_.set(new volScalarField(1.0/mechanical().bulkModulus()));
}


const volScalarField& nonLinGeomUpdatedLagSolid::rKappa() const
{
    if (rKappaPtr_.empty())
    {
        makeRKappa();
    }

    return rKappaPtr_();
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


#ifdef USE_PETSC

label nonLinGeomUpdatedLagSolid::initialiseJacobian(Mat& jac)
{
#ifndef FOAMEXTEND
    if (highOrderJacobian())
    {
        // Note: the non-zero structure is preallocated from the moving least
        // squares stencils on the initial mesh. As this solid model moves the
        // mesh, the stencils are re-calculated each time step (see
        // updateTotalFields) and, if they change, PETSc will report a new
        // non-zero allocation error. This has been verified for rigid motion,
        // where the stencils are unchanged by construction; a case with
        // sufficient mesh deformation to re-order a stencil would require the
        // Jacobian to be re-created after mesh motion
        return hofvm::initialiseJacobian
        (
            jac,
            *this,
            displacementMLS(),
            DD(),
            blockSize_
        );
    }
#endif

    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label nonLinGeomUpdatedLagSolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label nonLinGeomUpdatedLagSolid::formResidual
(
    Vec f,
    const Vec x
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
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    DD.correctBoundaryConditions();

    // Update total displacement
    D() = D().oldTime() + DD;

    if (solvePressure() && highOrderResidual())
    {
        FatalErrorInFunction
            << "solvePressure must be disabled when using high order "
            << "residual calculation. Mixed approach not yet supported!"
            << abort(FatalError);
    }

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Update cell-centre displacement increment gradient
        gradDD() = displacementMLS().grad(DD);

        // Update displacement increment gradient at the face quadrature points
        mechanical().grad(DD, gradDQuad());
#endif
    }
    else
    {
        // Update displacement increment gradient
        mechanical().grad(DD, gradDD());
    }

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

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Accumulate the total deformation gradient at the face quadrature
        // points, as the mechanical law expects the total displacement gradient
        updateQuadratureKinematics();

        // Calculate sigma at the face quadrature points
        mechanical().correct(gradDTotalQuadPtr_(), sigmaQuad());

        // Calculate the cell-centre stress using run-time selectable
        // mechanical law
        // The residual uses the quadrature point stress, but the cell-centre
        // stress is still required by tractionBoundarySnGrad
        mechanical().correct(sigma());
#endif
    }
    else
    {
        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        sigma() = dev(sigma()) - p*I;

        // Calculate the pressure gradient
        const volVectorField gradp(fvc::grad(p));

        // Re-calculate the pressure stabilisation parameter
        pressureStabilisation().updateScalar(p, &gradp);

        // Dimensional consistency factor
        const dimensionedScalar one
        (
            "one", dimensionSet(-2, 4, 4, 0, 0, 0, 0), 1.0
        );

        // Compute the positive face-interpolated reciprocal of the approximate
        // momentum equation diagonal. This is the solid analogue of rAUf in
        // pressure-velocity coupling and has units of [Pa].
        {
            fvVectorMatrix approxMomJ
            (
                fvm::laplacian(impKf_, DD, "laplacian(DDD,DD)")
              - rho()*fvm::d2dt2(DD)
            );
            approxMomJ.relax();
            rAUf() = -1.0/(fvc::interpolate(approxMomJ.A())*one);
        }

        // Calculate pressure equation residual
        scalarField pressureResidual
        (
          - p*rKappa()
          + pressureStabilisation().cellScalar(&rAUf(), true)*one
          - 0.5*(pow(J_, 2.0) - 1.0)/J_
        );

        // Make residual extensive
        pressureResidual *= mesh().V();

        // Copy the pressureResidual into the f field as the final equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1
        );
    }

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

#ifndef FOAMEXTEND
    if (highOrderResidual())
    {
        // Replace the Cauchy traction (force per unit deformed area) with the
        // relative nominal traction (force per unit area in the configuration
        // at the start of the time step, i.e. the current mesh), calculated by
        // integrating the first Piola-Kirchhoff stress over the face
        // quadrature points
        // Note: gradDQuad() holds the displacement increment gradient, so the
        // deformation gradient assembled within surfaceIntegrate is the
        // relative deformation gradient, as required here
        traction = hofvc::surfaceIntegrate(sigmaQuad(), gradDQuad(), mesh());
    }
#endif

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // Note: in the high-order case, the stabilisation is referred to the
    // configuration at the start of the time step; it vanishes on convergence,
    // so this does not change the converged answer
    momentumStabilisation().updateVector(DD, &gradDD());
    traction += impKf_*momentumStabilisation().faceVector();

    // Calculate the force at the faces
    surfaceVectorField force
    (
        highOrderResidual()
      ? mesh().magSf()*traction
      : magSfCurrent*traction
    );

    // Enforce traction boundary conditions
    enforceTractionBoundaries(force, DD, nCurrent, magSfCurrent);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual
    (
        fvc::div(force)
      + rho()
       *(
           g() - dampingCoeff()*fvc::ddt(D())
        )
    );

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        residual -= rho()*hofvc::d2dt2(D());
#endif
    }
    else
    {
        residual -= rho()*fvc::d2dt2(D());
    }

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
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    return 0;
}


label nonLinGeomUpdatedLagSolid::formJacobian
(
    Mat jac,
    const Vec x
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
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    DD.correctBoundaryConditions();

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();

        {
            // Dimensional consistency factor
            const dimensionedScalar one
            (
                "one", dimensionSet(-2, 4, 4, 0, 0, 0, 0), 1.0
            );

            // Compute the positive face-interpolated reciprocal of the approximate
            // momentum equation diagonal (solid analogue of rAUf), [Pa]
            {
                fvVectorMatrix approxMomJ
                (
                    fvm::laplacian(impKf_, DD, "laplacian(DDD,DD)")
                  - rho()*fvm::d2dt2(DD)
                );
                approxMomJ.relax();
                rAUf() = -1.0/(fvc::interpolate(approxMomJ.A())*one);
            }

            fvScalarMatrix approxPressureJ
            (
              - fvm::Sp(rKappa(), p)
              + one*pressureStabilisation().scalarJacobian(p, &rAUf())
            );

            // Insert the pressure equation
            foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            (
                approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            );
        }

        // Insert DD-in-p equation coeffs coming from tr(grad(DD)) == div(DD)
        foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
        (
            p,
            DD,
            jac,
            blockSize_ - 1,            // row offset
            0,                         // column offset
            solidModel::twoD() ? 2 : 3 // number of scalar components of DD
        );

        // Insert p-in-DD term
        foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
        (
            p,
            jac,
            0,                         // row offset
            blockSize_ - 1,            // column offset
            solidModel::twoD() ? 2 : 3 // number of scalar equations to insert
        );
    }

    if (highOrderJacobian())
    {
#ifndef FOAMEXTEND
        // Note: unlike the fallback fvVectorMatrix approxJ path below, we do
        // not currently apply matrix under-relaxation to the high-order
        // Jacobian assembled directly into PETSc. If this becomes important
        // for robustness, an equivalent relaxation step may need to be added.
        tmp<volScalarField> tK = mechanical().bulkModulus();
        const volScalarField& K = tK();

        tmp<volScalarField> tMu = (impK_ - K)*(3.0/4.0);
        const volScalarField& mu = tMu();

        tmp<volScalarField> tLambda = impK_ - 2.0*mu;
        const volScalarField& lambda = tLambda();

        const movingLeastSquares& mls = displacementMLS();

        hofvm::laplacianIntoPETScMatrix
        (
            jac,
            *this,
            mls,
            DD,
            mu
        );

        hofvm::laplacianTransposeIntoPETScMatrix
        (
            jac,
            *this,
            mls,
            DD,
            mu
        );

        hofvm::laplacianTraceIntoPETScMatrix
        (
            jac,
            *this,
            mls,
            DD,
            lambda
        );

        fvVectorMatrix transientJ
        (
          - rho()*hofvm::d2dt2(DD)
        );

        if (dampingCoeff().value() > SMALL)
        {
            transientJ -= dampingCoeff()*rho()*fvmDdtVectorCompat(DD);
        }

        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
        (
            transientJ, jac, 0, 0, solidModel::twoD() ? 2 : 3
        );
#endif
    }
    else
    {
        // Calculate a segregated approximation of the Jacobian
        fvVectorMatrix approxJ
        (
            momentumStabilisation().vectorJacobian(DD, &impKf_)
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
    }

    return 0;
}


#endif // USE_PETSC

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

#ifndef FOAMEXTEND
    if (highOrderResidual())
    {
        // Store the total deformation gradient at the face quadrature points:
        // it is the reference for the accumulation in the next time step
        storeQuadratureKinematicsOldTime();
    }
#endif

    // Move the mesh to the deformed configuration
#ifdef OPENFOAM_NOT_EXTEND
    const vectorField oldPoints = mesh().points();
#else
    const vectorField oldPoints = mesh().allPoints();
#endif
    moveMesh(oldPoints, pointDD());

#ifndef FOAMEXTEND
    if (highOrderResidual() || highOrderJacobian())
    {
        // The moving least squares stencils, quadrature points and
        // interpolation coefficients are all geometric, so they are
        // re-calculated here on the moved mesh
        //
        // This is the first solid model with a moving mesh to use the
        // high-order approach, so it is the first to face the choice between:
        //   1. building the stencils once on the initial mesh and re-using
        //      them, which is cheaper and keeps the Jacobian sparsity pattern
        //      fixed, but lets the stencils become progressively less
        //      appropriate as the mesh deforms; and
        //   2. re-calculating them after every mesh motion, as done here,
        //      which keeps the stencils consistent with the current mesh at
        //      the cost of rebuilding them each time step.
        //
        // Option 2 is used here because it is the conservative choice: the
        // stencil always matches the mesh it is used on. Note that for a rigid
        // motion the two options are equivalent, since the stencil is
        // distance-sorted and a rigid motion preserves all distances.
        //
        // The distinction is expected to matter more when this approach is
        // extended to fluids, where the mesh may deform far more than in a
        // typical solid case
        //
        // Note: FQuadOldPtr_ and gradDTotalQuadPtr_ are deliberately not
        // cleared, as the quadrature points move with the material
        clearMovingLeastSquaresData();
    }
#endif

    solidModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
