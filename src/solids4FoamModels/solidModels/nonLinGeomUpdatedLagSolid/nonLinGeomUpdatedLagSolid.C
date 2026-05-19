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

            // traction.boundaryFieldRef()[patchI] =
            //     tracPatch.traction() - nPatch*tracPatch.pressure();
            if (tracPatch.useUndeformedArea())
            {
                notImplemented("Not implemented for updated Lagrangian");

                // const scalarField& magSfPatch =
                //     D.mesh().boundary()[patchI].magSf();

                // forceP =
                //     (
                //         tracPatch.traction() - nPatch*tracPatch.pressure()
                //     )*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfCurrentPatch;
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
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for DD using PETSc SNES" << endl;

    // Update D boundary conditions
    DD().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Seed the PETSc solution vector from the predicted fields
        packSolution(foamPetscSnesHelper::solution());
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Map the PETSc solution back into the DD field (and p when active)
    unpackSolution(foamPetscSnesHelper::solution());

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
    rAUfTimeIndex_(-1),
    rAUfDeltaT_(0),
    scaleMixedPetScFields_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "scaleMixedPetScFields", true
        )
    ),
    pressureUnknownScaleType_
    (
        solidModelDict().lookupOrDefault<word>("pressureUnknownScale", "twoMu")
    ),
    pressureUnknownScale_(1.0),
    pressureScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("pressureScaleFactor", 1.0)
    ),
    pressureScaleByTwoMu_
    (
        solidModelDict().lookupOrDefault<Switch>("pressureScaleByTwoMu", true)
    ),
    twoMuRef_(1.0),
    pressureEqnScale_(pressureScaleFactor_),
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

        // Use a volume-weighted average of 2*mu as the physical scale
        // of the pressure equation. The pressure-row residual and
        // Jacobian are then multiplied by
        // pressureEqnScale_ = pressureScaleFactor_ * twoMuRef_ so that
        // their natural magnitude is comparable to the momentum block.
        const volScalarField twoMu(2.0*mechanical().shearModulus());
        scalar twoMuV = 0;
        scalar Vtot = 0;
        forAll(twoMu, cellI)
        {
            const scalar Vc = mesh().V()[cellI];
            twoMuV += twoMu[cellI]*Vc;
            Vtot += Vc;
        }
        reduce(twoMuV, sumOp<scalar>());
        reduce(Vtot, sumOp<scalar>());
        twoMuRef_ = twoMuV/Vtot;
        pressureEqnScale_ =
            pressureScaleFactor_*(pressureScaleByTwoMu_ ? twoMuRef_ : 1.0);

        pressureUnknownScale_ = 1.0;
        if (scaleMixedPetScFields_)
        {
            if
            (
                pressureUnknownScaleType_ == "twoMu"
             || pressureUnknownScaleType_ == "2mu"
            )
            {
                pressureUnknownScale_ = twoMuRef_;
            }
            else if
            (
                pressureUnknownScaleType_ == "user"
             || pressureUnknownScaleType_ == "scalar"
            )
            {
                pressureUnknownScale_ =
                    readScalar
                    (
                        solidModelDict().lookup("pressureUnknownScaleValue")
                    );
            }
            else if (pressureUnknownScaleType_ == "none")
            {
                pressureUnknownScale_ = 1.0;
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown pressureUnknownScale "
                    << pressureUnknownScaleType_ << nl
                    << "Valid options are twoMu, user, scalar, none"
                    << abort(FatalError);
            }

            if (pressureUnknownScale_ <= VSMALL)
            {
                FatalErrorInFunction
                    << "pressureUnknownScale must be positive, found "
                    << pressureUnknownScale_ << abort(FatalError);
            }
        }

        Info<< "pressureEqnScale = " << pressureEqnScale_
            << ", where pressureScaleFactor = " << pressureScaleFactor_
            << " and 2*mu = " << twoMuRef_ << endl;

        Info<< "PETSc pressure unknown scale = " << pressureUnknownScale_
            << " (scaleMixedPetScFields = " << scaleMixedPetScFields_
            << ", pressureUnknownScale = " << pressureUnknownScaleType_
            << ")" << endl;
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


void nonLinGeomUpdatedLagSolid::updateRAUfIfStale()
{
    if (!solvePressure())
    {
        // rAUf is unused unless we are solving the mixed system
        return;
    }

    const label tIdx = runTime().timeIndex();
    const scalar dt = runTime().deltaT().value();

    // rAUf depends on (impKf_, rho, mesh, deltaT) only. It is fresh
    // when both the timeIndex and deltaT match what we cached
    if
    (
        rAUfTimeIndex_ >= 0
     && rAUfTimeIndex_ == tIdx
     && mag(dt - rAUfDeltaT_) <= SMALL*max(mag(dt), SMALL)
    )
    {
        return;
    }

    // Build the approximate momentum diagonal. fvm::laplacian and
    // fvm::d2dt2 read only the BC structure of DD, not its values, so
    // the resulting diagonal is independent of the current Newton
    // iterate and any MFFD perturbation
    fvVectorMatrix approxMomJ
    (
        fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
      - rho()*fvm::d2dt2(DD())
    );
    approxMomJ.relax();
    rAUf() = -1.0/(fvc::interpolate(approxMomJ.A()));

    rAUfTimeIndex_ = tIdx;
    rAUfDeltaT_ = dt;
}


#ifdef USE_PETSC

void nonLinGeomUpdatedLagSolid::unpackSolution(const Vec x)
{
    // Copy x into the DD field
    volVectorField& DD = const_cast<volVectorField&>(this->DD());
    vectorField& DDI = DD;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DDI,
        0, // Location of first DDI component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions on DD
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

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Copy the scaled pressure unknown pHat from x into the
        // physical p field via p = pressureUnknownScale_ * pHat
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        scalarField pHat(pI.size(), 0.0);
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pHat, blockSize_ - 1
        );
        pI = pressureUnknownScale_*pHat;

        // Enforce the boundary conditions on p
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        sigma() = dev(sigma()) - p*I;
    }
}


void nonLinGeomUpdatedLagSolid::packSolution(Vec x)
{
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        primitiveField(DD()),
        x,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    if (solvePressure())
    {
        // Insert the scaled pressure unknown pHat = p/pressureUnknownScale_
        scalarField pHat(primitiveField(p()));
        pHat /= pressureUnknownScale_;
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pHat,
            x,
            blockSize_ - 1
        );
    }
}

#endif // USE_PETSC


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
    // Copy x into the DD field (and p when solvePressure() is active),
    // refresh dependent kinematic fields and correct boundary
    // conditions
    unpackSolution(x);

    // Take a non-const reference to DD for local use below
    volVectorField& DD = const_cast<volVectorField&>(this->DD());

    if (solvePressure())
    {
        // Pressure has already been unpacked from x by unpackSolution(x)
        // above, its BCs corrected and the pressure component of sigma
        // updated; take a reference for local use
        volScalarField& p = const_cast<volScalarField&>(this->p());

        // Calculate the pressure gradient
        const volVectorField gradp(fvc::grad(p));

        // Re-calculate the pressure stabilisation parameter
        pressureStabilisation().updateScalar(p, &gradp);

        // Refresh rAUf (the positive face-interpolated reciprocal of
        // the approximate momentum equation diagonal -- the solid
        // analogue of rAUf in pressure-velocity coupling, units [m^2/Pa])
        // only when the mesh or deltaT have changed. The diagonal is
        // value-independent of DD and p, so this is safe under PETSc
        // matrix-free finite-difference perturbations
        updateRAUfIfStale();

        // Calculate pressure equation residual. Keep the finite-strain
        // volumetric term -0.5*(J^2-1)/J as in the published nonlinear
        // updated-Lagrangian formulation; its linearisation about DD=0
        // yields -V*div(DD) which is what InsertFvmDivUIntoPETScMatrix
        // assembles in formJacobian.
        scalarField pressureResidual
        (
          - p*rKappa()
          + pressureStabilisation().cellScalar(&rAUf(), true)
          - 0.5*(pow(J_, 2.0) - 1.0)/J_
        );

        // Make residual extensive
        pressureResidual *= mesh().V();

        // Apply the physical row scaling. pressureEqnScale_ already
        // bakes in both the user-facing pressureScaleFactor and the
        // 2*mu physical scale.
        if (pressureEqnScale_ != 1.0)
        {
            pressureResidual *= pressureEqnScale_;
        }

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

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    momentumStabilisation().updateVector(DD, &gradDD());
    traction += impKf_*momentumStabilisation().faceVector();

    // Calculate the force at the faces
    surfaceVectorField force(magSfCurrent*traction);

    // Enforce traction boundary conditions
    enforceTractionBoundaries(force, DD, nCurrent, magSfCurrent);

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
    // Copy x into the DD field (and p when solvePressure() is active),
    // refresh dependent kinematic fields and correct boundary
    // conditions
    unpackSolution(x);

    // Take a non-const reference to DD for local use below
    volVectorField& DD = const_cast<volVectorField&>(this->DD());

    if (solvePressure())
    {
        // Pressure has already been unpacked from x by unpackSolution(x)
        // above and its BCs corrected; take a reference for local use
        volScalarField& p = const_cast<volScalarField&>(this->p());

        {
            // Refresh rAUf only when the mesh or deltaT have changed
            // (the diagonal of the approximate momentum Jacobian is
            // independent of DD and p values)
            updateRAUfIfStale();

            fvScalarMatrix approxPressureJ
            (
              - pressureEqnScale_*pressureUnknownScale_*fvm::Sp(rKappa(), p)
              + pressureEqnScale_*pressureUnknownScale_
               *pressureStabilisation().scalarJacobian(p, &rAUf())
            );

            // Insert the pressure equation
            foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            (
                approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            );
        }

        // Insert DD-in-p equation coefficients matching the
        // linearisation of -0.5*(J^2-1)/J about DD=0, which to leading
        // order on the current mesh equals -V*div(DD).
        //
        // InsertFvmDivU's sign convention: scale=+1 assembles
        // `-V*div(U)`. So we pass `+pressureEqnScale_` to get
        // J_pDD = -alpha*V*div(DD).
        foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
        (
            p,
            DD,
            jac,
            blockSize_ - 1,             // row offset (p row)
            0,                          // column offset (DD columns)
            solidModel::twoD() ? 2 : 3, // number of DD components
            pressureEqnScale_           // scale (helper returns -V*div with +1)
        );

        // Insert p-in-DD term. InsertFvmGrad's updated sign
        // convention: scale=+1 assembles `-V*grad(p)` (matching the
        // momentum equation contribution of -grad(p)). Apply the
        // pressure-unknown scale so the column corresponds to the
        // scaled unknown pHat = p/pressureUnknownScale_.
        foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
        (
            p,
            jac,
            0,                          // row offset
            blockSize_ - 1,             // column offset
            solidModel::twoD() ? 2 : 3, // number of DD components
            pressureUnknownScale_       // scale (helper returns -V*grad with +1)
        );
    }

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
