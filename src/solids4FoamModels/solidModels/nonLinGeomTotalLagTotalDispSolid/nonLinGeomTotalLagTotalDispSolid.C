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

#include "nonLinGeomTotalLagTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
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

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagTotalDispSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void nonLinGeomTotalLagTotalDispSolid::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    // D() = D().oldTime() + U()*runTime().deltaT();
    D() = D().oldTime() + U()*runTime().deltaT()
        + 0.5*sqr(runTime().deltaT())*A_;

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Predict p using the dp/dt field
        p() = p().oldTime() + autoPtrRef(dpdtPtr_)*runTime().deltaT();
        // p() = p().oldTime() + dpdt*runTime().deltaT()
        //     + 0.5*sqr(runTime().deltaT())*d2pdt2;

        sigma() = dev(sigma()) - p()*I;
    }
}


void nonLinGeomTotalLagTotalDispSolid::enforceTractionBoundaries
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
                const scalarField& magSfPatch =
                    D.mesh().boundary()[patchI].magSf();

                forceP = tracP*magSfPatch;
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


bool nonLinGeomTotalLagTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver using an implicit segregated approach"
        << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    if (predictor_ && newTimeStep())
    {
        predict();
    }

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    int iCorr = 0;
    scalar currentResidualNorm = 0;
    scalar initialResidualNorm = 0;
    scalar deltaXNorm = 0;
    scalar xNorm = 0;
    const convergenceParameters convParam =
        readConvergenceParameters(solidModelDict());

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    // Undeformed unit normal vectors at the faces
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Calculate deformed area vectors and normals
        const surfaceVectorField SfCurrent
        (
            fvc::interpolate(J_*Finv_.T()) & mesh().Sf()
        );
        const surfaceScalarField magSfCurrent(mag(SfCurrent));
        const surfaceVectorField nCurrent(SfCurrent/magSfCurrent);

        // Traction vectors at the faces
        surfaceVectorField traction(nCurrent & fvc::interpolate(sigma()));

        // Add stabilisation to the traction
        // We add this before enforcing the traction condition as the stabilisation
        // is set to zero on traction boundaries
        momentumStabilisation().updateVector(D(), &gradD());
        traction += impKf_*momentumStabilisation().faceVector();

        // Calculate the force at the faces
        surfaceVectorField force(magSfCurrent*traction);

        // Enforce traction boundary conditions
        enforceTractionBoundaries(force, D(), nCurrent, magSfCurrent);

        // Momentum equation total displacement total Lagrangian form
#ifndef OPENFOAM_COM
        // Assemble the RHS in stages.
        // The equivalent chained tmp fvMatrix expression is stable on OpenFOAM.com.
        tmp<fvVectorMatrix> tRhsEqn
        (
            fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
        );
        tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, D(), "laplacian(DD,D)");
        tmpRef(tRhsEqn) += fvc::div(force);
        tmpRef(tRhsEqn) += rho()*g();

        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == tRhsEqn
        );
#else
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(force)
          + rho()*g()
        );
#endif

        // Add damping
        if (dampingCoeff().value() > SMALL)
        {
            DEqn += dampingCoeff()*rho()*fvm::ddt(D());
        }

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
#ifdef OPENFOAM_NOT_EXTEND
                        D().primitiveField() - D().prevIter().primitiveField()
#else
                        D().internalField() - D().prevIter().internalField()
#endif
                    )
                )
            );

        // Norm of the solution
#ifdef OPENFOAM_NOT_EXTEND
        xNorm = sqrt(gSum(magSqr(D().primitiveField())));
#else
        xNorm = sqrt(gSum(magSqr(D().internalField())));
#endif

        // Store the initial residual
        if (iCorr == 0)
        {
            initialResidualNorm = currentResidualNorm;
            Info<< "Initial Residual Norm = " << initialResidualNorm << nl
                << "Initial Solution Norm = " << xNorm << endl;
        }

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
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

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

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


bool nonLinGeomTotalLagTotalDispSolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Use the segregated solver as a predictor
        //evolveImplicitSegregated();

        // Map the D field to the SNES solution vector
        foamPetscSnesHelper::InsertFieldComponents<vector>
        (
            primitiveFieldRef(D()),
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD()
          ? makeList<label>({0,1})
          : makeList<label>({0,1,2})
        );

        if (solvePressure())
        {
            // Map the p field to the SNES solution vector
            foamPetscSnesHelper::InsertFieldComponents<scalar>
            (
                p(),
                foamPetscSnesHelper::solution(),
                blockSize_ -1 // Location of first component
            );
        }
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
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

    if (solvePressure())
    {
        // Map the PETSc solution to the p field
        // p is located in the last ("blockSize - 1") component
        scalarField& pI = p();
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            foamPetscSnesHelper::solution(),
            pI,
            blockSize_ - 1 // Location of p component
        );

        p().correctBoundaryConditions();

        // Update dpdt
        autoPtrRef(dpdtPtr_) = fvc::ddt(p());
    }

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Update the kinematic fields using the high-order gradient
        gradD() = displacementMLS().grad(D());
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);

        // Calculate the cell centre stress using run-time selectable
        // mechanical law
        mechanical().correct(sigma());
#endif
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

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

nonLinGeomTotalLagTotalDispSolid::nonLinGeomTotalLagTotalDispSolid
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
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(F_)
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
        fvc::d2dt2(D())
    ),
    impK_
    (
        solvePressure()
      ? 2.0*mechanical().shearModulus()
      : mechanical().impK()
    ),
    impKf_(fvc::interpolate(impK_)),
    rImpK_(1.0/impK_),
    rKappaPtr_(),
    dpdtPtr_(),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    )
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // It is important to call the stress calculation procedure during the
    // constructor to allow it to correctly initialise fields
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        mechanical().correct(sigma());
    }

    Info<< "solvePressure = " << solvePressure() << endl;
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

        // Initialise dpdt field
        dpdtPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                fvc::ddt(p())
            )
        );
    }

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        DD() = D() - D().oldTime();
        mechanical().grad(D(), gradD());
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }

    // Check the gradScheme
    const word gradDScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh().gradScheme("grad(" + D().name() +')')
#else
        mesh().schemesDict().gradScheme("grad(" + D().name() +')')
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
        forAll(D().boundaryField(), patchI)
        {
            if
            (
                isA<solidTractionFvPatchVectorField>
                (
                    D().boundaryField()[patchI]
                )
            )
            {
                Info<< "    Setting `extrapolateValue` to `true` on the "
                    << mesh().boundary()[patchI].name() << " patch of the D "
                    << "field" << endl;

                solidTractionFvPatchVectorField& tracPatch =
                    refCast<solidTractionFvPatchVectorField>
                    (
                        boundaryFieldRef(D())[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
}


void nonLinGeomTotalLagTotalDispSolid::makeRKappa() const
{
    if (rKappaPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    rKappaPtr_.set(new volScalarField(1.0/mechanical().bulkModulus()));
}


const volScalarField& nonLinGeomTotalLagTotalDispSolid::rKappa() const
{
    if (rKappaPtr_.empty())
    {
        makeRKappa();
    }

    return rKappaPtr_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagTotalDispSolid::evolve()
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

label nonLinGeomTotalLagTotalDispSolid::initialiseJacobian(Mat& jac)
{
#ifndef FOAMEXTEND
    if (highOrderJacobian())
    {
        return hofvm::initialiseJacobian
        (
            jac,
            *this,
            displacementMLS(),
            D(),
            blockSize_
        );
    }
#endif

    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label nonLinGeomTotalLagTotalDispSolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label nonLinGeomTotalLagTotalDispSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();

    // Copy x into the D field
    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
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
    D.correctBoundaryConditions();

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
        // Update cell-centre gradient of displacement
        gradD() = displacementMLS().grad(D);

        // Update gradient of displacement at face quadrature points
        mechanical().grad(D, gradDQuad());
#endif
    }
    else
    {
        // Update gradient of displacement
        mechanical().grad(D, gradD());

        // Enforce the boundary conditions again for any conditions that use
        // gradD
        D.correctBoundaryConditions();
    }

    // Increment of displacement
    DD() = D - D.oldTime();

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Calculate sigma at the face quadrature points
        mechanical().correct(gradDQuad(), sigmaQuad());
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
                fvm::laplacian(impKf_, D, "laplacian(DD,D)")
              - rho()*fvm::d2dt2(D)
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
        pressureResidual *= mesh.V();

        // Copy the pressureResidual into the f field as the final equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1
        );
    }

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfCurrent
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
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
        // nominal traction (force per unit undeformed area), calculated by
        // integrating the first Piola-Kirchhoff stress over the face
        // quadrature points
        traction = hofvc::surfaceIntegrate(sigmaQuad(), gradDQuad(), mesh);
    }
#endif

    //fvc::div(J_*Finv_ & sigma(), "div(sigma)");

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // Note: in the high-order case, the stabilisation is referred to the
    // undeformed configuration; it vanishes on convergence, so this does not
    // change the converged answer
    momentumStabilisation().updateVector(D, &gradD());
    traction += impKf_*momentumStabilisation().faceVector();

    // Calculate the force at the faces
    surfaceVectorField force
    (
        highOrderResidual()
      ? mesh.magSf()*traction
      : magSfCurrent*traction
    );

    // Enforce traction boundary conditions
    enforceTractionBoundaries(force, D, nCurrent, magSfCurrent);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual
    (
        // fvc::div(magSfCurrent*traction)
        fvc::div(force)
      + rho()
       *(
            g() - dampingCoeff()*fvc::ddt(D)
        )
    );

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        residual -= rho()*hofvc::d2dt2(D);
#endif
    }
    else
    {
        residual -= rho()*fvc::d2dt2(D);
    }

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    //residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();

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


label nonLinGeomTotalLagTotalDispSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Copy x into the D field
    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
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
    D.correctBoundaryConditions();

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
                    fvm::laplacian(impKf_, D, "laplacian(DD,D)")
                  - rho()*fvm::d2dt2(D)
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

            // Insert D-in-p equation coeffs coming from tr(grad(D)) == div(D)
            foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
            (
                p,
                D,
                jac,
                blockSize_ - 1,            // row offset
                0,                         // column offset
                solidModel::twoD() ? 2 : 3 // number of scalar components of D
            );

            // Insert p-in-D term
            // Insert "-grad(p)" (equivalent to "-div(p*I)") into the D equation
            foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
            (
                p,
                jac,
                0,                         // row offset
                blockSize_ - 1,            // column offset
                solidModel::twoD() ? 2 : 3 // number of scalar equations to insert
            );
        }
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
            D,
            mu
        );

        hofvm::laplacianTransposeIntoPETScMatrix
        (
            jac,
            *this,
            mls,
            D,
            mu
        );

        hofvm::laplacianTraceIntoPETScMatrix
        (
            jac,
            *this,
            mls,
            D,
            lambda
        );

        fvVectorMatrix transientJ
        (
          - rho()*hofvm::d2dt2(D)
        );

        if (dampingCoeff().value() > SMALL)
        {
            transientJ -= dampingCoeff()*rho()*fvmDdtVectorCompat(D);
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
            momentumStabilisation().vectorJacobian(D, &impKf_)
          - rho()*fvm::d2dt2(D)
        );

        if (dampingCoeff().value() > SMALL)
        {
            approxJ -= dampingCoeff()*rho()*fvm::ddt(D);
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

tmp<vectorField> nonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
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
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
