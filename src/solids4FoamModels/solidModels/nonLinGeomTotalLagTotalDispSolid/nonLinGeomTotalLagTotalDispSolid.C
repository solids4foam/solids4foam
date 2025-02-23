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
        p() = p().oldTime() + dpdtPtr_.ref()*runTime().deltaT();
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
                const scalarField& magSfPatch =
                    D.mesh().boundary()[patchI].magSf();

                force.boundaryFieldRef()[patchI] =
                (
                    tracPatch.traction() - nPatch*tracPatch.pressure()
                )*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                force.boundaryFieldRef()[patchI] =
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
        )
        {
            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // Set shear traction to zero
            // traction.boundaryFieldRef()[patchI] =
                // sqr(nPatch) & traction.boundaryField()[patchI];
            force.boundaryFieldRef()[patchI] =
                sqr(nPatch) & force.boundaryField()[patchI];
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
        // To-do: add a stabilisation traction function to momentumStabilisation
        const scalar scaleFactor =
            readScalar(stabilisation().dict().lookup("scaleFactor"));
        const surfaceTensorField gradDf(fvc::interpolate(gradD()));
        traction += scaleFactor*impKf_*(fvc::snGrad(D()) - (n & gradDf));

        // Calculate the force at the faces
        surfaceVectorField force(magSfCurrent*traction);

        // Enforce traction boundary conditions
        // enforceTractionBoundaries(traction, D, nCurrent);
        enforceTractionBoundaries(force, D(), nCurrent, magSfCurrent);

        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(force)
          + rho()*g()
        );

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
            D().primitiveFieldRef(),
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
        );

        // Map the p field to the SNES solution vector
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            p().primitiveFieldRef(),
            foamPetscSnesHelper::solution(),
            blockSize_ -1 // Location of first component
        );
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the D field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        D().primitiveFieldRef(),
        0, // Location of first component
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    D().correctBoundaryConditions();

    if (solvePressure())
    {
        // Map the PETSc solution to the p field
        // p is located in the last ("blockSize - 1") component
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            foamPetscSnesHelper::solution(),
            p().primitiveFieldRef(),
            blockSize_ - 1 // Location of p component
        );

        p().correctBoundaryConditions();

        // Update dpdt
        dpdtPtr_.ref() = fvc::ddt(p());
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

    return true;
}


void nonLinGeomTotalLagTotalDispSolid::makePDiffusivity() const
{
    if (pDiffusivityPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const scalar pressureSmoothingCoeff
    (
        readScalar(solidModelDict().lookup("pressureSmoothingCoeff"))
    );

    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D())
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(D());
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    pDiffusivityPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "pDiffusivity",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -pressureSmoothingCoeff*impKf_/fvc::interpolate(approxJ.A())
        )
    );
}


const surfaceScalarField& nonLinGeomTotalLagTotalDispSolid::pDiffusivity() const
{
    if (pDiffusivityPtr_.empty())
    {
        makePDiffusivity();
    }

    return pDiffusivityPtr_.ref();
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
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh().nCells(),
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3),
        solvePressure()
      ? (
            solidModel::twoD()
          ? labelListList({ {0, 1}, {2} })    // 0-1: momentum. 2: pressure
          : labelListList({ {0, 1, 2}, {3} }) // 0-2: momentum. 3: pressure
        )
      : labelListList({}), // No need to assign fields
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
    pDiffusivityPtr_(),
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
        mesh().gradScheme("grad(" + D().name() +')')
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
                        D().boundaryFieldRef()[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
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


label nonLinGeomTotalLagTotalDispSolid::initialiseJacobian(Mat jac)
{
    // Initialise based on compact stencil fvMesh
    return Foam::initialiseJacobian(jac, mesh(), blockSize_);
}


label nonLinGeomTotalLagTotalDispSolid::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
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
        blockSize_, // Block size of x
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Enforce the boundary conditions
    D.correctBoundaryConditions();

    // Update gradient of displacement
    mechanical().grad(D, gradD());

    // Enforce the boundary conditions again for any conditions that use gradD
    D.correctBoundaryConditions();

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

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1, blockSize_
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        sigma() = dev(sigma()) - p*I;
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

    //fvc::div(J_*Finv_ & sigma(), "div(sigma)");

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD()));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Calculate the force at the faces
    surfaceVectorField force(magSfCurrent*traction);

    // Enforce traction boundary conditions
    // enforceTractionBoundaries(traction, D, nCurrent);
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
            g() - fvc::d2dt2(D) - dampingCoeff()*fvc::ddt(D)
        )
    );

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
        blockSize_, // Block size of x
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    if (solvePressure())
    {
        volScalarField& p = const_cast<volScalarField&>(this->p());

        // Divided by bulkModulus form
        const volScalarField kappa("kappa", mechanical().bulkModulus());
        const surfaceScalarField kappaf(fvc::interpolate(kappa));
        scalarField pressureResidual
        (
          - p/kappa
          + fvc::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
          - fvc::div
            (
                (pDiffusivity()/kappaf)*mesh.Sf()
              & fvc::interpolate(fvc::grad(p))
            )
          - 0.5*(pow(J_, 2.0) - 1.0)
        //- tr(gradD())
        );

        // Make residual extensive
        pressureResidual *= mesh.V();

        // Copy the pressureResidual into the f field as the final equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1, blockSize_
        );
    }

    return 0;
}


label nonLinGeomTotalLagTotalDispSolid::formJacobian
(
    Mat jac,
    const PetscScalar *x
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
        blockSize_, // Block size of x
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
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
            x, pI, blockSize_ - 1, blockSize_
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();
    }

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D, "laplacian(DD,D)")
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

    if (solvePressure())
    {
        const volScalarField& p = this->p();

        const volScalarField kappa("kappa", mechanical().bulkModulus());
        //const volScalarField rKappa(1.0/mechanical().bulkModulus());
        const volScalarField rKappa(1.0/kappa);
        const surfaceScalarField kappaf(fvc::interpolate(kappa));
        {
            // Calculate pressure equation matrix
            const dimensionedScalar one("one", dimless, 1);
            fvScalarMatrix approxPressureJ
            (
              - fvm::Sp(rKappa, p)
              + fvm::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
            );

            // Insert the pressure equation
            foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            (
                approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            );
        }

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

    return 0;
}


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
