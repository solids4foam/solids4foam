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

#include "linGeomTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumStabilisation.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"


// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

#ifdef USE_PETSC

// // User data "context" for PETSc functions
// typedef struct appCtx
// {
//     // Reference to the solid model object
//     Foam::solidModels::linGeomTotalDispSolid& solMod_;

//     // Constructor
//     appCtx
//     (
//         Foam::solidModels::linGeomTotalDispSolid& solMod
//     )
//     :
//         solMod_(solMod)
//     {}
// } appCtx;


PetscErrorCode formResidualLinGeomTotalDispSolid
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Vec f,        // residual
    void *ctx     // user context
)
{
    const PetscScalar *xx;
    PetscScalar       *ff;
    appCtx *user = (appCtx *)ctx;

    // Access x and f data
    CHKERRQ(VecGetArrayRead(x, &xx));
    CHKERRQ(VecGetArray(f, &ff));

    // Map the solution to an OpenFOAM field
    const bool twoD = user->solMod_.twoD();
    const int blockSize = twoD ? 2 : 3;
    Foam::volVectorField& D = user->solMod_.D();
    Foam::vectorField& DI = D;

    {
        int index = 0;
        forAll(DI, localCellI)
        {
            DI[localCellI][Foam::vector::X] = xx[index++];
            DI[localCellI][Foam::vector::Y] = xx[index++];

            if (!twoD)
            {
                DI[localCellI][Foam::vector::Z] = xx[index++];
            }
        }
    }

    // Restore the solution vector
    CHKERRQ(VecRestoreArrayRead(x, &xx));

    // Enforce the D boundary conditions and sync processor boundaries
    D.correctBoundaryConditions();

    // Compute the residual
    const Foam::vectorField res(user->solMod_.residualMomentum(D));

    // Map the data to f
    forAll(res, localBlockRowI)
    {
        const Foam::vector& resI = res[localBlockRowI];
        const int blockRowI = localBlockRowI;

        ff[blockRowI*blockSize] = resI.x();
        ff[blockRowI*blockSize + 1] = resI.y();

        if (!twoD)
        {
            ff[blockRowI*blockSize + 2] = resI.z();
        }
    }

    // Restore the source vector
    CHKERRQ(VecRestoreArray(f, &ff));

    return 0;
}


PetscErrorCode formJacobianLinGeomTotalDispSolid
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Mat jac,      // Jacobian
    Mat B,        // Preconditioner matrix (can be jac)
    void *ctx     // user context
)
{
    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Map the solution to an OpenFOAM field
    appCtx *user = (appCtx *)ctx;
    const bool twoD = user->solMod_.twoD();
    const int blockSize = twoD ? 2 : 3;
    Foam::volVectorField& D = user->solMod_.D();
    Foam::vectorField& DI = D;

    {
        int index = 0;
        forAll(DI, localCellI)
        {
            DI[localCellI][Foam::vector::X] = xx[index++];
            DI[localCellI][Foam::vector::Y] = xx[index++];

            if (!twoD)
            {
                DI[localCellI][Foam::vector::Z] = xx[index++];
            }
        }
    }

    // Enforce boundary conditions
    D.correctBoundaryConditions();

    // Restore solution vector
    CHKERRQ(VecRestoreArrayRead(x, &xx));


    // Compute Jacobian in OpenFOAM format
    Foam::sparseMatrix matrix;
    matrix += user->solMod_.JacobianMomentum(D)();

    // Initialise the matrix if it has yet to be allocated; otherwise zero all
    // entries
    MatInfo info;
    MatGetInfo(B, MAT_LOCAL, &info);
    if (info.nz_used)
    {
        // Foam::Info<< "Zeroing the Jacobian" << Foam::endl;
        // Zero the matrix but do not reallocate the space
        // The "-snes_lag_jacobian -2" PETSc option can be used to avoid
        // re-building the matrix
        CHKERRQ(MatZeroEntries(B));
    }
    else
    {
        Foam::Info<< "    Initialising the matrix" << Foam::endl;

        // Set the block size
        CHKERRQ(MatSetBlockSize(B, blockSize));

        // Number of vector unknowns on this processor
        const int blockn = user->solMod_.globalCells().localSize();

        // Count the number of non-zeros in the matrix

        // Number of on-processor non-zeros per row
        int* d_nnz = (int*)malloc(blockn*sizeof(int));

        // Number of off-processor non-zeros per row
        int* o_nnz = (int*)malloc(blockn*sizeof(int));

        // Initialise d_nnz and o_nnz to zero
        for (int i = 0; i < blockn; ++i)
        {
            d_nnz[i] = 1; // count diagonal cell
            o_nnz[i] = 0;
        }

        // Count neighbours sharing an internal face
        const Foam::labelUList& own = user->solMod_.mesh().owner();
        const Foam::labelUList& nei = user->solMod_.mesh().neighbour();
        forAll(own, faceI)
        {
            const Foam::label ownCellID = own[faceI];
            const Foam::label neiCellID = nei[faceI];
            d_nnz[ownCellID]++;
            d_nnz[neiCellID]++;
        }

        // Count off-processor neighbour cells
        forAll(user->solMod_.mesh().boundary(), patchI)
        {
            if (user->solMod_.mesh().boundary()[patchI].type() == "processor")
            {
                const Foam::unallocLabelList& faceCells =
                    user->solMod_.mesh().boundary()[patchI].faceCells();

                forAll(faceCells, fcI)
                {
                    const Foam::label cellID = faceCells[fcI];
                    o_nnz[cellID]++;
                }
            }
            else if (user->solMod_.mesh().boundary()[patchI].coupled())
            {
                // Other coupled boundaries are not implemented
                Foam::FatalError
                    << "Coupled boundary are not implemented, except for"
                    << " processor boundaries" << Foam::abort(Foam::FatalError);
            }
        }

        // Allocate parallel matrix
        //CHKERRQ(MatMPIAIJSetPreallocation(B, 0, d_nnz, 0, o_nnz));
        // Allocate parallel matrix with the same conservative stencil per node
        //CHKERRQ(MatMPIAIJSetPreallocation(B, d_nz, NULL, 0, NULL));
        CHKERRQ(MatMPIBAIJSetPreallocation(B, blockSize, 0, d_nnz, 0, o_nnz));

        // Raise an error if mallocs are required during matrix assembly
        MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }

    // Insert OpenFOAM matrix into PETSc matrix
    // Note: the sparseMatrix object contains global indices, which we use when
    // inserting coefficients
    const Foam::sparseMatrixData& data = matrix.data();
    PetscScalar values2d[4];
    PetscScalar values3d[9];
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        const Foam::tensor& coeff = iter();
        const Foam::label globalBlockRowI = iter.key()[0];
        const Foam::label globalBlockColI = iter.key()[1];

        if (twoD)
        {
            // Prepare values
            values2d[0] = coeff.xx();
            values2d[1] = coeff.xy();
            values2d[2] = coeff.yx();
            values2d[3] = coeff.yy();

            // Insert tensor coefficient
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    B, 1, &globalBlockRowI, 1, &globalBlockColI, values2d,
                    ADD_VALUES
                )
            );
        }
        else // 3-D
        {
            // Prepare values
            values3d[0] = coeff.xx();
            values3d[1] = coeff.xy();
            values3d[2] = coeff.xz();
            values3d[3] = coeff.yx();
            values3d[4] = coeff.yy();
            values3d[5] = coeff.yz();
            values3d[6] = coeff.zx();
            values3d[7] = coeff.zy();
            values3d[8] = coeff.zz();

            // Insert tensor coefficient
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    B, 1, &globalBlockRowI, 1, &globalBlockColI, values3d,
                    ADD_VALUES
                )
            );
        }
    }

    // Complete matrix assembly
    CHKERRQ(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
    if (jac != B)
    {
        CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
    }

    return 0;
}
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomTotalDispSolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomTotalDispSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void linGeomTotalDispSolid::predict()
{
    Info<< "Applying linear predictor to D" << endl;

    // Predict D using previous time steps
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


void linGeomTotalDispSolid::enforceTractionBoundaries
(
    surfaceVectorField& traction,
    const volVectorField& D,
    const surfaceVectorField& n
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

            const vectorField& nPatch = n.boundaryField()[patchI];

            traction.boundaryFieldRef()[patchI] =
                tracPatch.traction() - nPatch*tracPatch.pressure();
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
            // Unit normals
            const vectorField& nPatch = n.boundaryField()[patchI];

            // Set shear traction to zero
            traction.boundaryFieldRef()[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
        }
    }
}


bool linGeomTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver using an implicit segregated approach"
        << endl;

    if (predictor_ && newTimeStep())
    {
        predict();
    }

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        SolverPerformance<vector> solverPerfD;
        SolverPerformance<vector>::debug = 0;
#else
        lduSolverPerformance solverPerfD;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the momentum equation for D" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DEqn
            (
                rho()*fvm::d2dt2(D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(D(), gradD(), impK_)
            );

            // Add damping
            if (dampingCoeff().value() > SMALL)
            {
                DEqn += dampingCoeff()*rho()*fvm::ddt(D());
            }

            // Under-relaxation the linear system
            DEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DEqn);

            // Solve the linear system
            solverPerfD = DEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Update the momentum equation inverse diagonal field
            // This may be used by the mechanical law when calculating the
            // hydrostatic pressure
            const volScalarField DEqnA("DEqnA", DEqn.A());

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());
        }
        while
        (
            !converged
            (
                iCorr,
#ifdef OPENFOAM_NOT_EXTEND
                mag(solverPerfD.initialResidual()),
                cmptMax(solverPerfD.nIterations()),
#else
                solverPerfD.initialResidual(),
                solverPerfD.nIterations(),
#endif
                D()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), gradD(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


bool linGeomTotalDispSolid::evolveSnes()
{
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Map the D field to the SNES solution vector
        Vec& x = xPtr_();
        PetscScalar *xx;
        VecGetArray(x, &xx);
        const Foam::vectorField& DI = D();
        {
            int index = 0;
            forAll(DI, i)
            {
                xx[index++] = DI[i][vector::X];
                xx[index++] = DI[i][vector::Y];

                if (!solidModel::twoD())
                {
                    xx[index++] = DI[i][vector::Z];
                }
            }
        }

        // Restore solution vector
        CHKERRQ(VecRestoreArray(x, &xx));
    }

    // Take a reference to the SNES solver object
    SNES& snes = snesPtr_();

    // Take a reference to the solution vector
    Vec& x = xPtr_();

    // Solve the nonlinear system
    CHKERRQ(SNESSolve(snes, NULL, x));

    // Check if SNES converged
    SNESConvergedReason reason;
    CHKERRQ(SNESGetConvergedReason(snes, &reason));
    if (reason < 0)
    {
        WarningIn("bool linGeomTotalDispSolid::evolveSnes()")
            << "PETSc SNES solver return error check disabled" << endl
            << "The SNES nonlinear solver did not converge." << nl
            << " PETSc SNES convergence error code: " << reason << nl
            << " PETSc SNES convergence reason: "
            << SNESConvergedReasons[reason] << endl;

        if (stopOnPetscError_)
        {
            FatalErrorIn("bool linGeomTotalDispSolid::evolveSnes()")
                << "Stopping because of the PETSc SNES error"
                << "Set `stopOnPetscError` to `false` to continue on PETSc "
                << "SNES errors"
                << abort(FatalError);
        }
    }

    // Retrieve the solution
    {
        const PetscScalar *xx;
        VecGetArrayRead(x, &xx);
        Foam::vectorField& DI = D();
        {
            int index = 0;
            forAll(DI, i)
            {
                DI[i].x() = xx[index++];
                DI[i].y() = xx[index++];

                if (!solidModel::twoD())
                {
                    DI[i].z() = xx[index++];
                }
            }
        }
        D().correctBoundaryConditions();

        // Restore solution vector
        CHKERRQ(VecRestoreArrayRead(x, &xx));
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    return true;
}


bool linGeomTotalDispSolid::evolveExplicit()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the solid momentum equation for D using an explicit "
            << "approach" << nl
            << "Simulation Time, Clock Time, Max Stress" << endl;
    }

    physicsModel::printInfo() = bool
    (
        time().timeIndex() % infoFrequency() == 0
     || mag(time().value() - time().endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< time().value() << " " << time().elapsedClockTime()
            << " " << max(mag(sigma())).value() << endl;

        physicsModel::printInfo() = false;
    }

    // Take references for brevity and efficiency
    const fvMesh& mesh = solidModel::mesh();
    volVectorField& D = solidModel::D();
    volTensorField& gradD = solidModel::gradD();
    volVectorField& U = solidModel::U();
    volSymmTensorField& sigma = solidModel::sigma();
    const volScalarField& rho = solidModel::rho();

    // Central difference scheme

    // Take a reference to the current and previous time-step
    const dimensionedScalar& deltaT = time().deltaT();
    //const dimensionedScalar& deltaT0 = time().deltaT0();

    // Compute the velocity
    // Note: this is the velocity at the middle of the time-step
    //pointU_ = pointU_.oldTime() + 0.5*(deltaT + deltaT0)*pointA_.oldTime();
    U = U.oldTime() + deltaT*A_.oldTime();

    // Compute displacement
    D = D.oldTime() + deltaT*U;

    // Enforce boundary conditions on the displacement field
    D.correctBoundaryConditions();

    if (solidModel::twoD())
    {
        // Remove displacement in the empty directions
        forAll(mesh.geometricD(), dirI)
        {
            if (mesh.geometricD()[dirI] < 0)
            {
                D.primitiveFieldRef().replace(dirI, 0.0);
            }
        }
    }

    // Update gradient of displacement
    mechanical().grad(D, gradD);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma);

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Calculate the traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // Solve the momentum equation for acceleration
    A_ = fvc::div(mesh.magSf()*traction)/rho
       + g()
       - dampingCoeff()*fvc::ddt(D);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomTotalDispSolid::linGeomTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    A_
    (
        IOobject
        (
            "A",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength/pow(dimTime, 2), vector::zero)
    ),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    globalCells_(mesh().nCells()),
    stopOnPetscError_
    (
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true)
    ),
    snesPtr_(),
    xPtr_(),
    APtr_(),
    snesUserPtr_()
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // For consistent restarts, we will calculate the gradient field
    D().correctBoundaryConditions();
    D().storePrevIter();
    mechanical().grad(D(), gradD());

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

        // Lookup the PETSc options file
        fileName optionsFile(solidModelDict().lookup("optionsFile"));
        optionsFile.expand();

        // Initialise PETSc with an options file
        PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

        // Create the PETSc SNES object
        snesPtr_.set(new SNES());
        SNES& snes = snesPtr_();
        SNESCreate(PETSC_COMM_WORLD, &snes);

        // Create user data context
        snesUserPtr_.set(new appCtx(*this));
        appCtx& user = snesUserPtr_();

        // Set the user context
        SNESSetApplicationContext(snes, &user);

        // Set the residual function
        SNESSetFunction(snes, NULL, formResidualLinGeomTotalDispSolid, &user);

        // Coefficient block size, e.g. 2x2 in 2-D, 3x3 in 3-D
        const int blockSize = solidModel::twoD() ? 2 : 3;

        // Global system size
        const int blockN = globalCells_.size();
        const int N = blockSize*blockN;

        // Local (this processor) system size
        const int blockn = globalCells_.localSize();
        const int n = blockSize*blockn;

        // Create the Jacobian matrix
        APtr_.set(new Mat());
        Mat& A = APtr_();
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, n, n, N, N);
        MatSetFromOptions(A);
        MatSetType(A, MATMPIAIJ);
        //MatSetType(A, MATMPIBAIJ);

        // Set the Jacobian function
        SNESSetJacobian(snes, A, A, formJacobianLinGeomTotalDispSolid, &user);

        // Set solver options
        // Uses default options, can be overridden by command line options
        SNESSetFromOptions(snes);

        // Create the solution vector
        xPtr_.set(new Vec());
        Vec& x = xPtr_();
        VecCreate(PETSC_COMM_WORLD, &x);
        VecSetSizes(x, n, N);
        VecSetBlockSize(x, blockSize);
        VecSetType(x, VECMPI);
        PetscObjectSetName((PetscObject) x, "Solution");
        VecSetFromOptions(x);
        VecZeroEntries(x);
    }
    else if (solutionAlg() != solutionAlgorithm::EXPLICIT)
    {
        if (gradDScheme == "leastSquaresS4f")
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should only be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " and "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithms" << abort(FatalError);
        }

    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void linGeomTotalDispSolid::setDeltaT(Time& runTime)
{
    if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        // Max wave speed in the domain
        const scalar waveSpeed = max
        (
            Foam::sqrt(mechanical().impK()/mechanical().rho())
        ).value();

        // deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
        // In the current discretisation, information can move two cells per
        // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
        // calculating the required stable time-step
        // i.e. deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
        // For safety, we should use a time-step smaller than this e.g. Abaqus uses
        // stableTimeStep/sqrt(2): we will default to this value
        const scalar requiredDeltaT =
            1.0/
            gMax
            (
                DimensionedField<scalar, Foam::surfaceMesh>
                (
                    mesh().surfaceInterpolation::
                    deltaCoeffs().internalField()
                   *waveSpeed
                )
            );

        // Lookup the desired Courant number
        const scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

        const scalar newDeltaT = maxCo*requiredDeltaT;

        // Update print info
        physicsModel::printInfo() = bool
        (
            runTime.timeIndex() % infoFrequency() == 0
         || mag(runTime.value() - runTime.endTime().value()) < SMALL
        );

        physicsModel::printInfo() = false;

        if (time().timeIndex() == 1)
        {
            Info<< nl << "Setting deltaT = " << newDeltaT
                << ", maxCo = " << maxCo << endl;
        }

        runTime.setDeltaT(newDeltaT);
    }
}


bool linGeomTotalDispSolid::evolve()
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
    else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
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
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::EXPLICIT
               ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


tmp<vectorField> linGeomTotalDispSolid::residualMomentum
(
    const volVectorField& D
)
{
    // Prepare result
    tmp<vectorField> tresidual(new vectorField(D.size(), vector::zero));
    vectorField& residual = tresidual.ref();

    // Enforce the boundary conditions
    const_cast<volVectorField&>(D).correctBoundaryConditions();

    // Update gradient of displacement
    mechanical().grad(D, gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma()));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD()));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    residual =
        fvc::div(mesh().magSf()*traction)
      + rho()
       *(
            g() - fvc::d2dt2(D) - dampingCoeff()*fvc::ddt(D)
        );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh().V();

    return tresidual;
}


tmp<sparseMatrix> linGeomTotalDispSolid::JacobianMomentum
(
    const volVectorField& D
)
{
    // Count the number of non-zeros for a Laplacian discretisation
    // This equals the sum of one plus the number of internal faces for each,
    // which can be calculated as nCells + 2*nInternalFaces
    // Multiply by the blockSize since we will form the block matrix
    const int blockSize = solidModel::twoD() ? 2 : 3;
    const label numNonZeros =
        blockSize*returnReduce
        (
            mesh().nCells() + 2.0*mesh().nInternalFaces(), sumOp<label>()
        );

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

    // Convert fvMatrix matrix to sparseMatrix

    // Initialise matrix
    tmp<sparseMatrix> tmatrix(new sparseMatrix(numNonZeros));
    sparseMatrix& matrix = tmatrix.ref();

    // Insert the diagonal
    {
        const vectorField diag(approxJ.DD());
        forAll(diag, blockRowI)
        {
            const tensor coeff
            (
                diag[blockRowI][vector::X], 0, 0,
                0, diag[blockRowI][vector::Y], 0,
                0,  0, diag[blockRowI][vector::Z]
            );

            const label globalBlockRowI = globalCells_.toGlobal(blockRowI);

            matrix(globalBlockRowI, globalBlockRowI) = coeff;
        }
    }

    // Insert the off-diagonal
    {
        const labelUList& own = mesh().owner();
        const labelUList& nei = mesh().neighbour();
        const scalarField& upper = approxJ.upper();
        forAll(own, faceI)
        {
            const tensor coeff(upper[faceI]*I);

            const label blockRowI = own[faceI];
            const label blockColI = nei[faceI];

            const label globalBlockRowI = globalCells_.toGlobal(blockRowI);
            const label globalBlockColI = globalCells_.toGlobal(blockColI);

            matrix(globalBlockRowI, globalBlockColI) = coeff;
            matrix(globalBlockColI, globalBlockRowI) = coeff;
        }
    }

    // Collect the global cell indices from neighbours at processor boundaries
    // These are used to insert the off-processor coefficients
    // First, send the data
    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pD = D.boundaryField()[patchI];
        if (pD.type() == "processor")
        {
            // Take a copy of the faceCells (local IDs) and convert them to
            // global IDs
            labelList globalFaceCells(mesh().boundary()[patchI].faceCells());
            globalCells_.inplaceToGlobal(globalFaceCells);

            // Send global IDs to the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh().boundary()[patchI]);
            // procPatch.compressedSend
            procPatch.send
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );
        }
    }
    // Next, receive the data
    PtrList<labelList> neiProcGlobalIDs(D.boundaryField().size());
    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pD = D.boundaryField()[patchI];
        if (pD.type() == "processor")
        {
            neiProcGlobalIDs.set(patchI, new labelList(pD.size()));
            labelList& globalFaceCells = neiProcGlobalIDs[patchI];

            // Receive global IDs from the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh().boundary()[patchI]);
            // procPatch.compressedReceive
            procPatch.receive
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );
        }
    }

    // Insert the off-processor coefficients
    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pD = D.boundaryField()[patchI];

        if (pD.type() == "processor")
        {
            const vectorField& intCoeffs = approxJ.internalCoeffs()[patchI];
            const vectorField& neiCoeffs = approxJ.boundaryCoeffs()[patchI];
            const unallocLabelList& faceCells =
                mesh().boundary()[patchI].faceCells();
            const labelList& neiGlobalFaceCells = neiProcGlobalIDs[patchI];

            forAll(pD, faceI)
            {
                const label globalBlockRowI =
                    globalCells_.toGlobal(faceCells[faceI]);

                // On-proc diagonal coefficient
                {
                    const tensor coeff
                    (
                        intCoeffs[faceI][vector::X], 0, 0,
                        0, intCoeffs[faceI][vector::Y], 0,
                        0, 0, intCoeffs[faceI][vector::Z]
                    );

                    matrix(globalBlockRowI, globalBlockRowI) += coeff;
                }

                // Off-proc off-diagonal coefficient
                {
                    const tensor coeff
                    (
                        neiCoeffs[faceI][vector::X], 0, 0,
                        0, neiCoeffs[faceI][vector::Y], 0,
                        0, 0, neiCoeffs[faceI][vector::Z]
                    );

                    const label globalBlockColI = neiGlobalFaceCells[faceI];

                    matrix(globalBlockRowI, globalBlockColI) += coeff;
                }
            }
        }
        else if (pD.coupled()) // coupled but not a processor boundary
        {
            FatalErrorIn
            (
                "tmp<sparseMatrix> linGeomTotalDispSolid::JacobianMomentum"
            )   << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        // else non-coupled boundary contributions have already been added to
        // the diagonal
    }

    return tmatrix;
}


tmp<vectorField> linGeomTotalDispSolid::tractionBoundarySnGrad
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
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
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


void linGeomTotalDispSolid::end()
{
    // Destroy PETSc objects
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        SNESDestroy(&snesPtr_());
        snesPtr_.clear();

        VecDestroy(&xPtr_());
        xPtr_.clear();

        MatDestroy(&APtr_());
        APtr_.clear();

        snesUserPtr_.clear();

        PetscFinalize();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
