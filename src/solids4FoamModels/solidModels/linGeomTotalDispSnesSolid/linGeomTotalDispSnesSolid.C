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

#include "linGeomTotalDispSnesSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumStabilisation.H"
#include "backwardDdtScheme.H"
#include "sparseMatrixTools.H"
#ifdef USE_PETSC
    #include <petscksp.h>
    #include <petscsnes.h>
#endif

// TESTING
#include "linearElastic.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

#ifdef USE_PETSC

// User data "context" for PETSc functions
typedef struct appCtx
{
    // Reference to the solid model object
    Foam::solidModels::linGeomTotalDispSnesSolid& solMod_;

    // Constructor
    appCtx
    (
        Foam::solidModels::linGeomTotalDispSnesSolid& solMod
    )
    :
        solMod_(solMod)
    {}
} appCtx;


PetscErrorCode formResidualLinGeomTotalDispSnesSolid
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
    // VecRestoreArrayRead() and VecRestoreArray() must be called when access is
    // no longer needed
    CHKERRQ(VecGetArrayRead(x,&xx));
    CHKERRQ(VecGetArray(f,&ff));

    // Map the solution to an OpenFOAM field
    const bool twoD = user->solMod_.twoD();
    // Foam::volVectorField D("petsc_D", user->solMod_.D());
    // Foam::volVectorField D("D", user->solMod_.D());
    Foam::volVectorField& D = user->solMod_.D();
    Foam::vectorField& DI = D;
    // const Foam::boolList& ownedByThisProc =
    //     user->solMod_.gPointIndices().ownedByThisProc();
    if (Foam::Pstream::parRun())
    {
        // globalPoint should do what I need
        Foam::FatalError
            << "Fix parallel" << Foam::abort(Foam::FatalError);
    }
    {
        int index = 0;
        forAll(DI, i)
        {
            //if (ownedByThisProc[i])
            {
                DI[i].x() = xx[index++];
                DI[i].y() = xx[index++];

                if (!twoD)
                {
                    DI[i].z() = xx[index++];
                }
            }
        }
    }


    // Compute the residual
    const Foam::vectorField res(user->solMod_.residualMomentum(D));

    // Map the data to f
    // const Foam::labelList& localToGlobalPointMap =
    //     user->solMod_.gPointIndices().localToGlobalPointMap();
    const int blockSize = twoD ? 2 : 3;
    forAll(res, localBlockRowI)
    {
        const Foam::vector& resI = res[localBlockRowI];
        //const int blockRowI = localToGlobalPointMap[localBlockRowI];
        const int blockRowI = localBlockRowI; // serial only

        ff[blockRowI*blockSize] = resI.x();
        ff[blockRowI*blockSize + 1] = resI.y();

        if (!twoD)
        {
            ff[blockRowI*blockSize + 2] = resI.z();
        }
    }

    // Restore vectors
    CHKERRQ(VecRestoreArrayRead(x,&xx));
    CHKERRQ(VecRestoreArray(f,&ff));

    return 0;
}


PetscErrorCode formJacobianLinGeomTotalDispSnesSolid
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Mat jac,      // Jacobian
    Mat B,        // Jaconian precondioner (can be jac)
    void *ctx     // user context
)
{
    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Map the solution to an OpenFOAM field
    appCtx *user = (appCtx *)ctx;
    const bool twoD = user->solMod_.twoD();
    // Foam::volVectorField D("petsc_D", user->solMod_.D());
    // Foam::volVectorField D("D", user->solMod_.D());
    Foam::volVectorField& D = user->solMod_.D();
    Foam::vectorField& DI = D;
    // const Foam::boolList& ownedByThisProc =
    //     user->solMod_.gPointIndices().ownedByThisProc();
    if (Foam::Pstream::parRun())
    {
        // globalPoint should do what I need
        Foam::FatalError
            << "Fix parallel" << Foam::abort(Foam::FatalError);
    }
    {
        int index = 0;
        forAll(DI, i)
        {
            // if (ownedByThisProc[i])
            {
                DI[i].x() = xx[index++];
                DI[i].y() = xx[index++];

                if (!twoD)
                {
                    DI[i].z() = xx[index++];
                }
            }
        }
    }

    // Correct boundaries
    // Do I need a traction loop here?
    D.correctBoundaryConditions();

    // Restore solution vector
    CHKERRQ(VecRestoreArrayRead(x, &xx));


    // Compute Jacobian in OpenFOAM format
    Foam::sparseMatrix matrix;
    matrix += user->solMod_.JacobianMomentum(D)();


    // Set matrix coefficients, if any, to zero
    // TODO: only set the matrix once since it does not change!!!
    MatInfo info;
    MatGetInfo(B, MAT_LOCAL, &info);
    if (info.nz_used)
    {
        // Foam::Info<< "Zeroing the matrix" << Foam::endl;
        // If we don't change the matrix then we should return it
        CHKERRQ(MatZeroEntries(B));
    }
    else
    {
        const int blockSize = twoD ? 2 : 3;
        Foam::Info<< "Initialising the matrix" << Foam::endl;
        Foam::labelList nonZerosPerBlockRow(user->solMod_.mesh().nCells(), 1);
        forAll(user->solMod_.mesh().owner(), faceI)
        {
            nonZerosPerBlockRow[user->solMod_.mesh().owner()[faceI]]++;
            nonZerosPerBlockRow[user->solMod_.mesh().neighbour()[faceI]]++;
        }
        const int d_nz = blockSize*Foam::max(nonZerosPerBlockRow);
        Foam::Info<< "    d_nz = " << d_nz << Foam::endl;

        MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        // Allocate parallel matrix
        //CHKERRQ(MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz));
        // Allocate parallel matrix with the same conservative stencil per node
        CHKERRQ(MatMPIAIJSetPreallocation(B, d_nz, NULL, 0, NULL));
    }

    // Todo: do I need to set the number of non-zeros here?
    // This will majorly impact performance

    // Insert OpenFOAM matrix into PETSc matrix
    // Note: we use global indices when inserting coefficients
    const Foam::sparseMatrixData& data = matrix.data();
    // const Foam::labelList& localToGlobalPointMap =
    //     user->solMod_.gPointIndices().localToGlobalPointMap();
    PetscScalar values2d[4];
    PetscScalar values3d[9];
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        const Foam::tensor& coeff = iter();
        // const int blockRowI = localToGlobalPointMap[iter.key()[0]];
        // const int blockColI = localToGlobalPointMap[iter.key()[1]];
        const int blockRowI = iter.key()[0];
        const int blockColI = iter.key()[1];

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
                    B, 1, &blockRowI, 1, &blockColI, values2d, ADD_VALUES
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
                    B, 1, &blockRowI, 1, &blockColI, values3d, ADD_VALUES
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

defineTypeNameAndDebug(linGeomTotalDispSnesSolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomTotalDispSnesSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void linGeomTotalDispSnesSolid::predict()
{
    Info<< "Applying linear predictor to D" << endl;

    // Predict D using previous time steps
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomTotalDispSnesSolid::linGeomTotalDispSnesSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false))
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> linGeomTotalDispSnesSolid::residualMomentum
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
    }

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


tmp<sparseMatrix> linGeomTotalDispSnesSolid::JacobianMomentum
(
    const volVectorField& D
)
{
    // Count the number of non-zeros for a Laplacian discretisation
    // This equals the sum of one plus the number of internal faces for each,
    // which can be calculated as nCells + 2*nInternalFaces
    // Multiply by 3 since we will form the block matrix
    // Todo: look 2D vs 3D
    const label numNonZeros =
        3*returnReduce
        (
            mesh().nCells() + 2.0*mesh().nInternalFaces(), sumOp<label>()
        );
    // globalIndex globalCells(mesh.nCells());
    // label globalCellI = globalCells.toGlobal(cellI);

    // const surfaceVectorField Df(fvc::interpolate(D));
    // const dimensionedScalar dimLL("1", dimless/pow(dimLength, 2), 1);
    // surfaceScalarField scaleFac(1.0 + (Df & Df)*dimLL);

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

    // Optional under-relaxation of the linear system
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

            matrix(blockRowI, blockRowI) = coeff;
        }
    }

    // Insert the off-diagonal
    {
        const labelUList& own = mesh().owner();
        const labelUList& nei = mesh().neighbour();
        const scalarField& upper = approxJ.upper();
        forAll(own, faceI)
        {
            const label blockRowI = own[faceI];
            const label blockColI = nei[faceI];
            const tensor coeff(upper[faceI]*I);

            matrix(blockRowI, blockColI) = coeff;
            matrix(blockColI, blockRowI) = coeff;
        }
    }

    return tmatrix;
}


bool linGeomTotalDispSnesSolid::evolve()
{
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    if (predictor_)
    {
        predict();
    }

    // Store fields for under-relaxation and residual calculation
    D().storePrevIter();

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Create user data context
    appCtx user(*this);

    // Lookup the PETSc options file
    fileName optionsFile(solidModelDict().lookup("optionsFile"));
    optionsFile.expand();

    // Initialise PETSc with an options file
    PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

    // Create a SNES solver context
    SNES snes;
    SNESCreate(PETSC_COMM_WORLD, &snes);

    // Set the user context
    SNESSetApplicationContext(snes, &user);

    // Set the residual function
    SNESSetFunction(snes, NULL, formResidualLinGeomTotalDispSnesSolid, &user);

    // Initialise the Jacobian matrix

    const int blockSize = twoD_ ? 2 : 3;
    if (Pstream::parRun())
    {
        FatalError
            << "evolve to be fixed for parallel" << abort(FatalError);
    }
    // const Foam::labelList& localToGlobalPointMap =
    //     globalPointIndices_.localToGlobalPointMap();

    // Find size of global system, i.e. the highest global cell index + 1
    //const int blockN = Foam::gMax(localToGlobalPointMap) + 1;
    const int blockN = mesh().nCells(); // fix for parallel
    const int N = blockSize*blockN;

    // Find the start and end global point indices for this proc
    const int blockStartID = 0; // fix for parallel
    const int blockEndID = mesh().nCells() - 1; // fix for parallel
    // int blockStartID = N;
    // int blockEndID = -1;
    // const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    // forAll(ownedByThisProc, pI)
    // {
    //     if (ownedByThisProc[pI])
    //     {
    //         blockStartID = Foam::min(blockStartID, localToGlobalPointMap[pI]);
    //         blockEndID = Foam::max(blockEndID, localToGlobalPointMap[pI]);
    //     }
    // }
    //const int startID = blockSize*blockStartID;
    //const int endID = blockSize*(blockEndID + 1) - 1;

    // Find size of local system, i.e. the range of points owned by this proc
    const int blockn = blockEndID - blockStartID + 1;
    const int n = blockSize*blockn;

    // Create the matrix
    Mat A;
    CHKERRQ(MatCreate(PETSC_COMM_WORLD, &A));

    // Set matrix characteristics
    CHKERRQ(MatSetSizes(A, n, n, N, N));
    //CHKERRQ(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
    CHKERRQ(MatSetFromOptions(A));
    CHKERRQ(MatSetType(A, MATMPIAIJ));
    //CHKERRQ(MatSetType(A, MATMPIBAIJ));
    CHKERRQ(MatSetBlockSize(A, blockSize));

    // Count the number of non-zeros
    // TODO: fix for parallel
    // int* d_nnz = (int*)malloc(n*sizeof(int));
    // int* o_nnz = (int*)malloc(n*sizeof(int));
    // label d_nnz[n];
    // label o_nnz[n];
    // int d_nz = 0;
    // Foam::sparseMatrixTools::setNonZerosPerRow
    // (
    //     d_nnz,
    //     o_nnz,
    //     d_nz,
    //     n,
    //     blockSize,
    //     ownedByThisProc,
    //     globalPointIndices_.stencilSizeOwned(),
    //     globalPointIndices_.stencilSizeNotOwned()
    // );
    labelList nonZerosPerBlockRow(mesh().nCells(), 1);
    forAll(mesh().owner(), faceI)
    {
        nonZerosPerBlockRow[mesh().owner()[faceI]]++;
        nonZerosPerBlockRow[mesh().neighbour()[faceI]]++;
    }
    const int d_nz = blockSize*max(nonZerosPerBlockRow);

    // Allocate parallel matrix
    //CHKERRQ(MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz));
    // Allocate parallel matrix with the same conservative stencil per node
    CHKERRQ(MatMPIAIJSetPreallocation(A, d_nz, NULL, 0, NULL));
    // TODO: change d_nnz/o_nnz to block sizes!
    //CHKERRQ(MatMPIBAIJSetPreallocation(A, blockSize, 0, d_nnz, 0, o_nnz));

    // TO BE FIXED: some mallocs are still needed in parallel!
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    CHKERRQ(MatSetUp(A));

    // Do not call the matrix assembly as we have not inserted any values
    //CHKERRQ(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    //CHKERRQ(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    // Set the Jacobian function
    SNESSetJacobian(snes, A, A, formJacobianLinGeomTotalDispSnesSolid, &user);

    // Set solver options
    // Uses default options, can be overridden by command line options
    SNESSetFromOptions(snes);

    // // Find the start and end global point indices for this proc
    // forAll(ownedByThisProc, pI)
    // {
    //     if (ownedByThisProc[pI])
    //     {
    //         blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
    //         blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
    //     }
    // }
    //const label startID = blockSize*blockStartID;
    //const label endID = blockSize*(blockEndID + 1) - 1;

    // Create the solution vector
    Vec x;
    CHKERRQ(VecCreate(PETSC_COMM_WORLD, &x));
    CHKERRQ(VecSetSizes(x, n, N));
    CHKERRQ(VecSetBlockSize(x, blockSize));
    CHKERRQ(VecSetType(x, VECMPI));
    CHKERRQ(PetscObjectSetName((PetscObject) x, "Solution"));
    CHKERRQ(VecSetFromOptions(x));

    // Set initial guess
    CHKERRQ(VecSet(x, 0.0));

    // Solve the nonlinear system
    CHKERRQ(SNESSolve(snes, NULL, x));

    // Check if SNES converged
    SNESConvergedReason reason;
    CHKERRQ(SNESGetConvergedReason(snes, &reason));
    if (reason < 0)
    {
        Warning
            << "PETSc SNES solver return error check disabled" << endl
            << "The SNES nonlinear solver did not converge." << nl
            << " PETSc SNES convergence error code: " << reason << nl
            << " PETSc SNES convergence reason: "
            << SNESConvergedReasons[reason] << endl;
        // FatalErrorIn("vertexCentredLinGeomSolid::evolveSnes()")
        //     << "The SNES nonlinear solver did not converge."
        //     << " PETSc SNES error code = " << reason << abort(FatalError);
    }

    // Retrieve the solution
    const PetscScalar *xx;
    VecGetArrayRead(x,&xx);
    Foam::vectorField& DI = D();
    {
        int index = 0;
        forAll(DI, i)
        {
            // if (ownedByThisProc[i])
            {
                DI[i].x() = xx[index++];
                DI[i].y() = xx[index++];

                if (!twoD_)
                {
                    DI[i].z() = xx[index++];
                }
            }
        }
    }
    D().correctBoundaryConditions();

    // Destroy PETSc objects
    VecDestroy(&x);
    SNESDestroy(&snes);

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


tmp<vectorField> linGeomTotalDispSnesSolid::tractionBoundarySnGrad
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
