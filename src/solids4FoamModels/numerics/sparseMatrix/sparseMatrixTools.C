/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sparseMatrixTools.H"
#include "OFstream.H"
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <petscksp.h>

// * * * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * * * //

bool Foam::sparseMatrixTools::checkTwoD(const polyMesh& mesh)
{
    const Vector<label>& geomD = mesh.geometricD();

    label nDir = 0;
    forAll(geomD, dirI)
    {
        if (geomD[dirI] > SMALL)
        {
            nDir++;
        }
    }

    bool twoD = false;
    if (nDir == 2)
    {
        if (mesh.solutionD()[vector::Z] > -1)
        {
            FatalErrorIn("checkTwoD()")
                << "For 2-D models, the empty direction "
                << "must be z!" << abort(FatalError);
        }

        twoD = true;
    }
    else if (nDir != 3)
    {
        FatalErrorIn("checkTwoD()")
            << "Only implemented for 2-D and 3-D models!"
            << abort(FatalError);
    }

    return twoD;
}


void Foam::sparseMatrixTools::solveLinearSystemEigen
(
    const sparseMatrix& matrix,
    const vectorField& source,
    vectorField& solution,
    const bool twoD,
    const bool exportToMatlab,
    const bool debug
)
{
    // For now, we can directly use the Eigen direct solver to solve the
    // linear system

    // Define the number of degrees of freedom
    label nDof;
    if (twoD)
    {
        nDof = 2*solution.size();
    }
    else
    {
        nDof = 3*solution.size();
    }

    // Create Eigen matrix triplets to store coefficients
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(nDof);

    const sparseMatrixData& data = matrix.data();
    for
    (
        sparseMatrixData::const_iterator iter = data.begin();
        iter != data.end();
        ++iter
    )
    {
        const tensor& coeff = iter();

        if (twoD)
        {
            const label rowI = 2*iter.key()[0];
            const label colI = 2*iter.key()[1];

            coefficients.push_back(Eigen::Triplet<scalar>(rowI, colI, coeff.xx()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI, colI+1, coeff.xy()));

            coefficients.push_back(Eigen::Triplet<scalar>(rowI+1, colI, coeff.yx()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI+1, colI+1, coeff.yy()));
        }
        else // 3-D
        {
            const label rowI = 3*iter.key()[0];
            const label colI = 3*iter.key()[1];

            coefficients.push_back(Eigen::Triplet<scalar>(rowI, colI, coeff.xx()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI, colI+1, coeff.xy()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI, colI+2, coeff.xz()));

            coefficients.push_back(Eigen::Triplet<scalar>(rowI+1, colI, coeff.yx()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI+1, colI+1, coeff.yy()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI+1, colI+2, coeff.yz()));

            coefficients.push_back(Eigen::Triplet<scalar>(rowI+2, colI, coeff.zx()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI+2, colI+1, coeff.zy()));
            coefficients.push_back(Eigen::Triplet<scalar>(rowI+2, colI+2, coeff.zz()));
        }
    }

    // Create Eigen sparse matrix
    Eigen::SparseMatrix<scalar> A(nDof, nDof);

    // Insert triplets into the matrix
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Compressing matrix is meant to help performance
    A.makeCompressed();

    // Create source vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(nDof);
    {
        label index = 0;
        forAll(source, i)
        {
            b(index++) = source[i].x();
            b(index++) = source[i].y();

            if (!twoD)
            {
                b(index++) = source[i].z();
            }
        }
    }

    if (exportToMatlab)
    {
        Info<< "Exporting linear system to matlabSparseMatrix.txt and "
            << "matlabSource.txt" << endl;

        // Write matrix
        Eigen::saveMarket(A, "matlabSparseMatrix.txt");

        // Write source
        OFstream sourceFile("matlabSource.txt");
        for (int rowI = 0; rowI < A.rows(); rowI++)
        {
            sourceFile
                << b(rowI) << endl;
        }
    }

    // Construct the solver
    Eigen::SparseLU
    <
        Eigen::SparseMatrix<scalar>, Eigen::COLAMDOrdering<int>
    > solver(A);

    // Initialise the solution vector to zero
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(nDof);
    x.setZero();

    // Check initial residual
    const Eigen::Matrix<scalar, Eigen::Dynamic, 1> initResidual = A*x - b;

    // Exit early if the initial residual is small
    if (initResidual.squaredNorm() < 1e-12)
    {
        Info<< "    Linear solver initial residual is "
            << initResidual.squaredNorm() << ": exiting" << endl;
    }
    else
    {
        // Solve system
        x = solver.solve(b);
    }

    // Copy  to solution field
    {
        label index = 0;
        forAll(solution, i)
        {
            solution[i].x() = x(index++);
            solution[i].y() = x(index++);

            if (!twoD)
            {
                solution[i].z() = x(index++);
            }
        }
    }
}


Foam::BlockSolverPerformance<Foam::vector>
Foam::sparseMatrixTools::solveLinearSystemPETSc
(
    const sparseMatrix& matrix,
    const vectorField& source,
    vectorField& solution,
    const bool twoD,
    fileName& optionsFile,
    const pointField& points,
    const boolList& ownedByThisProc,
    const labelList& localToGlobalPointMap,
    const labelList& stencilSizeOwned,
    const labelList& stencilSizeNotOwned,
    const bool debug
)
{
    if (debug)
    {
        Info<< "BlockSolverPerformance<vector> "
            << "sparseMatrixTools::solveLinearSystemPETSc: start" << endl;
    }

    // Set the block coefficient size (2 for 2-D, 3 for 3-D)
    label blockSize = 3;
    if (twoD)
    {
        blockSize = 2;
    }

    // Find size of global system, i.e. the highest global point index + 1
    const label blockN = gMax(localToGlobalPointMap) + 1;
    const label N = blockSize*blockN;
    if (debug)
    {
        Pout<< "blockN = " << blockN << ", N = " << N << endl;
    }

    // Find the start and end global point indices for this proc
    label blockStartID = N;
    label blockEndID = -1;
    forAll(ownedByThisProc, pI)
    {
        if (ownedByThisProc[pI])
        {
            blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
            blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
        }
    }
    const label startID = blockSize*blockStartID;
    const label endID = blockSize*blockEndID;
    if (debug)
    {
        Pout<< "blockStartID = " << blockStartID
            << ", blockEndID = " << blockEndID
            << ", startID = " << startID << ", endID = " << endID << endl;
    }

    // Find size of local system, i.e. the range of points owned by this proc
    const label blockn = blockEndID - blockStartID + 1;
    const label n = blockSize*blockn;
    if (debug)
    {
        Pout<< "blockn = " << blockn << ", n = " << n << endl;
    }

    // Define PETSc variables
    static char help[] = "Solves a linear system with KSP.\n\n";
    Vec            x, b;         /* approx solution, RHS */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    PetscReal      norm;         /* norm of solution error */
    PetscErrorCode ierr;
    PetscInt       its;
    PetscScalar*   xArr;
    MatInfo        matinfo;
    //PetscBool      nonzeroguess = PETSC_FALSE;

    // To keep PETSc happy, we will give it dummy command-line arguments
    int argc = 1;
    char * args1[argc];
    args1[0] = (char*)"solids4Foam";
    char** args = args1;

    // Initialise PETSc with options file
    optionsFile.expand();
    PetscInitialize(&argc, &args, optionsFile.c_str(), help);

    //MPI_Comm_size(PETSC_COMM_WORLD,&size);

    //ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    //ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
    //ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);


    // Create PETSc matrix

    MatCreate(PETSC_COMM_WORLD, &A);

    // Set the local and global matrix size
    // (matrix, local rows, local cols, global rows, global cols)
    //MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetSizes(A, n, n, N, N);

    MatSetFromOptions(A);

    // Set matrix to parallel type
    MatSetType(A, MATMPIAIJ);

    // Set the block coefficient size
    MatSetBlockSize(A, blockSize);

    // Pre-allocate matrix memory: this is critical for performance

    // Set on-core (d_nnz) and off-core (o_nnz) non-zeros per row
    // o_nnz is currently not set correctly, as it distinguishes between on-core
    // and off-core instead of owned (all on-core) vs not-owned (on-core and
    // off-core). For now, we will just use the max on-core non-zeros to
    // initialise not-owned values
    label d_nnz[n];
    label o_nnz[n];
    label d_nz = 0;
    setNonZerosPerRow
    (
        d_nnz,
        o_nnz,
        d_nz,
        n,
        blockSize,
        ownedByThisProc,
        stencilSizeOwned,
        stencilSizeNotOwned
    );

    // Find max non-zeros in a row
    if (debug)
    {
        Pout<< "        Max non-zeros per row = " << d_nz << endl;
    }

    // Serial matrix
    // Set exact number of non-zeros per row
    // MatSeqAIJSetPreallocation(A, 0, d_nnz);
    // or conservatively as
    // MatSeqAIJSetPreallocation(A, nz, NULL);

    // Parallel matrix
    MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz);
    // MatMPIAIJSetPreallocation(A, 0, d_nnz, d_nz, NULL);
    // or conservatively as
    // MatMPIAIJSetPreallocation(A, d_nz, NULL, o_nz, NULL);
    // const label nz = 81; // way too much in 2-D!
    // MatMPIAIJSetPreallocation(A, nz, NULL, nz, NULL);

    // Not sure if this set is needed but it does not hurt
    MatSetUp(A);


    // Insert coefficients into the matrix
    // Note: we use global indices when inserting coefficients

    if (debug)
    {
        Info<< "    Inserting PETSc matrix coefficients: start" << endl;
    }

    const sparseMatrixData& data = matrix.data();
    for
    (
        sparseMatrixData::const_iterator iter = data.begin();
        iter != data.end();
        ++iter
    )
    {
        const tensor& coeff = iter();
        const label blockRowI = localToGlobalPointMap[iter.key()[0]];
        const label blockColI = localToGlobalPointMap[iter.key()[1]];

        if (twoD)
        {
            // Prepare values
            const PetscScalar values[4] =
            {
                coeff.xx(), coeff.xy(),
                coeff.yx(), coeff.yy()
            };

            // Insert tensor coefficient
            MatSetValuesBlocked
            (
                A, 1, &blockRowI, 1, &blockColI, values, ADD_VALUES
            );
            // Pout<< "Mat insert: blockRow = " << blockRowI
            //     << ", row = " << blockRowI*blockSize << endl;
        }
        else // 3-D
        {
            // Prepare values
            // Maybe I can use coeff.cdata() here?
            const PetscScalar values[9] =
            {
                coeff.xx(), coeff.xy(), coeff.xz(),
                coeff.yx(), coeff.yy(), coeff.yz(),
                coeff.zx(), coeff.zy(), coeff.zz()
            };

            // Insert tensor coefficient
            ierr = MatSetValuesBlocked
            (
                A, 1, &blockRowI, 1, &blockColI, values, ADD_VALUES
            );

            if (ierr > 0)
            {
                FatalErrorIn("sparseMatrixTools::solveLinearSystemPETSc(...)")
                    << "PETSc returned the error code "<< ierr << nl
                    << "    blockRowI = " << blockRowI
                    << ", blockColI = " << blockColI
                    << ", coeff = " << coeff
                    << abort(FatalError);
            }
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // Check pre-allocation effectiveness, i.e. were additional memory
    // allocations needed or was space left unused
    if (debug)
    {
        //MatGetInfo(A, MAT_LOCAL, &matinfo);
        MatGetInfo(A, MAT_GLOBAL_SUM, &matinfo);
        // MatGetInfo(A, MAT_GLOBAL_MAX, &matinfo);
        Pout<< "nz_allocated = " << matinfo.nz_allocated
            << ", nz_used = " << matinfo.nz_used
            << ", nz_unneeded = " << matinfo.nz_unneeded
            << ", memory = " << matinfo.memory
            << ", assemblies = " << matinfo.assemblies
            << ", mallocs = " << matinfo.mallocs << endl;
    }

    if (debug)
    {
        Info<< "    Inserting PETSc matrix coefficients: end" << endl;
    }

    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // Populate the PETSc source vector
    // Note: uses global indices

    if (debug)
    {
        Info<< "        Populating the source vector" << endl;
    }

    // Create PETSc vectors
    // x and b have local size n
    // Note: the sizes of x and b are, in general, not equal the number of
    // points on this process, as they only contain the points owned by this
    // proc. To acess the values not-owned by the proc, we will use an
    // xWithGhosts vector
    // x is the global solution vector
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, n, N);
    VecSetBlockSize(x, blockSize);
    VecSetType(x, VECMPI);
    PetscObjectSetName((PetscObject) x, "Solution");
    // VecSetSizes(x, PETSC_DECIDE, N);
    VecSetFromOptions(x);

    // Create the source (b) using the same settings as b
    VecDuplicate(x, &b);
    PetscObjectSetName((PetscObject) b, "Source");

    {
        forAll(source, localBlockRowI)
        {
            const vector& sourceI = source[localBlockRowI];
            const label blockRowI = localToGlobalPointMap[localBlockRowI];

            if (twoD)
            {
                // Prepare values
                const PetscScalar values[2] = {sourceI.x(), sourceI.y()};

                // Insert values
                VecSetValuesBlocked(b, 1, &blockRowI, values, ADD_VALUES);
            }
            else
            {
                // Prepare values
                const PetscScalar values[3] =
                {
                    sourceI.x(), sourceI.y(), sourceI.z()
                };

                // Insert values
                VecSetValuesBlocked(b, 1, &blockRowI, values, ADD_VALUES);
            }
        }
    }

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);


    // Create KSP linear solver
    // Info<< "    Creating the linear solver" << endl;
    KSPCreate(PETSC_COMM_WORLD, &ksp);


    // Set operators. Here the matrix that defines the linear system
    // also serves as the preconditioning matrix.
    KSPSetOperators(ksp, A, A);


    // Set linear solver defaults for this problem
    // This are overwritten by the options file
    // - By extracting the KSP and PC contexts from the KSP context,
    //   we can then directly call any KSP and PC routines to set
    //   various options.
    // - The following four statements are optional; all of these
    //   parameters could alternatively be specified at runtime via
    //   KSPSetFromOptions();
    // ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    // ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    // ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc);
    //ierr = KSPSetType(ksp, KSPFGMRES);
    // ierr = PCSetType(pc, PCJACOBI);
    //ierr = PCSetType(pc, PCILU);
    ierr = KSPSetTolerances
    (
        ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT
    );

    // Set runtime options, e.g.,
    //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    // These options will override those specified above as long as
    // KSPSetFromOptions() is called _after_ any other customization
    // routines.
    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);

    // if (nonzeroguess)
    // {
    //     PetscScalar p = .5;
    //     // ierr = VecSet(x,p);CHKERRQ(ierr);
    //     // ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
    //     ierr = VecSet(x,p);
    //     ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    // }


    // Pass the point coordinates to PETSc to allow multigrid
    {
        PC pc;
        void (*f)(void) = NULL;

        KSPGetPC(ksp, &pc);
        PetscObjectQueryFunction((PetscObject)pc, "PCSetCoordinates_C", &f);

        if (f)
        {
            PetscInt sdim = vector::nComponents;
            if (twoD)
            {
                sdim = 2;
            }

            List<PetscReal> petscPoints(points.size()*sdim);

            auto iter = petscPoints.data();
            for (const vector& v : points)
            {
                *(iter++) = v.x();
                *(iter++) = v.y();

                if (!twoD)
                {
                    *(iter++) = v.z();
                }
            }

            PCSetCoordinates(pc, sdim, n, petscPoints.data());
        }
    }


    // Solve linear system
    if (debug)
    {
        Info<< "        Solving the linear solver: start" << endl;
    }

    ierr = KSPSolve(ksp, b, x);
    if (ierr > 0)
    {
        FatalErrorIn("sparseMatrixTools::solveLinearSystemPETSc(...)")
            << "PETSc linear solver returned the error code "<< ierr
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "        Solving the linear solver: end" << endl;
    }


    // Copy the results from the PETSc vector into the foam field

    if (debug)
    {
        Info<< "        Copying the solution vector" << endl;
    }

    // TODO: this only copies in the values owned by this proc
    // I still need to sync the points not owned by this proc
    VecGetArray(x, &xArr);
    {
        label index = 0;
        forAll(solution, i)
        {
            if (ownedByThisProc[i])
            {
                solution[i].x() = xArr[index++];
                solution[i].y() = xArr[index++];

                if (!twoD)
                {
                    solution[i].z() = xArr[index++];
                }
            }
        }
    }
    VecRestoreArray(x, &xArr);

    // Sync values not owned by this proc

    {
        // Count points non-owned by this proc
        label nNotOwnedByThisProc = 0;
        forAll(ownedByThisProc, pI)
        {
            if (!ownedByThisProc[pI])
            {
                nNotOwnedByThisProc++;
            }
        }
        nNotOwnedByThisProc *= blockSize;

        PetscInt indices[nNotOwnedByThisProc];
        int index = 0;
        forAll(ownedByThisProc, pI)
        {
            if (!ownedByThisProc[pI])
            {
                indices[index++] = blockSize*localToGlobalPointMap[pI];
                indices[index++] = blockSize*localToGlobalPointMap[pI] + 1;

                if (!twoD)
                {
                    indices[index++] = blockSize*localToGlobalPointMap[pI] + 2;
                }
            }
        }

        IS indexSet;
        ISCreateGeneral
        (
            PETSC_COMM_WORLD,
            nNotOwnedByThisProc,
            indices,
            PETSC_COPY_VALUES,
            &indexSet
        );

        // Local vector for holding not-owned values
        Vec xNotOwned;
        VecCreate(PETSC_COMM_WORLD, &xNotOwned);
        VecSetSizes(xNotOwned, nNotOwnedByThisProc, PETSC_DECIDE);
        VecSetType(xNotOwned, VECMPI);
        VecSetUp(xNotOwned);

        // Context for syncing data
        VecScatter ctx;
        VecScatterCreate(x, indexSet, xNotOwned, NULL, &ctx);
        VecScatterBegin(ctx, x, xNotOwned, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, x, xNotOwned, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&ctx);

        // Populate not-owned values
        PetscScalar* xNotOwnedArr;
        VecGetArray(xNotOwned, &xNotOwnedArr);
        {
            label index = 0;
            forAll(solution, i)
            {
                if (!ownedByThisProc[i])
                {
                    solution[i].x() = xNotOwnedArr[index++];
                    solution[i].y() = xNotOwnedArr[index++];

                    if (!twoD)
                    {
                        solution[i].z() = xNotOwnedArr[index++];
                    }
                }
            }
        }
        VecRestoreArray(xNotOwned, &xNotOwnedArr);

        // Destroy the index set
        ISDestroy(&indexSet);
    }

    // View solver info; we could instead use the option -ksp_view to
    // print this info to the screen at the conclusion of KSPSolve().
    if (debug)
    {
        ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        //CHKERRQ(ierr);
    }


    // Check the error
    // ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
    // ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    // ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
    //ierr = VecAXPY(x,neg_one,u);
    ierr = VecNorm(x, NORM_2, &norm);
    ierr = KSPGetIterationNumber(ksp, &its);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);
    // Info<< "    Error norm: " << norm << ", nIters = " << its << endl;


    // Free work space.  All PETSc objects should be destroyed when they
    // are no longer needed.
    // Info<< "        Freeing memory" << endl;
    // ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
    // ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
    // ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    //ierr = VecDestroy(&u);
    ierr = VecDestroy(&b);
    ierr = MatDestroy(&A);
    ierr = KSPDestroy(&ksp);


    // I should not call this here otherwise I cannot call this function again
    // Always call PetscFinalize() before exiting a program.  This routine
    //   - finalizes the PETSc libraries as well as MPI
    //   - provides summary and diagnostic information if certain runtime
    //     options are chosen (e.g., -log_summary).
    //ierr = PetscFinalize();

    if (debug)
    {
        Info<< "BlockSolverPerformance<vector> "
            << "sparseMatrixTools::solveLinearSystemPETSc: end" << endl;
    }

    vector initRes(vector::one);
    vector finalRes(norm*vector::one);

    if (twoD)
    {
        initRes.z() = 0;
        finalRes.z() = 0;
    }

    return BlockSolverPerformance<vector>
    (
        "PETSc", // solver name
        "pointD", // field name
        initRes, // initial residual
        finalRes, // final residual
        its // nIteration
        //false, // converged
        //false // singular
    );
}


void Foam::sparseMatrixTools::setNonZerosPerRow
(
    label d_nnz[],
    label o_nnz[],
    label& d_nz,
    const int nRows,
    const int blockSize,
    const boolList& ownedByThisProc,
    const labelList& stencilSizeOwned,
    const labelList& stencilSizeNotOwned
)
{
    // Initialise d_nnz and o_nnz to zero
    for (int i = 0; i < nRows; ++i)
    {
        d_nnz[i] = 0;
        o_nnz[i] = 0;
    }

    // Initialise max on-core non-zeros to zero
    d_nz = 0;

    label rowI = 0;
    forAll(stencilSizeOwned, blockRowI)
    {
        if (ownedByThisProc[blockRowI])
        {
            const label nCompOwned = blockSize*stencilSizeOwned[blockRowI];
            const label nCompNotOwned =
                blockSize*stencilSizeNotOwned[blockRowI];

            d_nz = max(d_nz, nCompOwned);

            d_nnz[rowI] += nCompOwned;
            o_nnz[rowI++] += nCompNotOwned;

            d_nnz[rowI] += nCompOwned;
            o_nnz[rowI++] += nCompNotOwned;

            if (blockSize == 3)
            {
                d_nnz[rowI] += nCompOwned;
                o_nnz[rowI++] += nCompNotOwned;
            }
            else if (blockSize != 2)
            {
                FatalErrorIn
                (
                    "void Foam::sparseMatrixTools::setNonZerosPerRow(...)"
                )   << "Not implemented for blockSize = " << blockSize
                    << abort(FatalError);
            }
        }
    }
}


void Foam::sparseMatrixTools::enforceFixedDof
(
    sparseMatrix& matrix,
    vectorField& source,
    const boolList& fixedDofs,
    const symmTensorField& fixedDofDirections,
    const pointField& fixedDofValues,
    const scalar fixedDofScale
)
{
    // Loop though the matrix and overwrite the coefficients for fixed DOFs
    // To enforce the value we will set the diagonal to the identity and set
    // the source to zero. The reason the source is zero is that we are solving
    // for the correction and the correction is zero for fixed values.
    // Rather than setting the identity on the diagonal, we will scale it by
    // fixedDofScale to improve the condition number, although the
    // preconditioner should not care.
    // Secondly, for any non-fixed-DOF equations which refer to fixed DOFs, we
    // will eliminate these coeffs and add their contribution (which is known)
    // to the source.
    sparseMatrixData& data = matrix.data();
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        const label blockRowI = iter.key()[0];
        const label blockColI = iter.key()[1];

        if (fixedDofs[blockRowI])
        {
            tensor& coeff = iter();

            // Set the source to zero as the correction to the displacement
            // is zero
            source[blockRowI] =
                ((I - fixedDofDirections[blockRowI]) & source[blockRowI]);

            // Eliminate the fixed directions from the coeff
            coeff = ((I - fixedDofDirections[blockRowI]) & coeff);

            if (blockRowI == blockColI)
            {
                // Scale the fixed directions to the 1 with scaling
                coeff -= tensor(fixedDofScale*fixedDofDirections[blockRowI]);
            }
        }
        else if (fixedDofs[blockColI])
        {
            // This equation refers to a fixed direction
            // We will eliminate the coeff and add the contribution to the
            // source
            tensor& coeff = iter();

            // // Add contribution to the source: this is not correct!
            // source[blockRowI] -=
            //     ((fixedDofDirections[blockColI] & coeff) & fixedDofValues[blockColI]);

            // Eliminate the fixed directions
            coeff = ((I - fixedDofDirections[blockColI]) & coeff);
        }
    }
}


// void Foam::sparseMatrixTools::addFixedDofToSource
// (
//     vectorField& source,
//     const boolList& fixedDofs,
//     const symmTensorField& fixedDofDirections,
//     const scalar fixedDofScale
// )
// {
//     forAll(fixedDofs, pointI)
//     {
//         // For now, we will only allow the full vector to be set
//         // To-do: allow fixed components
//         if (fixedDofs[pointI])
//         {
//             // Set the source to zero as the correction to the displacement
//             // is zero
//             source[pointI] = ((I - fixedDofDirections[pointI]) & source[pointI]);
//         }
//     }
// }


// void Foam::sparseMatrixTools::addFixedDofToMatrix
// (
//     sparseMatrix& matrix,
//     const boolList& fixedDofs,
//     const symmTensorField& fixedDofDirections,
//     const scalar fixedDofScale
// )
// {
//     // Loop though the matrix and overwrite the coefficients for fixed DOFs
//     // To enforce the value we will set the diagonal to the identity
//     // The preconditioner will be able to improve the conditioning for
//     // iterative solvers
//     sparseMatrixData& data = matrix.data();

//     for (auto iter = data.begin(); iter != data.end(); ++iter)
//     {
//         const label blockRowI = iter.key()[0];

//         if (fixedDofs[blockRowI])
//         {
//             const label blockColI = iter.key()[1];
//             tensor& coeff = iter();

//             if (blockRowI == blockColI)
//             {
//                 // Set the block diagonal to the identity matrix with scaling
//                 // coeff = -fixedDofScale*tensor(I);
//                 // Eliminate the fixed directions
//                 coeff = ((I - fixedDofDirections[blockRowI]) & coeff);
//                 // Set the fixed directions with scaling
//                 coeff = coeff - fixedDofScale*fixedDofDirections[blockRowI];
//             }
//             else
//             {
//                 // Eliminate the off-diagonals
//                 // coeff = tensor::zero;
//                 // Eliminate the fixed directions
//                 coeff = ((I - fixedDofDirections[blockRowI]) & coeff);
//             }
//         }
//     }
// }

// ************************************************************************* //
