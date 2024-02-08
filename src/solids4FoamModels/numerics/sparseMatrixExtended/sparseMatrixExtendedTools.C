/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "sparseMatrixExtendedTools.H"
#include "OFstream.H"
#ifdef USE_PETSC
    #include <petscksp.h>
#endif

// * * * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * * * //

bool Foam::sparseMatrixExtendedTools::checkTwoD(const polyMesh& mesh)
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

#ifdef USE_PETSC

#ifdef OPENFOAM_NOT_EXTEND
    Foam::SolverPerformance<Foam::vector>
#else
    Foam::BlockSolverPerformance<Foam::vector>
#endif

Foam::sparseMatrixExtendedTools::solveLinearSystemPETSc
(
    const sparseMatrixExtended& matrix,
    const Field<RectangularMatrix<scalar>>& source,
    Field<RectangularMatrix<scalar>>& solution,
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
            << "sparseMatrixExtendedTools::solveLinearSystemPETSc: start" << endl;
    }

    // Set the block coefficient size (3 for 2-D, 4 for 3-D)
    label blockSize = 4;
    if (twoD)
    {
        blockSize = 3;
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
    const label endID = blockSize*(blockEndID + 1) - 1;
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


    // To keep PETSc happy, we will give it dummy command-line arguments
    // int argc = 1;
    // char * args1[argc];
    // args1[0] = (char*)"solids4Foam";
    // char** args = args1;

    // Initialise PETSc with options file
    optionsFile.expand();
    // static char help[] = "Solves a linear system with KSP.\n\n";
    PetscErrorCode ierr;
    // ierr = PetscInitialize(&argc, &args, optionsFile.c_str(), help); checkErr(ierr);
    if (debug)
    {
        Pout<< "PetscInitialize: start" << endl;
    }
    ierr = PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL); checkErr(ierr);
    if (debug)
    {
        Pout<< "PetscInitialize: end" << endl;
    }


    //MPI_Comm_size(PETSC_COMM_WORLD,&size);

    //ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    //ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
    //ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);


    // Create PETSc matrix

    Mat A;
    ierr = MatCreate(PETSC_COMM_WORLD, &A); checkErr(ierr);

    // Set the local and global matrix size
    // (matrix, local rows, local cols, global rows, global cols)
    //MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    ierr = MatSetSizes(A, n, n, N, N); checkErr(ierr);
    //ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N); checkErr(ierr);

    ierr = MatSetFromOptions(A); checkErr(ierr);

    // Set matrix to parallel type
    ierr = MatSetType(A, MATMPIAIJ); checkErr(ierr);

    // Set the block coefficient size
    ierr = MatSetBlockSize(A, blockSize); checkErr(ierr);

    // Pre-allocate matrix memory: this is critical for performance

    // Set on-core (d_nnz) and off-core (o_nnz) non-zeros per row
    // o_nnz is currently not set correctly, as it distinguishes between on-core
    // and off-core instead of owned (all on-core) vs not-owned (on-core and
    // off-core). For now, we will just use the max on-core non-zeros to
    // initialise not-owned values

    int* d_nnz = (int*)malloc(n*sizeof(int));
    int* o_nnz = (int*)malloc(n*sizeof(int));
    // label d_nnz[n];
    // label o_nnz[n];
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
    ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz); checkErr(ierr);
    //ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, d_nz, NULL); checkErr(ierr);
    // or conservatively as
    // ierr = MatMPIAIJSetPreallocation(A, d_nz, NULL, o_nz, NULL);
    // ierr = MatMPIAIJSetPreallocation(A, d_nz, NULL, d_nz, NULL); checkErr(ierr);
    // const label nz = 81; // way too much in 2-D!
    // ierr = MatMPIAIJSetPreallocation(A, nz, NULL, nz, NULL);

    // Optional: no error if additional memory allocation is required
    // If false, then an error is thrown for additional allocations
    // If preallocation was correct (or conservative) then an error should never
    // be thrown
    // For now, we will disable this check in debug mode so we can see how many
    // mallocs were made
    if (debug)
    {
        MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // Not sure if this set is needed but it does not hurt
    ierr = MatSetUp(A); checkErr(ierr);


    // Insert coefficients into the matrix
    // Note: we use global indices when inserting coefficients

    if (debug)
    {
        Pout<< "    Inserting PETSc matrix coefficients: start" << endl;
    }

    const sparseMatrixExtendedData& data = matrix.data();
    for
    (
        sparseMatrixExtendedData::const_iterator iter = data.begin();
        iter != data.end();
        ++iter
    )
    {
        const RectangularMatrix<scalar>& coeff = iter();
        const label blockRowI = localToGlobalPointMap[iter.key()[0]];
        const label blockColI = localToGlobalPointMap[iter.key()[1]];

        if (twoD)
        {
            // Prepare values
            const PetscScalar values[9] =
            {
                coeff(1,1), coeff(1,2), coeff(1,3),
                coeff(2,1), coeff(2,2), coeff(2,3),
                coeff(3,1), coeff(3,2), coeff(3,3)
            };

            // Insert tensor coefficient
            ierr = MatSetValuesBlocked
            (
                A, 1, &blockRowI, 1, &blockColI, values, ADD_VALUES
            ); checkErr(ierr);
        }
        else // 3-D
        {
            // Prepare values
            // Maybe I can use coeff.cdata() here?
            const PetscScalar values[16] =
            {
                coeff(1,1), coeff(1,2), coeff(1,3), coeff(1,4),
                coeff(2,1), coeff(2,2), coeff(2,3), coeff(2,4),
                coeff(3,1), coeff(3,2), coeff(3,3), coeff(3,4),
                coeff(4,1), coeff(4,2), coeff(4,3), coeff(4,4)
                
            };

            // Insert tensor coefficient
            ierr = MatSetValuesBlocked
            (
                A, 1, &blockRowI, 1, &blockColI, values, ADD_VALUES
            );
            if (ierr > 0)
            {
                Pout<< "MatSetValuesBlocked returned ierr = " << ierr
                    << " for " << blockRowI << " " << blockColI << ": "
                    << coeff << endl;
            }
            checkErr(ierr);
        }
    }
    if (debug)
    {
        Pout<< "    Inserting PETSc matrix coefficients: end" << endl;
    }

    if (debug)
    {
        Pout<< "        Assembling the matrix: start" << endl;
    }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); checkErr(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); checkErr(ierr);
    if (debug)
    {
        Pout<< "        Assembling the matrix: end" << endl;
    }

    // Check pre-allocation effectiveness, i.e. were additional memory
    // allocations needed or was space left unused
    if (debug)
    {
        MatInfo        matinfo;
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

    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // Populate the PETSc source vector
    // Note: uses global indices

    if (debug)
    {
        Pout<< "        Populating the source vector" << endl;
    }

    // Create PETSc vectors
    // x and b have local size n
    // Note: the sizes of x and b are, in general, not equal the number of
    // points on this process, as they only contain the points owned by this
    // proc. To acess the values not-owned by the proc, we will use an
    // xWithGhosts vector
    // x is the global solution vector
    Vec x;
    ierr = VecCreate(PETSC_COMM_WORLD, &x); checkErr(ierr);
    ierr = VecSetSizes(x, n, N); checkErr(ierr);
    ierr = VecSetBlockSize(x, blockSize); checkErr(ierr);
    ierr = VecSetType(x, VECMPI); checkErr(ierr);
    ierr =  PetscObjectSetName((PetscObject) x, "Solution"); checkErr(ierr);
    // VecSetSizes(x, PETSC_DECIDE, N);
    ierr = VecSetFromOptions(x); checkErr(ierr);

    // Create the source (b) using the same settings as b
    Vec b;
    ierr =  VecDuplicate(x, &b); checkErr(ierr);
    ierr = PetscObjectSetName((PetscObject) b, "Source"); checkErr(ierr);

    {
        forAll(source, localBlockRowI)
        {
            const RectangularMatrix<scalar>& sourceI = source[localBlockRowI];
            const label blockRowI = localToGlobalPointMap[localBlockRowI];

            if (twoD)
            {
                // Prepare values
                const PetscScalar values[3] = {sourceI(1,1), sourceI(2,1), sourceI(3,1)};

                // Insert values
                ierr = VecSetValuesBlocked
                (
                    b, 1, &blockRowI, values, ADD_VALUES
                ); checkErr(ierr);
            }
            else
            {
                // Prepare values
                const PetscScalar values[4] =
                {
                    sourceI(1,1), sourceI(2,1), sourceI(3,1), sourceI(4,1)
                };

                // Insert values
                ierr = VecSetValuesBlocked
                (
                    b, 1, &blockRowI, values, ADD_VALUES
                ); checkErr(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(b); checkErr(ierr);
    ierr = VecAssemblyEnd(b); checkErr(ierr);


    // Create KSP linear solver
    // Pout<< "    Creating the linear solver" << endl;
    KSP            ksp;          /* linear solver context */
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); checkErr(ierr);


    // Set operators. Here the matrix that defines the linear system
    // also serves as the preconditioning matrix.
    ierr = KSPSetOperators(ksp, A, A); checkErr(ierr);


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
    PC pc;
    ierr = KSPGetPC(ksp, &pc); checkErr(ierr);
    //ierr = KSPSetType(ksp, KSPFGMRES);
    // ierr = PCSetType(pc, PCJACOBI);
    //ierr = PCSetType(pc, PCILU);
    ierr = KSPSetTolerances
    (
        ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT
    ); checkErr(ierr);

    // Set runtime options, e.g.,
    //     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    // These options will override those specified above as long as
    // KSPSetFromOptions() is called _after_ any other customization
    // routines.
    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); checkErr(ierr);

    // PetscBool      nonzeroguess = PETSC_FALSE;
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

        ierr = KSPGetPC(ksp, &pc); checkErr(ierr);
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

            ierr = PCSetCoordinates(pc, sdim, n, petscPoints.data());
            checkErr(ierr);
        }
    }

    // Solve linear system
    if (debug)
    {
        Pout<< "        Solving the linear solver: start" << endl;
    }
	
    ierr = KSPSolve(ksp, b, x); checkErr(ierr);

    if (debug)
    {
        Pout<< "        Solving the linear solver: end" << endl;
    }


    // Copy the results from the PETSc vector into the foam field

    if (debug)
    {
        Pout<< "        Copying the solution vector" << endl;
    }

    // Insert local process values
    PetscScalar*   xArr;
    ierr = VecGetArray(x, &xArr); checkErr(ierr);
    {
        label index = 0;
        forAll(solution, i)
        {
            if (ownedByThisProc[i])
            {
            	if (twoD)
                {
		            solution[i](1,1) = xArr[index++];
		            solution[i](2,1) = xArr[index++];
		            solution[i](3,1) = xArr[index++];
                }
                else
                {
		            solution[i](1,1) = xArr[index++];
		            solution[i](2,1) = xArr[index++];
		            solution[i](3,1) = xArr[index++];
		            solution[i](4,1) = xArr[index++];
		        }
            }
        }
    }
    ierr = VecRestoreArray(x, &xArr); checkErr(ierr);

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
        ierr = ISCreateGeneral
        (
            PETSC_COMM_WORLD,
            nNotOwnedByThisProc,
            indices,
            PETSC_COPY_VALUES,
            &indexSet
        ); checkErr(ierr);

        // Local vector for holding not-owned values
        Vec xNotOwned;
        ierr = VecCreate(PETSC_COMM_WORLD, &xNotOwned); checkErr(ierr);
        ierr = VecSetSizes(xNotOwned, nNotOwnedByThisProc, PETSC_DECIDE);
        checkErr(ierr);
        ierr = VecSetType(xNotOwned, VECMPI); checkErr(ierr);
        ierr = VecSetUp(xNotOwned); checkErr(ierr);

        // Context for syncing data
        VecScatter ctx;
        ierr = VecScatterCreate(x, indexSet, xNotOwned, NULL, &ctx);
        checkErr(ierr);
        ierr = VecScatterBegin(ctx, x, xNotOwned, INSERT_VALUES, SCATTER_FORWARD);
        checkErr(ierr);
        ierr = VecScatterEnd(ctx, x, xNotOwned, INSERT_VALUES, SCATTER_FORWARD);
        checkErr(ierr);
        ierr = VecScatterDestroy(&ctx); checkErr(ierr);

        // Populate not-owned values
        PetscScalar* xNotOwnedArr;
        ierr = VecGetArray(xNotOwned, &xNotOwnedArr); checkErr(ierr);
        {
            label index = 0;
            forAll(solution, i)
            {
                if (!ownedByThisProc[i])
                {
		        	if (twoD)
		            {
				        solution[i](1,1) = xNotOwnedArr[index++];
				        solution[i](2,1) = xNotOwnedArr[index++];
				        solution[i](3,1) = xNotOwnedArr[index++];
		            }
		            else
		            {
				        solution[i](1,1) = xNotOwnedArr[index++];
				        solution[i](2,1) = xNotOwnedArr[index++];
				        solution[i](3,1) = xNotOwnedArr[index++];
				        solution[i](4,1) = xNotOwnedArr[index++];
				    }
                }
            }
        }
        ierr = VecRestoreArray(xNotOwned, &xNotOwnedArr); checkErr(ierr);

        // Destroy the index set
        ierr = ISDestroy(&indexSet); checkErr(ierr);
    }

    // View solver info; we could instead use the option -ksp_view to
    // print this info to the screen at the conclusion of KSPSolve().
    if (debug)
    {
        ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); checkErr(ierr);
    }


    // ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
    // ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    // ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
    //ierr = VecAXPY(x,neg_one,u);
    PetscInt       its;
    PetscReal      norm;         /* norm of solution error */
    ierr = VecNorm(x, NORM_2, &norm); checkErr(ierr);
    ierr = KSPGetIterationNumber(ksp, &its); checkErr(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);
    // Pout<< "    Error norm: " << norm << ", nIters = " << its << endl;


    // Free work space.  All PETSc objects should be destroyed when they
    // are no longer needed.
    // Pout<< "        Freeing memory" << endl;
    // ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
    // ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
    // ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x); checkErr(ierr);
    //ierr = VecDestroy(&u);
    ierr = VecDestroy(&b); checkErr(ierr);
    ierr = MatDestroy(&A); checkErr(ierr);
    ierr = KSPDestroy(&ksp); checkErr(ierr);


    // I should not call this here otherwise I cannot call this function again
    // Always call PetscFinalize() before exiting a program.  This routine
    //   - finalizes the PETSc libraries as well as MPI
    //   - provides summary and diagnostic information if certain runtime
    //     options are chosen (e.g., -log_summary).
    //ierr = PetscFinalize();

    if (debug)
    {
        Pout<< "BlockSolverPerformance<vector> "
            << "sparseMatrixExtendedTools::solveLinearSystemPETSc: end" << endl;
    }

    vector initRes(vector::one);
    vector finalRes(norm*vector::one);

    if (twoD)
    {
        initRes.z() = 0;
        finalRes.z() = 0;
    }

#ifdef OPENFOAM_NOT_EXTEND
    return SolverPerformance<vector>
    (
        "PETSc", // solver name
        "pointD", // field name
        initRes, // initial residual
        finalRes, // final residual
        vector(its, its, its) // nIteration
        //false, // converged
        //false // singular
    );
#else
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
#endif
}
#endif


void Foam::sparseMatrixExtendedTools::setNonZerosPerRow
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
                    "void Foam::sparseMatrixExtendedTools::setNonZerosPerRow(...)"
                )   << "Not implemented for blockSize = " << blockSize
                    << abort(FatalError);
            }
        }
    }
}


void Foam::sparseMatrixExtendedTools::checkErr(const int ierr)
{
    if (ierr > 0)
    {
        FatalError
            << "PETSc returned the error code "<< ierr
            << abort(FatalError);
    }
}


void Foam::sparseMatrixExtendedTools::enforceFixedDof
(
    sparseMatrixExtended& matrix,
    Field<RectangularMatrix<scalar>>& source,
    const bool twoD,
    const boolList& fixedDofs,
    const symmTensorField& fixedDofDirections,
    const pointField& fixedDofValues,
    const scalar fixedDofScale
)
{
	const bool debug = 0;
	
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
    // This is only done for the displacement coefficients and source terms
    // of the momentum equation.
    sparseMatrixExtendedData& data = matrix.data();
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        const label blockRowI = iter.key()[0];
        const label blockColI = iter.key()[1];

        if (fixedDofs[blockRowI])
        {
            RectangularMatrix<scalar>& coeff = iter();
            
            // Extract the displacement coefficients of the momentum equation
            // on the left-hand side
            tensor momDispCoeff(tensor::zero);
            
            if (twoD)
            {
		        momDispCoeff.xx() = coeff(1,1);
		        momDispCoeff.xy() = coeff(1,2);
		        momDispCoeff.yx() = coeff(2,1);
		        momDispCoeff.yy() = coeff(2,2);
            }
            else
            {
                momDispCoeff.xx() = coeff(1,1);
		        momDispCoeff.xy() = coeff(1,2);
		        momDispCoeff.xz() = coeff(1,3);
		        momDispCoeff.yx() = coeff(2,1);
		        momDispCoeff.yy() = coeff(2,2);
		        momDispCoeff.yz() = coeff(2,3);
		        momDispCoeff.zx() = coeff(3,1);
		        momDispCoeff.zy() = coeff(3,2);
		        momDispCoeff.zz() = coeff(3,3);
		    }
		
            // Extract the source terms of the momentum equation
            vector sourceTerms(vector::zero);
            
            if (twoD)
            {
		        sourceTerms.x() = source[blockRowI](1,1); 
		        sourceTerms.y() = source[blockRowI](2,1);
            }
            else
            {
		        sourceTerms.x() = source[blockRowI](1,1); 
		        sourceTerms.y() = source[blockRowI](2,1);
		        sourceTerms.z() = source[blockRowI](3,1);
		    }     
            
            if (debug)
            {
                Info<< "blockRow fixed: " << blockRowI << nl
                    << "    row,col: " << blockRowI << "," << blockColI << nl
                    << "    fixedDir: " << fixedDofDirections[blockRowI] << nl
                    << "    coeff before: " << momDispCoeff << endl;
            }

            // Free direction
            const tensor freeDir(I - fixedDofDirections[blockRowI]);

            // Set the source to zero as the correction to the displacement
            // is zero
            //source[blockRowI] = (freeDir & source[blockRowI]);
            sourceTerms = (freeDir & sourceTerms); 

            // Eliminate the fixed directions from the coeff
            momDispCoeff = (freeDir & momDispCoeff);

            if (blockRowI == blockColI)
            {
                // Remove the fixed component from the free component equation
                momDispCoeff = (freeDir & momDispCoeff & freeDir);

                // Fixed direction
                const tensor& fixedDir = fixedDofDirections[blockRowI];

                // Set the fixed direction diagonal to enforce a zero correction
                momDispCoeff -= tensor(fixedDofScale*fixedDir);
            }
            
            //Insert the changed coefficients back into the matrix
            if (twoD)
            {
            	coeff(1,1) = momDispCoeff.xx();
            	coeff(1,2) = momDispCoeff.xy();
            	coeff(2,1) = momDispCoeff.yx();
            	coeff(2,2) = momDispCoeff.yy();
            }
            else
            {
            	coeff(1,1) = momDispCoeff.xx();
            	coeff(1,2) = momDispCoeff.xy();
            	coeff(1,3) = momDispCoeff.xz();
            	coeff(2,1) = momDispCoeff.yx();
            	coeff(2,2) = momDispCoeff.yy();
            	coeff(2,3) = momDispCoeff.yz();
            	coeff(3,1) = momDispCoeff.zx();
            	coeff(3,2) = momDispCoeff.zy();
            	coeff(3,3) = momDispCoeff.zz();
            }  
            
            //Insert the changed source terms back into the source
            if (twoD)
            {
            	source[blockRowI](1,1) = sourceTerms.x(); 
            	source[blockRowI](2,1) = sourceTerms.y(); 
            }
            else
            {
            	source[blockRowI](1,1) = sourceTerms.x(); 
            	source[blockRowI](2,1) = sourceTerms.y(); 
            	source[blockRowI](3,1) = sourceTerms.z(); 
            } 

            if (debug)
            {
                Info<< "    coeff after: " << momDispCoeff << nl << endl;
            }
             
        }
        else if (fixedDofs[blockColI])
        {
            // This equation refers to a fixed direction
            // We will eliminate the coeff and add the contribution to the
            // source
            RectangularMatrix<scalar>& coeff = iter();
            
            if (debug)
            {
                Info<< "blockCol fixed: " << blockColI << nl
                    << "    row,col: " << blockRowI << "," << blockColI << nl
                    << "    fixedDir: " << fixedDofDirections[blockColI] << nl
                    << "    coeff before: " << coeff << endl;
            }
            
            // Extract the displacement coefficients of the momentum equation
            // on the left-hand side
            tensor momDispCoeff(tensor::zero);
            
            if (twoD)
            {
		        momDispCoeff.xx() = coeff(1,1);
		        momDispCoeff.xy() = coeff(1,2);
		        momDispCoeff.yx() = coeff(2,1);
		        momDispCoeff.yy() = coeff(2,2);
            }
            else
            {
                momDispCoeff.xx() = coeff(1,1);
		        momDispCoeff.xy() = coeff(1,2);
		        momDispCoeff.xz() = coeff(1,3);
		        momDispCoeff.yx() = coeff(2,1);
		        momDispCoeff.yy() = coeff(2,2);
		        momDispCoeff.yz() = coeff(2,3);
		        momDispCoeff.zx() = coeff(3,1);
		        momDispCoeff.zy() = coeff(3,2);
		        momDispCoeff.zz() = coeff(3,3);
		    }

            // Directions where the DOFs are unknown
            const tensor freeDir(I - fixedDofDirections[blockColI]);

            // Eliminate the fixed directions    
            momDispCoeff = (momDispCoeff & freeDir);
            
            //Insert the changed coefficients back into the matrix
            if (twoD)
            {
            	coeff(1,1) = momDispCoeff.xx();
            	coeff(1,2) = momDispCoeff.xy();
            	coeff(2,1) = momDispCoeff.yx();
            	coeff(2,2) = momDispCoeff.yy();
            }
            else
            {
            	coeff(1,1) = momDispCoeff.xx();
            	coeff(1,2) = momDispCoeff.xy();
            	coeff(1,3) = momDispCoeff.xz();
            	coeff(2,1) = momDispCoeff.yx();
            	coeff(2,2) = momDispCoeff.yy();
            	coeff(2,3) = momDispCoeff.yz();
            	coeff(3,1) = momDispCoeff.zx();
            	coeff(3,2) = momDispCoeff.zy();
            	coeff(3,3) = momDispCoeff.zz();
            }

            if (debug)
            {
                Info<< "    coeff after: " << momDispCoeff << nl << endl;
            }            

        }
    }
}


// void Foam::sparseMatrixExtendedTools::addFixedDofToSource
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


// void Foam::sparseMatrixExtendedTools::addFixedDofToMatrix
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
//     sparseMatrixExtendedData& data = matrix.data();

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
