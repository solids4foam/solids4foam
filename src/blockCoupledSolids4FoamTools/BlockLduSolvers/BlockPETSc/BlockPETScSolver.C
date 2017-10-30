/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "blockLduSolvers.H"
#include "BlockPETScSolver.H"
#include "BlockSolverPerformance.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "OFstream.H"
#include "polyMesh.H"

// PETSc header files
#include <petscksp.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockPETScSolver, 0);

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockPETScSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockPETScSolver, asymMatrix
    );

    // For PETSc
    static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockPETScSolver::BlockPETScSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector>(fieldName, matrix, dict),
    PETScArgs_(dict.lookup("PETScArgs")),
    tol_(readScalar(dict.lookup("tolerance")))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector> Foam::BlockPETScSolver::solve
(
    Field<Foam::vector>& U,
    const Field<Foam::vector>& blockB
)
{
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector>"
            "Foam::BlockPETScSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "SparesLU direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<vector>& matrix = this->matrix_;

    // Matrix addressing
    const unallocLabelList& lowerAddr = matrix.mesh().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix.mesh().lduAddr().upperAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor>& d = matrix.diag().asSquare();
    const Field<tensor>& l = matrix.lower().asSquare();
    const Field<tensor>& u = matrix.upper().asSquare();

    // Check if mesh is 2-D or 3-D
    // We lookup the mesh by name: we need a better way for multiple regions
    // models; I don't think it is possible from just examining the matrix and
    // source as the discretisation is performed in 3-D
    // We could set diag empty components to zero in the solver and then check
    // here... maybe
    const Vector<label>& geomD =
        matrix.mesh().thisDb().parent().lookupObject<polyMesh>
        (
            "region0" // polyMesh::defaultRegion,
        ).geometricD();

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
        Warning
            << this->typeName << ": for 2-D models, it is assumed that z is the"
            << " empty direction" << endl;
        twoD = true;
    }
    else if (nDir != 3)
    {
        FatalErrorIn(this->typeName + " solve")
            << "solver only implemented for 2-D and 3-D models"
            << abort(FatalError);
    }

    // Degrees of freedom
    int m = 0;
    if (twoD)
    {
        // 2-D
        m = 2*d.size();
    }
    else
    {
        // 3-D
        m = 3*d.size();
    }


    // Define PETSc variables
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": initialising PETSc"
            << endl;
    }
    //Vec            x, b, u;      /* approx solution, RHS, exact solution */
    Vec            x, b;      /* approx solution, RHS */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    PetscReal      norm;         /* norm of solution error */
    PetscErrorCode ierr;
    //PetscInt       id,n = 10,col[3],its;
    PetscInt       n = m,its;
    PetscMPIInt    size;
    //PetscScalar    neg_one      = -1.0,one = 1.0,value[3];
    //PetscScalar    value[1];
    PetscScalar *testVec;
    PetscBool      nonzeroguess = PETSC_FALSE;

    // To keep PETSc happy, we will give it dummy command-line arguments
    //int argc = 1;
    //char* args1 = NULL;
    //char** args = &args1;
    Info<< "Creating args" << endl;
    // Convert PETScArgs to argc and args char arrays
    int argc = PETScArgs_.size() + 1;
    //char** args;
    //word solverName = "solidFoam";
    // List< char* > PETScArgsChar(PETScArgs_.size());
    // //args[0] = const_cast<char*>(solverName.c_str());
    // PETScArgsChar[0] = new char[solverName.length() + 1];
    // std::strcpy(PETScArgsChar[0], solverName.c_str());
    // args[0] = PETScArgsChar[0];
    // for (int i = 1; i < argc; i++)
    // {
    //     //args[i] = const_cast<char*>(PETScArgs_[i - 1].c_str());
    //     //PETScArgsChar[i] = new char[PETScArgs_[i - 1].length() + 1];
    //     //std::strcpy(PETScArgsChar[i], PETScArgs_[i - 1].c_str());
    //     //args[i] = PETScArgsChar[i];
    //     //char* temp =  new char(PETScArgs_[i - 1].c_str());
    //     //args[i] = temp;
    // }
    //char* args1 = "-ksp_type";
    //char* args2 = "MY_TYPE";
    //args[0] = args0;
    //args[1] = args1;
    //args[2] = args2;
    // char ** argsTest = petscArgs;
    // char * args1[] = {
    //     "solidFoam",
    //     "-ksp_type",
    //     "bicg",
    // };
    char * args1[argc];
    args1[0] = "solidFoam";
    for (int i = 1; i < argc; i++)
    {
        char *cstr = new char[PETScArgs_[i - 1].length() + 1];
        strcpy(cstr, PETScArgs_[i - 1].c_str());
        args1[i] = cstr;
    }
    //args1[1] = "-ksp_type";
    //args1[2] = "gmres";
    char** args = args1;
    //char* args[] = {"solidFoam", "-ksp_type", "MYTYPE"};
    Info<< "print args: " << endl;
    for (int i = 0; i < argc; i++)
    {
        Info<< "argv: " << args[i] << endl;
    }
    PetscInitialize(&argc, &args, (char*)0, help);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    // Hmnn the compiler gives out about these macros, they seem related to error
    // catching; it may be OK without them
    //if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
    //CHKERRQ(ierr);
    //if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    //ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
    //ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
    ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);

    // Create vectors, where we create one vector and then duplicate it as needed
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": creating the PETSc vectors"
            << endl;
    }
    VecCreate(PETSC_COMM_WORLD,&x);
    PetscObjectSetName((PetscObject) x, "Solution");
    VecSetSizes(x,PETSC_DECIDE,n);
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
    //VecDuplicate(x,&u);

    // Create the matrix
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": creating the PETSc matrix"
            << endl;
    }
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetFromOptions(A);
    MatSetUp(A);

    // Assemble the matrix
    //value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
    //for (i=1; i<n-1; i++)
    //{
    //    col[0] = i-1; col[1] = i; col[2] = i+1;
    //    MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);
    //}
    //i    = n - 1; col[0] = n - 2; col[1] = n - 1;
    //MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    //i    = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
    //MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    // MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": assembling the matrix..."
            << endl;
    }

    // Insert diagonal
    label ID = -1;
    label ID_p1 = -1;
    label ID_p2 = -1;
    forAll(d, dI)
    {
        if (twoD)
        {
            // 2-D
            ID = 2*dI;
        }
        else
        {
            // 3-D
            ID = 3*dI;
        }

        ID_p1 = ID + 1;
        ID_p2 = ID + 2;

        //MatSetValues(A,1,&id,1,&id,value,INSERT_VALUES);
        MatSetValues(A, 1, &ID,    1, &ID,    &d[dI].xx(), INSERT_VALUES);
        MatSetValues(A, 1, &ID,    1, &ID_p1, &d[dI].xy(), INSERT_VALUES);
        MatSetValues(A, 1, &ID_p1, 1, &ID,    &d[dI].yx(), INSERT_VALUES);
        MatSetValues(A, 1, &ID_p1, 1, &ID_p1, &d[dI].yy(), INSERT_VALUES);

        if (!twoD)
        {
            MatSetValues(A, 1, &ID,    1, &ID_p2, &d[dI].xz(), INSERT_VALUES);
            MatSetValues(A, 1, &ID_p1, 1, &ID_p2, &d[dI].yz(), INSERT_VALUES);
            MatSetValues(A, 1, &ID_p2, 1, &ID,    &d[dI].zx(), INSERT_VALUES);
            MatSetValues(A, 1, &ID_p2, 1, &ID_p1, &d[dI].zy(), INSERT_VALUES);
            MatSetValues(A, 1, &ID_p2, 1, &ID_p2, &d[dI].zz(), INSERT_VALUES);
        }
    }

    // Insert off-diagonal
    label rowI = -1;
    label rowI_p1 = -1;
    label rowI_p2 = -1;
    label colI = -1;
    label colI_p1 = -1;
    label colI_p2 = -1;
    forAll(u, uI)
    {
        const label own = lowerAddr[uI];
        const label nei = upperAddr[uI];

        const tensor& upper = u[uI];
        const tensor& lower = l[uI];

        if (twoD)
        {
            rowI = 2*own;
            colI = 2*nei;
        }
        else
        {
            rowI = 3*own;
            colI = 3*nei;
        }

        rowI_p1 = rowI + 1;
        rowI_p2 = rowI + 2;
        colI_p1 = colI + 1;
        colI_p2 = colI + 2;

        // Upper
        MatSetValues(A, 1, &rowI,    1, &colI,    &upper.xx(), INSERT_VALUES);
        MatSetValues(A, 1, &rowI,    1, &colI_p1, &upper.xy(), INSERT_VALUES);
        MatSetValues(A, 1, &rowI_p1, 1, &colI,    &upper.yx(), INSERT_VALUES);
        MatSetValues(A, 1, &rowI_p1, 1, &colI_p1, &upper.yy(), INSERT_VALUES);

        // Lower
        MatSetValues(A, 1, &colI,    1, &rowI,    &lower.xx(), INSERT_VALUES);
        MatSetValues(A, 1, &colI,    1, &rowI_p1, &lower.xy(), INSERT_VALUES);
        MatSetValues(A, 1, &colI_p1, 1, &rowI,    &lower.yx(), INSERT_VALUES);
        MatSetValues(A, 1, &colI_p1, 1, &rowI_p1, &lower.yy(), INSERT_VALUES);

        if (!twoD)
        {
            // Upper
            MatSetValues(A, 1, &rowI,    1, &colI_p2, &upper.xz(), INSERT_VALUES);
            MatSetValues(A, 1, &rowI_p1, 1, &colI_p2, &upper.yz(), INSERT_VALUES);
            MatSetValues(A, 1, &rowI_p2, 1, &colI,    &upper.zx(), INSERT_VALUES);
            MatSetValues(A, 1, &rowI_p2, 1, &colI_p1, &upper.zy(), INSERT_VALUES);
            MatSetValues(A, 1, &rowI_p2, 1, &colI_p2, &upper.zz(), INSERT_VALUES);

            // Lower
            MatSetValues(A, 1, &colI,    1, &rowI_p2, &lower.xz(), INSERT_VALUES);
            MatSetValues(A, 1, &colI_p1, 1, &rowI_p2, &lower.yz(), INSERT_VALUES);
            MatSetValues(A, 1, &colI_p2, 1, &rowI,    &lower.zx(), INSERT_VALUES);
            MatSetValues(A, 1, &colI_p2, 1, &rowI_p1, &lower.zy(), INSERT_VALUES);
            MatSetValues(A, 1, &colI_p2, 1, &rowI_p2, &lower.zz(), INSERT_VALUES);
        }
    }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": assembling the matrix -> done"
            << endl;
    }

    // Copy the source from foam to PETSc
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": assembling the source..."
            << endl;
    }
    label index = 0;
    forAll(blockB, rowI)
    {
        VecSetValues(b, 1, &index, &blockB[rowI].x(), INSERT_VALUES);
        index++;
        VecSetValues(b, 1, &index, &blockB[rowI].y(), INSERT_VALUES);
        index++;

        if (!twoD)
        {
            VecSetValues(b, 1, &index, &blockB[rowI].z(), INSERT_VALUES);
            index++;
        }
    }

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": assembling the source -> done"
            << endl;
    }

    // Set exact solution and compute RHS
    //VecSet(u,one);
    //MatMult(A,u,b);

    // Create linear solver
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": creating the linear solver"
            << endl;
    }
    KSPCreate(PETSC_COMM_WORLD,&ksp);

    // Set operators. Here the matrix that defines the linear system
    // also serves as the preconditioning matrix.
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": setting the preconditioner"
            << endl;
    }
    KSPSetOperators(ksp,A,A);

    /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
        */
    // ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    // ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    // ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": setting linear solver options"
            << endl;
    }

    ierr = KSPGetPC(ksp,&pc);
    ierr = PCSetType(pc,PCJACOBI);
    ierr = KSPSetTolerances(ksp,tol_,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);


    /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
        */
    //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);

    //KSPGMRESSetRestart(ksp,500);

    if (nonzeroguess) {
        PetscScalar p = .5;
        // ierr = VecSet(x,p);CHKERRQ(ierr);
        // ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
        ierr = VecSet(x,p);
        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    }


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Solve linear system
        */
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": solving the linear system"
            << endl;
    }

    //ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,x);



    // Copy the results from the PETSc vector into the foam field
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": copying results into foam format"
            << endl;
    }

    VecGetArray(x, &testVec);
    index = 0;
    forAll(U, cellI)
    {
        U[cellI].x() = testVec[index++];
        U[cellI].y() = testVec[index++];

        if (!twoD)
        {
            U[cellI].z() = testVec[index++];
        }
    }

    VecRestoreArray(x, &testVec);


    /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
        */
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": checking solver info"
            << endl;
    }

    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Check the error
        */
    // ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
    // ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    // ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
    //ierr = VecAXPY(x,neg_one,u);
    ierr = VecNorm(x,NORM_2,&norm);
    ierr = KSPGetIterationNumber(ksp,&its);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);

    /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
        */
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": cleaning up"
            << endl;
    }

    // ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
    // ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
    // ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    //ierr = VecDestroy(&u);
    ierr = VecDestroy(&b);
    ierr = MatDestroy(&A);
    ierr = KSPDestroy(&ksp);

    /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
        */
    ierr = PetscFinalize();

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": done"
            << endl;
    }


    // Return solver info
    Warning
        << "Solver info not returned to foam yet" << endl;
    return BlockSolverPerformance<Foam::vector>
    (
        this->typeName,
        this->fieldName(),
        pTraits<Foam::vector>::one,
        pTraits<Foam::vector>::zero,
        0,
        true,
        false
    );
}


// ************************************************************************* //
