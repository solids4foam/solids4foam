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
#include "BlockMUMPS.H"
#include "BlockSolverPerformance.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "OFstream.H"
#include "polyMesh.H"

// MUMPS files
#include "dmumps_c.h"
#include <vector>

// Eigen files
//#include <Eigen/Sparse>
//#include <unsupported/Eigen/SparseExtra>
//#include <Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>
//#include "unsupported/Eigen/src/IterativeSolvers/Scaling.h"
//#include <Eigen/MetisSupport>
//#include "Eigen/src/MetisSupport/MetisSupport.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockMUMPS, 0);
    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockMUMPS, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVectorSolver, BlockMUMPS, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockMUMPS::BlockMUMPS
(
    const word& fieldName,
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector>(fieldName, matrix, dict),
    exportSystemMatlab_
    (
        dict.lookupOrDefault<Switch>("exportSystemMatlab", false)
    ),
    scale_(dict.lookupOrDefault<Switch>("scale", false))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector>
Foam::BlockMUMPS::solve
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
            "Foam::BlockMUMPS::solve"
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

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into Eigen format"
            << endl;
    }

    // Degrees of freedom
    int m = 0;
    int nNonZeros = 0;
    if (twoD)
    {
        // 2-D
        m = 2*d.size();
        nNonZeros = 4*d.size() + 4*u.size() + 4*l.size();
    }
    else
    {
        // 3-D
        m = 3*d.size();
        nNonZeros = 9*d.size() + 9*u.size() + 9*l.size();
    }



    Warning
        << "FIX include files for MUMPS: use environment variables"
        << endl;

    // Set up MUMPS parameters

    DMUMPS_STRUC_C id;

    id.job = -1;
    id.par = 1;

    // id.sym = 1; // for positive definite (faster)
    id.sym = 2; // for symmetric
    // id.sym = 0; // for unsymmetric

    id.comm_fortran = 1;
    //Timing t1;
    dmumps_c(&id);
    //t1.write("Initialization");

    // Streams
    id.icntl[0] = 6;
    id.icntl[1] = 6;
    id.icntl[2] = 6;
    id.icntl[3] = 3;

    // Ordering metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)
    id.icntl[6] = 5;
    List<int> irn;
    List<int> jcn;
    List<double> a;


    // Convert lduMatrix to MUMPS representation
    //Timing r1;
    id.n = m;             // number of rows
    id.nz = nNonZeros;    // number of non-zeros
    irn.resize(id.nz);    // i-th row number
    jcn.resize(id.nz);    // j-th column number
    a.resize(id.nz);      // values

    // Insert diagonal
    int ID = -1;
    int k = 0;
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

        irn[k] = ID;
        jcn[k] = ID;
        a[k] = d[dI].xx();
        ++k;

        irn[k] = ID;
        jcn[k] = ID + 1;
        a[k] = d[dI].xy();
        ++k;

        irn[k] = ID + 1;
        jcn[k] = ID;
        a[k] = d[dI].yx();
        ++k;

        irn[k] = ID + 1;
        jcn[k] = ID + 1;
        a[k] = d[dI].yy();
        ++k;

        if (!twoD)
        {
            irn[k] = ID;
            jcn[k] = ID + 2;
            a[k] = d[dI].xz();
            ++k;

            irn[k] = ID + 2;
            jcn[k] = ID;
            a[k] = d[dI].zx();
            ++k;

            irn[k] = ID + 2;
            jcn[k] = ID + 1;
            a[k] = d[dI].zy();
            ++k;

            irn[k] = ID + 2;
            jcn[k] = ID + 2;
            a[k] = d[dI].zz();
            ++k;
        }
    }

    // Insert off-diagonal
    int rowI = -1;
    int columnI = -1;
    forAll(u, uI)
    {
        const label own = lowerAddr[uI];
        const label nei = upperAddr[uI];

        const tensor& upper = u[uI];
        const tensor& lower = l[uI];

        if (twoD)
        {
            rowI = 2*own;
            columnI = 2*nei;
        }
        else
        {
            rowI = 3*own;
            columnI = 3*nei;
        }

        // Upper

        irn[k] = rowI;
        jcn[k] = columnI;
        a[k] = upper.xx();
        ++k;

        irn[k] = rowI;
        jcn[k] = columnI + 1;
        a[k] = upper.xy();
        ++k;

        irn[k] = rowI + 1;
        jcn[k] = columnI;
        a[k] = upper.yx();
        ++k;

        irn[k] = rowI + 1;
        jcn[k] = columnI + 1;
        a[k] = upper.yy();
        ++k;

        // Lower

        irn[k] = columnI;
        jcn[k] = rowI;
        a[k] = lower.xx();
        ++k;

        irn[k] = columnI;
        jcn[k] = rowI + 1;
        a[k] = lower.xy();
        ++k;

        irn[k] = columnI + 1;
        jcn[k] = rowI;
        a[k] = lower.yx();
        ++k;

        irn[k] = columnI + 1;
        jcn[k] = rowI + 1;
        a[k] = lower.yy();
        ++k;

        if (!twoD)
        {
            // Upper

            irn[k] = rowI;
            jcn[k] = columnI + 2;
            a[k] = upper.xz();
            ++k;

            irn[k] = rowI + 1;
            jcn[k] = columnI + 2;
            a[k] = upper.yz();
            ++k;

            irn[k] = rowI + 2;
            jcn[k] = columnI;
            a[k] = upper.zx();
            ++k;

            irn[k] = rowI + 2;
            jcn[k] = columnI + 1;
            a[k] = upper.zy();
            ++k;

            irn[k] = rowI + 2;
            jcn[k] = columnI + 2;
            a[k] = upper.zz();
            ++k;

            // Lower

            irn[k] = columnI;
            jcn[k] = rowI + 2;
            a[k] = lower.xz();
            ++k;

            irn[k] = columnI + 1;
            jcn[k] = rowI + 2;
            a[k] = lower.yz();
            ++k;

            irn[k] = columnI + 2;
            jcn[k] = rowI;
            a[k] = lower.zx();
            ++k;

            irn[k] = columnI + 2;
            jcn[k] = rowI + 1;
            a[k] = lower.zy();
            ++k;

            irn[k] = columnI + 2;
            jcn[k] = rowI + 2;
            a[k] = lower.zz();
            ++k;
        }
    }

    id.irn = &*irn.begin();
    id.jcn = &*jcn.begin();
    id.a = &*a.begin();

    //r1.write("Converting matrix");
    //Timing t2;


    // LU factorization

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": direct solution of sparse system" << endl;
    }

    id.job = 1;
    dmumps_c(&id);
    //t2.write("Analysis");
    //Timing t3;
    id.job = 2;
    dmumps_c(&id);
    //t3.write("Factorization");


    // Copy source (right-hand side) vector into List

    //std::vector<scalar> b(m);
    scalar b[m];

    label index = 0;
    forAll(blockB, rowI)
    {
        b[index++] = blockB[rowI].x();
        b[index++] = blockB[rowI].y();

        if (!twoD)
        {
            b[index++] = blockB[rowI].z();
        }
    }


    // Reading RHS and setting up MUMPS

    id.rhs = &*b;
    id.nrhs = 1;        // Number of right-hand sides
    id.lrhs = id.n;     // Dimension of matrix


    // Back substitution: this puts the solution into the right-hand side i.e. b
    id.job = 3;
    dmumps_c(&id);
    //t4.write("Back substitution");


    // Clean up
    //Timing t5;
    id.job = -2;
    dmumps_c(&id);
    //t5.write("Cleaning MUMPS");

    // if (BlockLduSolver::debug)
    // {
    //     Info<< "    solver time = "
    //         << (matrix.mesh().thisDb().time().elapsedCpuTime() - startT)
    //         << endl;
    // }

    // We copy the results from the std::vector into the geometric field
    // This can be costly

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": copying results into foam format" << nl
            << endl;
    }

    index = 0;
    forAll(U, cellI)
    {
        U[cellI].x() = b[index++];
        U[cellI].y() = b[index++];

        if (!twoD)
        {
            U[cellI].z() = b[index++];
        }
    }

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
