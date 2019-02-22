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

#include "error.H"
#include "blockLduPrecons.H"
#include "BlockEigenILUTPrecon.H"
#include "blockLduPrecons.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"

// Eigen files
#include <Eigen/IterativeLinearSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenILUTPrecon, 0);
    addToRunTimeSelectionTable
    (
        blockVectorPrecon, BlockEigenILUTPrecon, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::BlockEigenILUTPrecon::calcPrecon
(
    const BlockLduMatrix<vector>& matrix
)
{
    Info<< "void Foam::BlockEigenILUTPrecon::calcPrecon(matrix)" << endl;

    // Tolerance for the inclusion of coefficients
    const scalar tol = VSMALL;

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

    //bool twoD = false;
    if (nDir == 2)
    {
        Warning
            << this->typeName << ": for 2-D models, it is assumed that z is the"
            << " empty direction" << endl;
        twoD_ = true;
    }
    else if (nDir != 3)
    {
        FatalErrorIn("void Foam::BlockEigenILUTPrecon::calcPrecon")
            << "solver only implemented for 2-D and 3-D models"
            << abort(FatalError);
    }

    Info<< this->typeName << ": copying matrix coefficients into Eigen format"
        << endl;

    // Degrees of freedom
    int m = 0;
    if (twoD_)
    {
        // 2-D
        m = 2*d.size();
    }
    else
    {
        // 3-D
        m = 3*d.size();
    }

    // Create coefficient matrix: we must copy coeffs from ldu storage
    // This can be costly
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(m);

    // Insert coefficients

    // Insert diagonal
    label ID = -1;
    forAll(d, dI)
    {
        if (twoD_)
        {
            // 2-D
            ID = 2*dI;
        }
        else
        {
            // 3-D
            ID = 3*dI;
        }

        if (mag(d[dI].xx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID, ID, d[dI].xx())
                );
        }

        if (mag(d[dI].xy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID, ID+1, d[dI].xy())
                );
        }

        if (mag(d[dI].yx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID+1, ID, d[dI].yx())
                );
        }

        if (mag(d[dI].yy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(ID+1, ID+1, d[dI].yy())
                );
        }

        if (!twoD_)
        {
            if (mag(d[dI].xz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID, ID+2, d[dI].xz())
                    );
            }

            if (mag(d[dI].yz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+1, ID+2, d[dI].yz())
                    );
            }

            if (mag(d[dI].zx()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+2, ID, d[dI].zx())
                    );
            }

            if (mag(d[dI].zy()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+2, ID+1, d[dI].zy())
                    );
            }

            if (mag(d[dI].zz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(ID+2, ID+2, d[dI].zz())
                    );
            }
        }
    }

    // Insert off-diagonal
    label rowI = -1;
    label columnI = -1;
    forAll(u, uI)
    {
        const label own = lowerAddr[uI];
        const label nei = upperAddr[uI];

        const tensor& upper = u[uI];
        const tensor& lower = l[uI];

        if (twoD_)
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

        if (mag(upper.xx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI, columnI, upper.xx())
                );
        }

        if (mag(upper.xy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI, columnI+1, upper.xy())
                );
        }

        if (mag(upper.yx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI+1, columnI, upper.yx())
                );
        }

        if (mag(upper.yy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rowI+1, columnI+1, upper.yy())
                );
        }

        // Lower

        if (mag(lower.xx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI, rowI, lower.xx())
                );
        }

        if (mag(lower.xy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI, rowI+1, lower.xy())
                );
        }

        if (mag(lower.yx()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI+1, rowI, lower.yx())
                );
        }

        if (mag(lower.yy()) > tol)
        {
            coefficients.push_back
                (
                    Eigen::Triplet<scalar>(columnI+1, rowI+1, lower.yy())
                );
        }

        if (!twoD_)
        {
            // Upper

            if (mag(upper.xz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI, columnI+2, upper.xz())
                    );
            }

            if (mag(upper.yz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+1, columnI+2, upper.yz())
                    );
            }

            if (mag(upper.zx()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+2, columnI, upper.zx())
                    );
            }

            if (mag(upper.zy()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+2, columnI+1, upper.zy())
                    );
            }

            if (mag(upper.zz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(rowI+2, columnI+2, upper.zz())
                    );
            }

            // Lower

            if (mag(lower.xz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI, rowI+2, lower.xz())
                    );
            }

            if (mag(lower.yz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+1, rowI+2, lower.yz())
                    );
            }

            if (mag(lower.zx()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+2, rowI, lower.zx())
                    );
            }

            if (mag(lower.zy()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+2, rowI+1, lower.zy())
                    );
            }

            if (mag(lower.zz()) > tol)
            {
                coefficients.push_back
                    (
                        Eigen::Triplet<scalar>(columnI+2, rowI+2, lower.zz())
                    );
            }
        }
    }


    // Create Eigen sparse matrix and set coeffs
    Info<< "Creating Eigen::SparseMatrix" << endl;
    APtr_.set(new Eigen::SparseMatrix<scalar>(m, m));
    Eigen::SparseMatrix<scalar>& A = APtr_();

    Info<< "    setFromTriplets" << endl;
    A.setFromTriplets(coefficients.begin(), coefficients.end());


    // Compressing matrix is meant to help performance
    Info<< "    Compressing the matrix" << endl;
    A.makeCompressed();


    // Create the preconditioner
    Info<< "    Creating Eigen::IncompleteLUT" << endl;
    //preconPtr_.set(new Eigen::IncompleteLUT<scalar>(A));
    preconPtr_.set(new Eigen::IncompleteLUT<scalar>(A, VSMALL, fillFactor_));

    Info<< "void Foam::BlockEigenILUTPrecon::calcPrecon(matrix) done" << endl;
}


const Eigen::SparseMatrix<Foam::scalar>& Foam::BlockEigenILUTPrecon::A() const
{
    if (APtr_.empty())
    {
        FatalErrorIn
        (
            "const Eigen::SparseMatrix<scalar>& Foam::BlockEigenILUTPrecon::A()"
        )   << "A pointer not set: calcPrecon must be called during "
            << "construction" << abort(FatalError);
    }

    return APtr_();
}


const Eigen::IncompleteLUT<Foam::scalar>&
Foam::BlockEigenILUTPrecon::precon() const
{
    if (preconPtr_.empty())
    {
        FatalErrorIn
        (
            "const Eigen::IncompleteLUT<scalar>& "
            "Foam::BlockEigenILUTPrecon::precon() const"
        )   << "precon pointer not set: calcPrecon must be called during "
            << "construction" << abort(FatalError);
    }

    return preconPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BlockEigenILUTPrecon::BlockEigenILUTPrecon
(
    const BlockLduMatrix<vector>& matrix
)
:
    BlockLduPrecon<vector>(matrix),
    APtr_(NULL),
    preconPtr_(NULL),
    twoD_(false),
    fillFactor_(0)
{
    calcPrecon(matrix);
}


Foam::BlockEigenILUTPrecon::BlockEigenILUTPrecon
(
    const BlockLduMatrix<vector>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<vector>(matrix),
    APtr_(NULL),
    preconPtr_(NULL),
    twoD_(false),
    fillFactor_(readInt(dict.lookup("fillFactor")))
{
    Info<< "    fill factor: " << fillFactor_ << endl;

    calcPrecon(matrix);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BlockEigenILUTPrecon::~BlockEigenILUTPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BlockEigenILUTPrecon::precondition
(
    Field<vector>& blockX,
    const Field<vector>& blockB
) const
{
    // Copy source and solution vectors into std::vector

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(A().rows());
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(A().rows());

    label index = 0;

    forAll(blockB, rowI)
    {
        b(index) = blockB[rowI].x();
        x(index) = blockX[rowI].x();
        index++;

        b(index) = blockB[rowI].y();
        x(index) = blockX[rowI].y();
        index++;

        if (!twoD_)
        {
            b(index) = blockB[rowI].z();
            x(index) = blockX[rowI].z();
            index++;
        }
    }


    // Perform the preconditioning
    precon()._solve_impl(b, x);
    //x = precon().solve(b);

    // Copy the solution vector back into the foam format

    index = 0;

    forAll(blockX, cellI)
    {
        blockX[cellI].x() = x(index++);
        blockX[cellI].y() = x(index++);

        if (!twoD_)
        {
            blockX[cellI].z() = x(index++);
        }
    }
}


void Foam::BlockEigenILUTPrecon::preconditionT
(
    Field<vector>& blockX,
    const Field<vector>& blockB
) const
{
    FatalErrorIn("void Foam::BlockEigenILUTPrecon::preconditionT")
        << "not implemented" << abort(FatalError);
}


// ************************************************************************* //
