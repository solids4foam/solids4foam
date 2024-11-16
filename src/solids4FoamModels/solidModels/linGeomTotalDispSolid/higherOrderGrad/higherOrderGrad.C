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

#include "higherOrderGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "emptyPolyPatch.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(higherOrderGrad, 0);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void higherOrderGrad::makeStencils() const
{
    Info<< "higherOrderGrad::makeStencils()" << endl;

    if (stencilsPtr_ || cellBoundaryFacesPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const labelListList& cellCells = mesh.cellCells();

    stencilsPtr_.set(new List<DynamicList<label>>(mesh.nCells()));
    List<DynamicList<label>>& stencils = *stencilsPtr_;

    cellBoundaryFacesPtr_.set(new List<DynamicList<label>>(mesh.nCells()));
    List<DynamicList<label>>& cellBoundaryFaces = *cellBoundaryFacesPtr_;

    forAll(stencils, cellI)
    {
       DynamicList<label>& curStencil = stencils[cellI];
       curStencil.setCapacity(maxStencilSize_);
       const labelList& curCellCells = cellCells[cellI];

       labelHashSet stencilCells;
       labelHashSet prevLayer;

       // Add first layer of cells
       forAll(curCellCells, cI)
       {
           stencilCells.insert(curCellCells[cI]);
           prevLayer.insert(curCellCells[cI]);
       }

       // Remaining layers of cells
       for (int layerI = 1; layerI < nLayers_; layerI++)
       {
           labelList prevLayerCells(prevLayer.toc());
           labelHashSet curLayer;

           // Loop over previous layer and add one level of
           // layer neighbours
           for (const label cI : prevLayerCells)
           {
               const labelList& cellINei= mesh.cellCells()[cI];
               forAll (cellINei, nei)
               {
                   if (!stencilCells.found(cellINei[nei]))
                   {
                       curLayer.insert(cellINei[nei]);
                   }
               }
           }

           // Now we have curent layer which we need to add to the stencil
           // and the current layer will now be the previous layer for next
           // loop
           prevLayer.clear();
           prevLayer = curLayer;
           stencilCells.merge(curLayer);
       }

       if (Pstream::parRun())
       {
           notImplemented("not implemented for parallel run");
       }

       // Neighbours of first layer will store cellI in stencil so we
       // need to remove it
       stencilCells.erase(cellI);

       curStencil.append(stencilCells.toc());
    }

    forAll(stencils, cellI)
    {
       // Note: stencils are sorted but I do not see problem in that
       stencils[cellI].shrink();
    }


    // Make cell boundary face stencils
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    forAll(boundaryMesh, patchI)
    {
        if
        (
            includePatchInStencils_[patchI]
         && boundaryMesh[patchI].type() != emptyPolyPatch::typeName
         && !boundaryMesh[patchI].coupled()
        )
        {
            const labelUList& faceCells = boundaryMesh[patchI].faceCells();

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                const label gI = faceI + boundaryMesh[patchI].start();
                cellBoundaryFaces[cellID].append(gI);
            }
        }
    }

    forAll(cellBoundaryFaces, cellI)
    {
        cellBoundaryFaces[cellI].shrink();
    }

    Info<< "higherOrderGrad::makeStencils(): end" << endl;
}


const List<DynamicList<label>> higherOrderGrad::stencils() const
{
    if (!stencilsPtr_)
    {
        makeStencils();
    }

    return *stencilsPtr_;
}


void higherOrderGrad::generateExponents
(
    const label N,
    DynamicList<FixedList<label, 3>>& exponents
) const
{
    // Estimate the number of terms to set the capacity
    const label estimatedSize = (N + 1)*(N + 2)*(N + 3)/6;
    exponents.setCapacity(estimatedSize);

    // Add the constant term first
    exponents.append(FixedList<label, 3>{0, 0, 0});

    for (label n = 1; n <= N; ++n)
    {
        for (label i = n; i >= 0; --i)
        {
            for (label j = n - i; j >= 0; --j)
            {
                label k = n - i - j;
                if (i == 0 && j == 0 && k == 0)
                {
                    // Skip the constant term as it's already added
                    continue;
                }
                FixedList<label, 3> exponent = {i, j, k};
                exponents.append(exponent);
            }
        }
    }

    // Adjust capacity to actual size
    exponents.shrink();
}


void higherOrderGrad::calcCoeffs() const
{
    Info<< "higherOrderGrad::calcCoeffs()" << endl;

    if (interpCoeffsPtr_ || interpGradCoeffsPtr_)
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    interpCoeffsPtr_.set(new List<DynamicList<scalar>>(mesh.nCells()));
    List<DynamicList<scalar>>& interpCoeffs = *interpCoeffsPtr_;

    interpGradCoeffsPtr_.set(new List<DynamicList<vector>>(mesh.nCells()));
    List<DynamicList<vector>>& interpGradCoeffs = *interpGradCoeffsPtr_;

    if (!useQRDecomposition_)
    {
        choleskyPtr_.set(new List<Eigen::LLT<Eigen::MatrixXd>>(mesh.nCells()));
        QhatPtr_.set(new List<Eigen::MatrixXd>(mesh.nCells()));
    }

    sqrtWPtr_.set
    (
        new List<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>(mesh.nCells())
    );

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();
    const vectorField& CfI = mesh.Cf();

    // Calculate Taylor series exponents
    // 1 for zero order, 4 for 1 order, 10 for second order, etc.
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();
    if (debug)
    {
        Info<< "Np = " << Np << endl;
    }

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    if (calcConditionNumber_)
    {
        // Cells condition numbers
        condNumberPtr_.set
        (
            new volScalarField
            (
               IOobject
               (
                   "condNumber",
                   mesh.time().timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh,
               dimensionedScalar("0", dimless, Zero)
            )
        );
    }

    List<DynamicList<scalar>> c(mesh.nCells());
    List<DynamicList<scalar>> cx(mesh.nCells());
    List<DynamicList<scalar>> cy(mesh.nCells());
    List<DynamicList<scalar>> cz(mesh.nCells());

    const List<DynamicList<label>>& cellBoundaryFaces = this->cellBoundaryFaces();
    const List<DynamicList<label>> stencils = this->stencils();

    if (N_ < 1)
    {
        FatalErrorInFunction
            << "N must be at least 1!" << exit(FatalError);
    }

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiCellID = curStencil[cI];
            const scalar d = mag(CI[neiCellID] - CI[cellI]);
            if (d > maxDist)
            {
                maxDist = d;
            }
        }

        // Loop over neighbours and construct matrix Q
        const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

        // Use matrix format from Eigen/Dense library
        // Avoid initialisation to zero as we will set every entry below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < Nn; ++cI)
        {
            vector dx;
            if (cI < curStencil.size())
            {
                const label neiCellID = curStencil[cI];
                const vector& neiC = CI[neiCellID];
                dx = neiC - CI[cellI];
            }
            else
            {
                const label i = cI - curStencil.size();
                const label globalFaceID = cellBoundaryFaces[cellI][i];
                const vector& neiC = CfI[globalFaceID];
                dx = neiC - CI[cellI];
            }

            // Compute monomial values for each exponent
            for (label p = 0; p < Np; ++p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = exponent[2];

               // Compute factorial denominator
               const scalar factorialDenominator =
                   factorials[i]*factorials[j]*factorials[k];

               // Compute and assign monomial value with factorials
               // Note: the order of the quadratic and higher terms may not be
               // the same as the previous manual approach
               Q(p, cI) =
                   pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                  /factorialDenominator;
            }
        }

        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);
        //W.setZero(); // no need to waste time initialising

        for (label cI = 0; cI < Nn; cI++)
        {
            scalar d;

            if (cI < curStencil.size())
            {
                const vector& C = CI[cellI];
                const label neiCellID = curStencil[cI];
                const vector& neiC = CI[neiCellID];
                d = mag(neiC - C);
            }
            else
            {
                // For boundary cells we need to add boundary face as
                // neigbour
                const vector& C = CI[cellI];
                const label i = cI - curStencil.size();
                const label globalFaceID = cellBoundaryFaces[cellI][i];
                const vector& neiC = CfI[globalFaceID];
                d = mag(neiC - C);
            }

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W.diagonal()[cI] = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        Eigen::DiagonalMatrix<double, Eigen::Dynamic>& sqrtW =
            sqrtWPtr_()[cellI];
        sqrtW = W.diagonal().cwiseSqrt().asDiagonal();
        const Eigen::MatrixXd Qhat =
            Q.array().rowwise()*sqrtW.diagonal().transpose().array();

        // B hat
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
            sqrtW.diagonal().asDiagonal();

        // Declare A outside the if-else scope, but do not initialise
        Eigen::MatrixXd A;

        if (useQRDecomposition_)
        {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

            // Q and R matrices
            const Eigen::MatrixXd O = qr.householderQ();
            const Eigen::MatrixXd& R = qr.matrixQR().triangularView<Eigen::Upper>();

            // Slice Rbar and Qbar, as we do not need full matrix
            // Note: auto is a reference type here (Rbar, Qbar are not copied)
            const auto Rbar = R.topLeftCorner(Np, Np);
            const auto Qbar = O.leftCols(Np);

            // Perform element-wise multiplication and convert to MatrixXd
            const Eigen::MatrixXd QbarBhat =
                (
                    Qbar.transpose().array().rowwise()
                   *Bhat.diagonal().transpose().array()
                ).matrix();

            // Solve to get A
            // const Eigen::MatrixXd A =
            //     Rbar.colPivHouseholderQr().solve(Qbar.transpose()*Bhat);
            // Solve using the modified QbarBhat
            A = Rbar.colPivHouseholderQr().solve(QbarBhat);

            // To be aware of interpolation accuracy we need to control the
            // condition number
            if (calcConditionNumber_)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd
                (
                    Rbar, Eigen::ComputeFullU | Eigen::ComputeFullV
                );
                Eigen::VectorXd singularValues = svd.singularValues();

                volScalarField& condNumber = condNumberPtr_();
                condNumber[cellI] =
                    singularValues(0)
                   /(singularValues(singularValues.size() - 1) + VSMALL);
            }

            c[cellI].setCapacity(A.cols());
            cx[cellI].setCapacity(A.cols());
            cy[cellI].setCapacity(A.cols());
            cz[cellI].setCapacity(A.cols());

            Eigen::RowVectorXd cRow = A.row(0);
            Eigen::RowVectorXd cxRow = A.row(1);
            Eigen::RowVectorXd cyRow = A.row(2);
            Eigen::RowVectorXd czRow = A.row(3);

            for (label i = 0; i < A.cols(); ++i)
            {
                c[cellI].append(cRow(i));
                cx[cellI].append(cxRow(i));
                cy[cellI].append(cyRow(i));
                cz[cellI].append(czRow(i));
            }

            c[cellI].shrink();
            cx[cellI].shrink();
            cy[cellI].shrink();
            cz[cellI].shrink();
        }
        else // Cholesky decomposition of the "normal equations"
        {
            // Transpose Q to follow the standard convention
            // TODO: avoid this by assigning Q correctly from the start!
            // It may be clear to seperate QR and Cholesky into their own
            // functions
            Q = Q.transpose().eval();

            // Compute Q_hat = Q * W^{1/2}
            Eigen::MatrixXd& Qhat = QhatPtr_()[cellI];
            Qhat = Q.array().colwise()*sqrtW.diagonal().array();

            // Compute N = Q_hat^T * Q_hat = Q^T W Q
            Eigen::MatrixXd N = Qhat.transpose()*Qhat;

            if (debug)
            {
                Eigen::FullPivLU<Eigen::MatrixXd> lu(Q);
                int rank = lu.rank();
                Info<< "Rank of Q: " << rank << nl
                    << "Q rows: " << Q.rows() << nl
                    << "Q cols: " << Q.cols() << endl;

                if (rank < Q.cols())
                {
                    Warning
                        << "Design matrix Q is rank-deficient!" << endl;
                }

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(N);
                if (eigensolver.info() != Eigen::Success)
                {
                    Warning
                        << "Eigenvalue computation failed!" << endl;
                    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
                    std::cout
                        << "Eigenvalues of N: " << eigenvalues.transpose()
                        << std::endl;
                }
            }

            // Perform Cholesky decomposition
            Eigen::LLT<Eigen::MatrixXd>& cholesky = choleskyPtr_()[cellI];
            cholesky.compute(N);
            if (cholesky.info() != Eigen::Success)
            {
                FatalErrorInFunction
                    << "Cholesky decomposition failed; "
                    << "matrix is not positive definite."
                    << exit(FatalError);
            }
        }
    }

    if (useQRDecomposition_)
    {
        forAll(interpCoeffs, cellI)
        {
           const DynamicList<label>& curStencil = stencils[cellI];
           const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

           interpCoeffs[cellI].setCapacity(Nn);
           interpGradCoeffs[cellI].setCapacity(Nn);

           for (label I = 0; I < Nn; I++)
           {
               interpCoeffs[cellI].append(c[cellI][I]);
               interpGradCoeffs[cellI].append
               (
                   vector(cx[cellI][I], cy[cellI][I], cz[cellI][I])
               );
           }

           interpCoeffs[cellI].shrink();
           interpGradCoeffs[cellI].shrink();
        }

        // We can clear sqrtW since it is no longer needed for the QR
        // decomposition approach
        sqrtWPtr_.clear();
    }

    Info<< "higherOrderGrad::calcCoeffs(): end" << endl;
}


const List<DynamicList<scalar>>& higherOrderGrad::interpCoeffs() const
{
    if (!interpCoeffsPtr_)
    {
        calcCoeffs();
    }

    return *interpCoeffsPtr_;
}


const List<DynamicList<vector>>& higherOrderGrad::interpGradCoeffs() const
{
    if (!interpGradCoeffsPtr_)
    {
        calcCoeffs();
    }

    return *interpGradCoeffsPtr_;
}


const List<DynamicList<label>>& higherOrderGrad::cellBoundaryFaces() const
{
    if (!cellBoundaryFacesPtr_)
    {
        makeStencils();
    }

    return *cellBoundaryFacesPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

higherOrderGrad::higherOrderGrad
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
:
    mesh_(mesh),
    includePatchInStencils_(includePatchInStencils),
    N_(readInt(dict.lookup("N"))),
    nLayers_(readInt(dict.lookup("nLayers"))),
    k_(readScalar(dict.lookup("k"))),
    maxStencilSize_(readInt(dict.lookup("maxStencilSize"))),
    useQRDecomposition_(dict.lookup("useQRDecomposition")),
    calcConditionNumber_(dict.lookup("calcConditionNumber")),
    condNumberPtr_(),
    stencilsPtr_(),
    cellBoundaryFacesPtr_(),
    interpCoeffsPtr_(),
    interpGradCoeffsPtr_(),
    choleskyPtr_(),
    QhatPtr_(),
    sqrtWPtr_()
{
    if (calcConditionNumber_)
    {
        if (!useQRDecomposition_)
        {
            FatalErrorInFunction
                << "useQRDecomposition must be 'on' when `calcConditionNumber` is 'on'"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

higherOrderGrad::~higherOrderGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volTensorField> higherOrderGrad::grad(const volVectorField& D)
{
    // Info<< "higherOrderGrad::grad(...)" << endl;

    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<volTensorField> tgradD
    (
        new volTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero),
           "zeroGradient"
        )
    );
    volTensorField& gradD = tgradD.ref();

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const List<DynamicList<label>> stencils = this->stencils();
    const List<DynamicList<label>> cellBoundaryFaces =
        this->cellBoundaryFaces();
    //const List<DynamicList<scalar>>& interpCoeffs = this->interpCoeffs();
    const List<DynamicList<vector>>& interpGradCoeffs =
        this->interpGradCoeffs();

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];
        const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

        if (useQRDecomposition_)
        {
            // Loop over stencil and multiply stencil cell values with
            // corresponding interpolation coefficient
            for(label cI = 0; cI < Nn; cI++)
            {
                if (cI < curStencil.size())
                {
                    gradD[cellI] += interpGradCoeffs[cellI][cI]*D[curStencil[cI]];
                }
                else
                {
                    const label i = cI - curStencil.size();
                    const label globalFaceID = cellBoundaryFaces[cellI][i];

                    vector boundaryD = vector::zero;

                    forAll(boundaryMesh, patchI)
                    {
                        if
                        (
                            includePatchInStencils_[patchI]
                         && boundaryMesh[patchI].type() != emptyPolyPatch::typeName
                         && !boundaryMesh[patchI].coupled()
                        )
                        {
                            const label start = boundaryMesh[patchI].start();
                            const label nFaces = boundaryMesh[patchI].nFaces();

                            if (globalFaceID >= start && globalFaceID < start + nFaces)
                            {
                                const label k = globalFaceID - start;
                                boundaryD = D.boundaryField()[patchI][k];
                            }
                        }
                    }

                    gradD[cellI] += interpGradCoeffs[cellI][cI]*boundaryD;
                }
            }
        }
        else // Cholesky decomposition
        {
            for (label cmptI = 0; cmptI < vector::nComponents; ++cmptI)
            {
                // Prepare right-hand side vector (y) for Cholesky decomposition
                Eigen::VectorXd y(Nn);

                for (label cI = 0; cI < Nn; cI++)
                {
                    if (cI < curStencil.size())
                    {
                        y[cI] = D[curStencil[cI]][cmptI];
                    }
                    else
                    {
                        const label i = cI - curStencil.size();
                        const label globalFaceID = cellBoundaryFaces[cellI][i];

                        forAll(boundaryMesh, patchI)
                        {
                            if
                            (
                                includePatchInStencils_[patchI]
                                && boundaryMesh[patchI].type() != emptyPolyPatch::typeName
                                && !boundaryMesh[patchI].coupled()
                            )
                            {
                                const label start = boundaryMesh[patchI].start();
                                const label nFaces = boundaryMesh[patchI].nFaces();

                                if (globalFaceID >= start && globalFaceID < start + nFaces)
                                {
                                    const label k = globalFaceID - start;
                                    y[cI] = D.boundaryField()[patchI][k][cmptI];
                                }
                            }
                        }
                    }
                }

                // Compute y_hat = W^{1/2} * y
                const Eigen::VectorXd yhat =
                    sqrtWPtr_()[cellI].diagonal().array()*y.array();

                // Compute Q^T W y = Q_hat^T * y_hat
                const Eigen::VectorXd QTWy = QhatPtr_()[cellI].transpose()*yhat;

                // Solve for A
                const Eigen::VectorXd z = choleskyPtr_()[cellI].matrixL().solve(QTWy);
                const Eigen::VectorXd A = choleskyPtr_()[cellI].matrixU().solve(z);

                // Extract gradient components from A and assign the gradient field
                // Careful: they are column-wise
                // gradD[cellI][3*cmptI] = A[1];
                // gradD[cellI][3*cmptI + 1] = A[2];
                // gradD[cellI][3*cmptI + 2] = A[3];
                gradD[cellI][3*0 + cmptI] = A[1];
                gradD[cellI][3*1 + cmptI] = A[2];
                gradD[cellI][3*2 + cmptI] = A[3];
            }
        }
    }

    gradD.correctBoundaryConditions();

    // Info<< "higherOrderGrad::grad(...): end" << endl;

    return tgradD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
