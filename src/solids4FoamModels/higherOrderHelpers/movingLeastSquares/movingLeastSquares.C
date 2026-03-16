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

#include "movingLeastSquares.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"

#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedGradientFvPatchFields.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(movingLeastSquares, 0);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

const CompactListList<scalar>& movingLeastSquares::cellInterpCoeffs() const
{
    if (!cellInterpCoeffsPtr_)
    {
        calcCellCoeffs();
    }

    return cellInterpCoeffsPtr_();
}


const CompactListList<vector>& movingLeastSquares::cellGradCoeffs() const
{
    if (!cellGradCoeffsPtr_)
    {
        calcCellCoeffs();
    }

    return cellGradCoeffsPtr_();
}


const CompactListList<symmTensor>&
movingLeastSquares::cellSecondGradCoeffs() const
{
    if (!cellSecondGradCoeffsPtr_)
    {
        calcCellCoeffs();
    }

    return cellSecondGradCoeffsPtr_();
}


const CompactListList<symmTensor3rdOrder>&
movingLeastSquares::cellThirdGradCoeffs() const
{
    if (!cellThirdGradCoeffsPtr_)
    {
        calcCellCoeffs();
    }

    return cellThirdGradCoeffsPtr_();
}


const List<CompactListList<vector>>&
movingLeastSquares::faceGradCoeffs() const
{
    if (!faceGradCoeffsPtr_)
    {
        calcFaceCoeffs();
    }

    return faceGradCoeffsPtr_();
}


volScalarField& movingLeastSquares::cellConditionNumber() const
{
    if (!cellConditionNumberPtr_)
    {
        makeCellConditionNumber();
    }

    return cellConditionNumberPtr_();
}


void movingLeastSquares::makeCellConditionNumber() const
{
    if (cellConditionNumberPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    cellConditionNumberPtr_.set
    (
        new volScalarField
        (
           IOobject
           (
               "cellConditionNumber",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh_,
           dimensionedScalar("0", dimless, Zero),
           "zeroGradient"
        )
    );
}


surfaceScalarField& movingLeastSquares::faceConditionNumber() const
{
    if (!faceConditionNumberPtr_)
    {
        makeFaceConditionNumber();
    }

    return faceConditionNumberPtr_();
}


void movingLeastSquares::makeFaceConditionNumber() const
{
    if (faceConditionNumberPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    faceConditionNumberPtr_.set
    (
        new surfaceScalarField
        (
           IOobject
           (
               "faceConditionNumber",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh_,
           dimensionedScalar("0", dimless, Zero)
        )
    );
}


void movingLeastSquares::generateExponents
(
    const label N,
    DynamicList<FixedList<label, 3>>& exponents
) const
{
    if (N < 1)
    {
        FatalErrorInFunction
            << "N must be at least 1!" << abort(FatalError);
    }

    // Get the number of Taylor terms to set the capacity
    const label estimatedSize = minNn();
    exponents.setCapacity(estimatedSize);

    // Add the constant term first
    exponents.append(FixedList<label, 3>{0, 0, 0});

    // 2D and 3D cases have different number of exponents in Taylor series
    if (mesh_.nGeometricD() == 2)
    {
        for (label n = 1; n <= N; ++n)
        {
            for (label i = n; i >= 0; --i)
            {
                const label j = n - i;
                FixedList<label, 3> exponent  = {i, j, 1};
                exponents.append(exponent);
            }
        }
    }
    else
    {
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
    }

    // Adjust capacity to actual size
    exponents.shrink();
}


label movingLeastSquares::minNn() const
{
    // Taylor order in 2D case does not have terms related to z coordinate
    if (mesh_.nGeometricD() == 2)
    {
        return ((N_+1)*(N_+2)/2);
    }
    else
    {
        return ((N_+1)*(N_+2)*(N_+3)/6);
    }
}


label movingLeastSquares::rowOf
(
    const List<FixedList<label, 3>>& exponents,
    label i,
    label j,
    label k
) const
{
    const bool twoD = mesh_.nGeometricD() == 2;

    for (label r = 0; r < exponents.size(); ++r)
    {
        const FixedList<label, 3>& e = exponents[r];
        if (e[0] == i && e[1] == j && (twoD || e[2] == k))
        {
            return r;
        }
    }

    FatalErrorInFunction
        << "Missing exponent (" << i << "," << j << "," << k << ") in basis!"
        << abort(FatalError);

    //  Keeps compiler happy
    return -1;
}


movingLeastSquares::derivativeRows
movingLeastSquares::calcDerivativeRows
(
    const DynamicList<FixedList<label, 3>>& exponents,
    const bool twoD
) const
{
    derivativeRows rows;

    if (order() >= 2)
    {
        rows.ixx = rowOf(exponents, 2, 0, 0);
        rows.iyy = rowOf(exponents, 0, 2, 0);
        rows.ixy = rowOf(exponents, 1, 1, 0);

        if (!twoD)
        {
            rows.izz = rowOf(exponents, 0, 0, 2);
            rows.ixz = rowOf(exponents, 1, 0, 1);
            rows.iyz = rowOf(exponents, 0, 1, 1);
        }
    }

    if (order() >= 3)
    {
        rows.ixxx = rowOf(exponents, 3, 0, 0);
        rows.ixxy = rowOf(exponents, 2, 1, 0);
        rows.ixxz = twoD ? -1 : rowOf(exponents, 2, 0, 1);
        rows.ixyy = rowOf(exponents, 1, 2, 0);
        rows.ixyz = twoD ? -1 : rowOf(exponents, 1, 1, 1);
        rows.ixzz = twoD ? -1 : rowOf(exponents, 1, 0, 2);

        rows.iyyy = rowOf(exponents, 0, 3, 0);
        rows.iyyz = twoD ? -1 : rowOf(exponents, 0, 2, 1);
        rows.iyzz = twoD ? -1 : rowOf(exponents, 0, 1, 2);
        rows.izzz = twoD ? -1 : rowOf(exponents, 0, 0, 3);
    }

    return rows;
}


movingLeastSquares::QRSolution movingLeastSquares::QRSolve
(
    const Eigen::MatrixXd& Q,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& W
) const
{
     QRSolution result;

     const label Np = Q.rows();

     const Eigen::DiagonalMatrix<double, Eigen::Dynamic> sqrtW =
         W.diagonal().cwiseSqrt().asDiagonal();

     const Eigen::MatrixXd Qhat =
         Q.array().rowwise()*sqrtW.diagonal().transpose().array();

     const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& Bhat =
         sqrtW.diagonal().asDiagonal();

     // QR decomposition
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
           * Bhat.diagonal().transpose().array()
         ).matrix();

     result.A = Rbar.colPivHouseholderQr().solve(QbarBhat);

     if(calcConditionNumber_)
     {
         Eigen::JacobiSVD<Eigen::MatrixXd> svd
         (
             Rbar,
             Eigen::ComputeFullU | Eigen::ComputeFullV
         );

         const Eigen::VectorXd singularValues = svd.singularValues();

         result.cond =
             singularValues(0)
             /(singularValues(singularValues.size() - 1) + VSMALL);
     }

     return result;
}


void movingLeastSquares::calcCellCoeffs() const
{

    if
    (
        cellInterpCoeffsPtr_.valid()
     || cellGradCoeffsPtr_.valid()
     || cellSecondGradCoeffsPtr_.valid()
     || cellThirdGradCoeffsPtr_.valid()
    )
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    // References for brevity
    const fvMesh& mesh = mesh_;
    const bool twoD = mesh.nGeometricD() == 2;
    const globalIndex& globalCells = stencil_.globalCells();
    const vectorField& CI = mesh.C().internalField();
    const CompactListList<label>& stencils = stencil_.cellsStencil();

    // Allocate CompactListList, size is stencil plus cell centre
    labelList rowSizes(mesh.nCells());
    forAll(stencils, cI)
    {
        rowSizes[cI] = stencils[cI].size() + 1;
    }

    // Allocate CompactListList for interpolation coefficients
    cellInterpCoeffsPtr_.set(new CompactListList<scalar>(rowSizes));
    CompactListList<scalar>& cellInterpCoeffs = *cellInterpCoeffsPtr_;

    // Allocate CompactListList for gradient interpolation coefficients
    cellGradCoeffsPtr_.set(new CompactListList<vector>(rowSizes));
    CompactListList<vector>& cellGradCoeffs = *cellGradCoeffsPtr_;

    if (order() >= 2)
    {
        cellSecondGradCoeffsPtr_.set
        (
            new CompactListList<symmTensor>(rowSizes)
        );
    }
    if (order() >= 3)
    {
        cellThirdGradCoeffsPtr_.set
        (
            new CompactListList<symmTensor3rdOrder>(rowSizes)
        );
    }

    // Calculate Taylor series exponents, exponents differs for 2D and 3D case
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();

    // Compute and store derivatives rows position in matrix A
    const derivativeRows dRows = calcDerivativeRows(exponents, twoD);

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    // Reference to remote centres
    const Map<vector>& remoteCI = stencil_.remoteCentresMap();

    // Definition of Lambda function which will be used to get cell centres
    const auto cellCentre =
    [&](const label globalCellID) -> vector
    {
        if (globalCells.isLocal(globalCellID))
        {
            return CI[globalCells.toLocal(globalCellID)];
        }

        return remoteCI[globalCellID];
    };

    forAll(stencils, cellI)
    {
        const labelList& stencil = stencils[cellI];
        const vector& cellC = CI[cellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(stencil, cI)
        {
            const label neiGlobalCellID = stencil[cI];
            scalar d = mag(cellCentre(neiGlobalCellID) - cellC);

            maxDist = max(maxDist, d);
        }

        // Scaling factor for Taylor series
        const scalar h = 2.0 * maxDist;

        // In number of neigbours we will add cell centre itself
        const label Nn = stencil.size() + 1;

        // Initialise Q size, every entry in Q is set below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in aylor order!" << nl
                << "Stencil size = " << Nn << ", Taylor elements = " << Np
                << abort(FatalError);
        }

        // Loop over stencil points
        for (label cI = 0; cI < stencil.size(); ++cI)
        {
            const label neiGlobalCellID = stencil[cI];
            vector dx = cellCentre(neiGlobalCellID) - cellC;

            // Normalise dx to improve conditioning
            dx /= h;

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

                // Compute monomial values for each exponent
                if (twoD)
                {
                    Q(p, cI) =
                        pow(dx.x(), i)*pow(dx.y(), j)/factorialDenominator;
                }
                else
                {
                    Q(p, cI) =
                        pow(dx.x(), i)*pow(dx.y(), j)*pow(dx.z(), k)
                        /factorialDenominator;
                }
            }
        }

        // Add the final column for the cell-center point explicitly
        // monomials at dx=0 → [1, 0, 0, 0, ...]
        Q.col(Nn-1).setZero();
        Q(0, Nn-1) = 1.0;

        // Build W matrix
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);

        for (label cI = 0; cI < stencil.size(); cI++)
        {
            const label neiGlobalCellID = stencil[cI];
            scalar d = mag(cellCentre(neiGlobalCellID) - cellC);

            W.diagonal()[cI] = weightFunc_.weight(d, maxDist);
        }

        // Cell-centre point weight
        W.diagonal()(Nn-1) = 1.0;

        QRSolution qrs = QRSolve(Q, W);

        const Eigen::MatrixXd& A = qrs.A;

        if (calcConditionNumber_)
        {
            cellConditionNumber()[cellI] = qrs.cond;
        }

        Eigen::RowVectorXd cRow = A.row(0);
        Eigen::RowVectorXd cxRow = A.row(1)/h;
        Eigen::RowVectorXd cyRow = A.row(2)/h;
        Eigen::RowVectorXd czRow =
            twoD ? Eigen::RowVectorXd::Zero(A.cols()) : (A.row(3)/h).eval();

        for (label i = 0; i < A.cols(); ++i)
        {
            cellInterpCoeffs[cellI][i] = scalar(cRow(i));
            cellGradCoeffs[cellI][i] = vector(cxRow(i), cyRow(i), czRow(i));
        }

        if (order() >= 2)
        {
            const scalar invh2 = 1.0/(h*h);

            Eigen::RowVectorXd cxxRow = A.row(dRows.ixx) * invh2;
            Eigen::RowVectorXd cyyRow = A.row(dRows.iyy) * invh2;
            Eigen::RowVectorXd cxyRow = A.row(dRows.ixy) * invh2;

            Eigen::RowVectorXd czzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.izz) * invh2).eval();

            Eigen::RowVectorXd cxzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.ixz) * invh2).eval();

            Eigen::RowVectorXd cyzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols())
                : (A.row(dRows.iyz) * invh2).eval();

            // Store Hessian tensor
            for (label i = 0; i < A.cols(); ++i)
            {
                (*cellSecondGradCoeffsPtr_)[cellI][i] =
                    symmTensor
                    (
                        cxxRow(i),
                        cxyRow(i),
                        cxzRow(i),
                        cyyRow(i),
                        cyzRow(i),
                        czzRow(i)
                    );
            }
        }
        if (order() >= 3)
        {
            const scalar invh3 = 1.0/(h*h*h);

            Eigen::RowVectorXd cxxxRow = A.row(dRows.ixxx) * invh3;
            Eigen::RowVectorXd cxxyRow = A.row(dRows.ixxy) * invh3;
            Eigen::RowVectorXd cyyyRow = A.row(dRows.iyyy) * invh3;
            Eigen::RowVectorXd cxyyRow = A.row(dRows.ixyy) * invh3;

            Eigen::RowVectorXd cxyzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.ixyz) * invh3).eval();
            Eigen::RowVectorXd cxzzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.ixzz) * invh3).eval();
            Eigen::RowVectorXd cxxzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.ixxz) * invh3).eval();
            Eigen::RowVectorXd cyyzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.iyyz) * invh3).eval();
            Eigen::RowVectorXd cyzzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.iyzz) * invh3).eval();
            Eigen::RowVectorXd czzzRow =
                twoD ?
                Eigen::RowVectorXd::Zero(A.cols()) :
                (A.row(dRows.izzz) * invh3).eval();

            for (label i = 0; i < A.cols(); ++i)
            {
                symmTensor3rdOrder t
                (
                    cxxxRow(i),
                    cxxyRow(i),
                    cxxzRow(i),
                    cxyyRow(i),
                    cxyzRow(i),
                    cxzzRow(i),
                    cyyyRow(i),
                    cyyzRow(i),
                    cyzzRow(i),
                    czzzRow(i)
                );
                (*cellThirdGradCoeffsPtr_)[cellI][i] = t;
            }
        }
    }

    if (calcConditionNumber_)
    {
        Info<< "Writing " << cellConditionNumber().name() << " to time = "
            << mesh.time().value() << endl;
        cellConditionNumber().write();
    }
}


void movingLeastSquares::calcFaceCoeffs() const
{
    if (faceGradCoeffsPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    // References for brevity
    const fvMesh& mesh = mesh_;
    const bool twoD = mesh.nGeometricD() == 2;
    const vectorField& CI = mesh.C().internalField();
    const vectorField& Cf = mesh.faceCentres();
    const surfaceVectorField& Sf = mesh.Sf();
    const globalIndex& globalCells = stencil_.globalCells();

    faceGradCoeffsPtr_.set
    (
        new List<CompactListList<vector>>(mesh.nFaces())
    );
    List<CompactListList<vector>>& faceGradCoeffs = *faceGradCoeffsPtr_;

    // Calculate Taylor series exponents, exponents differs for 2D and 3D case
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(N_, exponents);
    const label Np = exponents.size();

    // Precompute factorials up to N
    List<scalar> factorials(N_ + 1, 1.0);
    for (label n = 1; n <= N_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    // Reference to face stencils, quadrature points and remote centres
    const CompactListList<label>& stencils = stencil_.facesStencil();
    const CompactListList<point>& faceQuadPts = quadrature_.faceQuadPoints();
    const Map<vector>& remoteCI = stencil_.remoteCentresMap();

    // Definition of Lambda function which will be used to get cell centres
    const auto cellCentre =
    [&](const label globalCellID) -> vector
    {
        if (globalCells.isLocal(globalCellID))
        {
            return CI[globalCells.toLocal(globalCellID)];
        }

        return remoteCI[globalCellID];
    };

    // Loop over all faces
    forAll(stencils, faceI)
    {
        const labelList& stencil = stencils[faceI];

        // Centre of current face
        const vector& faceCentre = Cf[faceI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(stencil, cI)
        {
            const label neiGlobalCellID = stencil[cI];

            scalar d = mag(cellCentre(neiGlobalCellID) - faceCentre);

            maxDist = max(maxDist, d);
        }

        // Scaling factor for Taylor series
        const scalar h = 2.0 * maxDist;

        // We need to extend stencil for ghost point at boundary. In the case
        // of the symmetru plane we need to reflect complete stencil
        bool ghostPoint = false;
        bool symmetryFace = false;

        // Face normal initialised and calculated below, used only for symmetry
        // planes
        vector faceNormal = vector::zero;

        if (!mesh.isInternalFace(faceI))
        {
            const label patchID = mesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh.boundaryMesh()[patchID];

            ghostPoint = includePatchInStencils_[patchID];

            if
            (
                isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchID])
             || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchID])
            )
            {
                symmetryFace = true;
                if (ghostPoint)
                {
                    FatalErrorInFunction
                        << "Face " << faceI << " is on symmetry plane but it is"
                        << " set to fix value" << abort(FatalError);
                }

                const label localFaceI = faceI - pp.start();
                faceNormal = Sf.boundaryField()[patchID][localFaceI];
                faceNormal /= (mag(faceNormal) + VSMALL);
            }
        }

        const label stencilSize = stencil.size();

        // Number of neighbours in stencil
        const label Nn =
            stencilSize
         + (ghostPoint ? 1 : 0)
         + (symmetryFace ? stencilSize : 0);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!" << nl
                << "Stencil size = " << Nn << ", Taylor elements = " << Np << nl
                << "Face centre = " << faceCentre << ", face = " << faceI
                << abort(FatalError);
        }

        // Face quadrature points
        const List<point> curFaceQuadPts = faceQuadPts[faceI];
        const label nbOfQuadPts = curFaceQuadPts.size();

        // Allocate CompactListList for this face
        labelList rowSizes(nbOfQuadPts, Nn);
        faceGradCoeffs[faceI] = CompactListList<vector>(rowSizes);

        // Average face condition number
        scalar avgCond = 0.0;

        // Loop over face quadrature points
        forAll(curFaceQuadPts, qpI)
        {
            const point& quadPoint = curFaceQuadPts[qpI];

            // Initialise Q size, every entry in Q is set below
            Eigen::MatrixXd Q(Np, Nn);

            // Loop over stencil points and compute Q
            for (label cI = 0; cI < Nn; ++cI)
            {
                label id = cI;

                // Stencil mirroring for symmetry plane face
                if (symmetryFace && cI >= stencilSize)
                {
                    id = cI - stencilSize;
                }

                // Ghost point is appended as last entry
                if (ghostPoint && cI == Nn - 1)
                {
                    Q(0, cI) = 1.0;
                    for (label p = 1; p < Np; ++p)
                    {
                        Q(p, cI) = 0.0;
                    }
                    continue;
                }

                const label neiGlobalCellID = stencil[id];

                vector dx = cellCentre(neiGlobalCellID) - quadPoint;

                // Mirror dx for symmetry plane ghost stencil part
                if (symmetryFace && cI >= stencilSize)
                {
                    dx = transform(I - 2.0*sqr(faceNormal), dx);
                }

                // Normalise dx to improve conditioning
                dx /= h;

                // Compute monomial values for each exponent
                for (label p = 0; p < Np; ++p)
                {
                    const FixedList<label, 3>& exponent = exponents[p];
                    const label i = exponent[0];
                    const label j = exponent[1];
                    const label k = exponent[2];

                    const scalar factorialDenominator =
                        factorials[i]*factorials[j]*factorials[k];

                    if (twoD)
                    {
                        Q(p, cI) =
                            pow(dx.x(), i)*pow(dx.y(), j)/factorialDenominator;
                    }
                    else
                    {
                        Q(p, cI) =
                            pow(dx.x(), i)
                          * pow(dx.y(), j)
                          * pow(dx.z(), k)
                          / factorialDenominator;
                    }
                }
            }
            // Build W matrix
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);

            for (label cI = 0; cI < Nn; ++cI)
            {
                // Mirrored half reuses weights
                if (symmetryFace && cI == stencilSize)
                {
                    W.diagonal().bottomRows(stencilSize)
                      = W.diagonal().topRows(stencilSize);
                    break;
                }

                // Ghost point gets unit weight
                if (ghostPoint && cI == Nn - 1)
                {
                    W.diagonal()[cI] = 1.0;
                    continue;
                }

                const label neiGlobalCellID = stencil[cI];
                scalar d = mag(cellCentre(neiGlobalCellID) - quadPoint);

                W.diagonal()[cI] = weightFunc_.weight(d, maxDist);
            }

            QRSolution qrs = QRSolve(Q, W);

            const Eigen::MatrixXd& A = qrs.A;

            if (calcConditionNumber_)
            {
                avgCond += qrs.cond*(1/scalar(nbOfQuadPts));
            }

            // Extract gradient rows and store them
            Eigen::RowVectorXd cxRow = A.row(1)/h;
            Eigen::RowVectorXd cyRow = A.row(2)/h;
            Eigen::RowVectorXd czRow =
                twoD ? Eigen::RowVectorXd::Zero(A.cols()) : (A.row(3)/h).eval();

            for (label i = 0; i < A.cols(); ++i)
            {
                 faceGradCoeffs[faceI][qpI][i] =
                     vector(cxRow(i), cyRow(i), czRow(i));
            }
        }

        if (calcConditionNumber_)
        {
            // Store face average condition number
            if (faceI < mesh_.nInternalFaces())
            {
                faceConditionNumber().internalFieldRef()[faceI] = avgCond;
            }
            else
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(faceI);
                const polyPatch& pp = mesh_.boundaryMesh()[patchID];
                const label localFaceI = faceI - pp.start();

                faceConditionNumber().boundaryFieldRef()[patchID][localFaceI]
                    = avgCond;
            }
        }
    }

    if (calcConditionNumber_)
    {
        Info<< "Writing " << faceConditionNumber().name() << " to time = "
            << mesh.time().value() << endl;
        faceConditionNumber().write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingLeastSquares::movingLeastSquares
(
    const fvMesh& mesh,
    const movingLeastSquaresStencil& stencil,
    const fvMeshQuadrature& quadrature,
    const weightFunction& weightFunc,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
:
    mesh_(mesh),
    stencil_(stencil),
    quadrature_(quadrature),
    weightFunc_(weightFunc),
    N_(readInt(dict.lookup("order"))),
    faceNn_(minNn() + readInt(dict.lookup("faceStencilExtraCells"))),
    cellNn_(minNn() + readInt(dict.lookup("cellStencilExtraCells"))),
    includePatchInStencils_(includePatchInStencils),
    calcConditionNumber_
    (
        dict.lookupOrDefault<Switch>("calcConditionNumber", false)
    ),
    cellInterpCoeffsPtr_(),
    cellGradCoeffsPtr_(),
    cellSecondGradCoeffsPtr_(),
    cellThirdGradCoeffsPtr_(),
    faceGradCoeffsPtr_(),
    cellConditionNumberPtr_(),
    faceConditionNumberPtr_()
{
    if (mesh_.nGeometricD() == 2)
    {
        const Vector<label>& solD = mesh_.solutionD();

        if (solD[vector::Z] != -1)
        {
            FatalErrorInFunction
                << "Empty direction should be vector::Z. "
                << "Please change mesh orientation!"
                << abort(FatalError);
        }
    }

    const scalar faceNnRatio = scalar(faceNn_) / scalar(minNn());
    const scalar cellNnRatio = scalar(cellNn_) / scalar(minNn());

    if (faceNnRatio < 1.0 || faceNnRatio > 3.0)
    {
        WarningInFunction
            << "Face stencil size is outside recommended range." << nl
            << "faceNn/minNn = " << faceNnRatio << nl
            << "Recommended range is approximately 1.5–2.5." << nl
            << "Check 'faceStencilExtraCells' entry in the dictionary."
            << endl;
    }

    if (cellNnRatio < 1.0 || cellNnRatio > 3.0)
    {
        WarningInFunction
        << "Cell stencil size is outside recommended range." << nl
        << "cellNn/minNn = " << cellNnRatio << nl
        << "Recommended range is approximately 1.5–2.5." << nl
        << "Check 'cellStencilExtraCells' entry in the dictionary."
        << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

movingLeastSquares::~movingLeastSquares()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingLeastSquares::clear() const
{
    cellInterpCoeffsPtr_.clear();
    cellGradCoeffsPtr_.clear();
    cellSecondGradCoeffsPtr_.clear();
    cellThirdGradCoeffsPtr_.clear();
    faceGradCoeffsPtr_.clear();
    cellConditionNumberPtr_.clear();
    faceConditionNumberPtr_.clear();
}


template<>
autoPtr<CompactListList<vector>>
movingLeastSquares::patchFaceQuadValues
(
    const GeometricField<vector, fvPatchField, volMesh>& vf,
    const label patchI
) const
{
    const fvPatchField<vector>& pf = vf.boundaryField()[patchI];

    // The problem is that evaluateQuadrature is not pure virtual function,
    // so we need to handle in this way
    if (!isA<fixedDisplacementFvPatchVectorField>(pf))
    {
        FatalErrorInFunction
            << "Patch " << patchI << " fixes value, but quadrature-point "
            << "evaluation is only implemented for "
            << fixedDisplacementFvPatchVectorField::typeName
            << abort(FatalError);
    }

    const fixedDisplacementFvPatchVectorField& patchField =
        refCast<const fixedDisplacementFvPatchVectorField>(pf);

    return patchField.evaluateQuadrature(quadrature_.faceQuadPoints());
}


tmp<volSymmTensorField> movingLeastSquares::secondGrad
(
    const volScalarField& s
) const
{
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil_.globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil_.remoteCellLocation();

    const CompactListList<label>& stencils = stencil_.cellsStencil();
    const CompactListList<symmTensor>& secondGradCoeffs =
        this->cellSecondGradCoeffs();

    const scalarField& sI = s.internalField();
    const List<Field<scalar>> remoteField = stencil_.remoteFieldPerProc(sI);

    tmp<volSymmTensorField> tSecondGrad
    (
        new volSymmTensorField
        (
            IOobject
            (
                "hessian(" + s.name() + ")",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedSymmTensor
            (
                "zero",
                s.dimensions()/sqr(dimLength),
                symmTensor::zero
            ),
            "zeroGradient"
        )
    );

    volSymmTensorField& secondGrad = tSecondGrad.ref();

    forAll(stencils, cellI)
    {
        const auto stencil = stencils[cellI];
        const auto coeffs = secondGradCoeffs[cellI];

        // Neighbour contribution
        forAll(stencil, cI)
        {
            secondGrad[cellI] +=
                coeffs[cI]
              * fieldValue
                (
                    stencil[cI],
                    globalCells,
                    remoteLoc,
                    sI,
                    remoteField
                );
        }

        // Cell-centre contribution
        secondGrad[cellI] += coeffs[stencil.size()] * sI[cellI];
    }

    secondGrad.correctBoundaryConditions();

    return tSecondGrad;
}


autoPtr<List<symmTensor3rdOrder>> movingLeastSquares::thirdGrad
(
    const volScalarField& s
) const
{
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil_.globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil_.remoteCellLocation();

    const CompactListList<label>& stencils = stencil_.cellsStencil();
    const CompactListList<symmTensor3rdOrder>& thirdGradCoeffs =
        cellThirdGradCoeffs();

    const scalarField& sI = s.internalField();
    const List<Field<scalar>> remoteField = stencil_.remoteFieldPerProc(sI);

    autoPtr<List<symmTensor3rdOrder>> tThirdGrad
    (
        new List<symmTensor3rdOrder>(mesh.nCells(), symmTensor3rdOrder::zero)
    );
    List<symmTensor3rdOrder>& thirdGrad = tThirdGrad.ref();

    forAll(stencils, cellI)
    {
        const auto stencil = stencils[cellI];
        const auto coeffs = thirdGradCoeffs[cellI];

        // Neighbour contribution
        forAll(stencil, cI)
        {
            thirdGrad[cellI] +=
                coeffs[cI]
              * fieldValue
                (
                    stencil[cI],
                    globalCells,
                    remoteLoc,
                    sI,
                    remoteField
                );
        }

        // Cell-centre contribution
        thirdGrad[cellI] += coeffs[stencil.size()] * sI[cellI];
    }

    return tThirdGrad;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
