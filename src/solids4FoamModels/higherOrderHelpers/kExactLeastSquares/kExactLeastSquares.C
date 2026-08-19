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

#include "kExactLeastSquares.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kExactLeastSquares, 0);
addToRunTimeSelectionTable
(
    leastSquaresScheme,
    kExactLeastSquares,
    dictionary
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const CompactListList<vector>& kExactLeastSquares::cellGradCoeffs() const
{
    if (cellGradCoeffsPtr_.empty())
    {
        calcCellCoeffs();
    }

    return autoPtrRef(cellGradCoeffsPtr_);
}


const CompactListList<symmTensor>&
kExactLeastSquares::cellSecondGradCoeffs() const
{
    if (cellSecondGradCoeffsPtr_.empty())
    {
        calcCellCoeffs();
    }

    return autoPtrRef(cellSecondGradCoeffsPtr_);
}


const CompactListList<symmTensor3rdOrder>&
kExactLeastSquares::cellThirdGradCoeffs() const
{
    if (cellThirdGradCoeffsPtr_.empty())
    {
        calcCellCoeffs();
    }

    return autoPtrRef(cellThirdGradCoeffsPtr_);
}


void kExactLeastSquares::cellGradCoeffsAtPoint
(
    const label cellI,
    const point& x,
    UList<vector>& coeffs
) const
{
    const CompactListList<label>& stencils = stencil().cellsStencil();
    const labelUList& cellStencil = stencils[cellI];

    if (coeffs.size() != cellStencil.size())
    {
        FatalErrorInFunction
            << "Coefficient list size " << coeffs.size()
            << " does not match stencil size " << cellStencil.size()
            << " for cell " << cellI << abort(FatalError);
    }

    const vector r = x - mesh_.C()[cellI];
    const UList<vector>& gradCoeffs = cellGradCoeffs()[cellI];

    const CompactListList<symmTensor>* secondGradCoeffsPtr = nullptr;
    const CompactListList<symmTensor3rdOrder>* thirdGradCoeffsPtr = nullptr;

    if (polynomialOrder() >= 2)
    {
        secondGradCoeffsPtr = &cellSecondGradCoeffs();
    }
    if (polynomialOrder() >= 3)
    {
        thirdGradCoeffsPtr = &cellThirdGradCoeffs();
    }

    forAll(coeffs, cI)
    {
        coeffs[cI] = gradCoeffs[cI];

        if (secondGradCoeffsPtr)
        {
            coeffs[cI] += (*secondGradCoeffsPtr)[cellI][cI] & r;
        }
        if (thirdGradCoeffsPtr)
        {
            coeffs[cI] +=
                0.5*(((*thirdGradCoeffsPtr)[cellI][cI] & r) & r);
        }
    }
}


label kExactLeastSquares::minNn() const
{
    // Taylor polynomial order in 2D case does not have terms related to z
    // coordinate
    if (mesh_.nGeometricD() == 2)
    {
        return ((polynomialOrder_ + 1)*(polynomialOrder_ + 2)/2) - 1;
    }
    else
    {
        return
        (
            (polynomialOrder_ + 1)
          * (polynomialOrder_ + 2)
          * (polynomialOrder_ + 3)
          / 6
          - 1
        );
    }
}


void kExactLeastSquares::generateExponents
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

    // 2D and 3D cases have different number of exponents in Taylor series
    if (mesh_.nGeometricD() == 2)
    {
        for (label n = 1; n <= N; ++n)
        {
            for (label i = n; i >= 0; --i)
            {
                const label j = n - i;
                FixedList<label, 3> exponent = {i, j, 0};
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
                    const label k = n - i - j;
                    FixedList<label, 3> exponent = {i, j, k};
                    exponents.append(exponent);
                }
            }
        }
    }

    // Adjust capacity to actual size
    exponents.shrink();
}


label kExactLeastSquares::rowOf
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


kExactLeastSquares::derivativeRows
kExactLeastSquares::calcDerivativeRows
(
    const DynamicList<FixedList<label, 3>>& exponents,
    const bool twoD
) const
{
    derivativeRows rows;

    if (polynomialOrder() >= 2)
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

    if (polynomialOrder() >= 3)
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


volScalarField& kExactLeastSquares::cellConditionNumber() const
{
    if (cellConditionNumberPtr_.empty())
    {
        makeCellConditionNumber();
    }

    return autoPtrRef(cellConditionNumberPtr_);
}


void kExactLeastSquares::makeCellConditionNumber() const
{
    if (cellConditionNumberPtr_.valid())
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
               "kExactCellConditionNumber",
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


void kExactLeastSquares::makeFaceGradStencil() const
{
    if (faceGradStencilPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const CompactListList<label>& cellStencils =
        stencil().cellsStencil();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    List<labelList> faceStencils(mesh.nFaces());

    // Loop over internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        // Index of face owner and neighbour
        const label ownCellI = owner[faceI];
        const label neiCellI = neighbour[faceI];

        // Owner and neighbour stencils
        const labelUList& ownStencil = cellStencils[ownCellI];
        const labelUList& neiStencil = cellStencils[neiCellI];

        // We don't need to explicitly add  owner and neighbour but it does not
        // matter
        labelHashSet mergedStencil
        (
            ownStencil.size() + neiStencil.size() + 2
        );

        mergedStencil.insert(globalCells.toGlobal(ownCellI));
        forAll(ownStencil, cI)
        {
            mergedStencil.insert(ownStencil[cI]);
        }

        mergedStencil.insert(globalCells.toGlobal(neiCellI));
        forAll(neiStencil, cI)
        {
            mergedStencil.insert(neiStencil[cI]);
        }

        labelList mergedIDs(mergedStencil.toc());
        sort(mergedIDs);
        faceStencils[faceI].transfer(mergedIDs);
    }

    // TO_DO: paralle boundaries? patches?

    // Store face stencils to faceGradStencilPtr_
    labelList rowSizes(mesh.nFaces(), 0);
    forAll(faceStencils, faceI)
    {
        rowSizes[faceI] = faceStencils[faceI].size();
    }

    faceGradStencilPtr_.set(new CompactListList<label>(rowSizes));
    CompactListList<label>& faceGradStencils = *faceGradStencilPtr_;

    forAll(faceStencils, faceI)
    {
        const labelList& src = faceStencils[faceI];
        labelUList dst = faceGradStencils[faceI];

        forAll(src, cI)
        {
            dst[cI] = src[cI];
        }
    }
}


void kExactLeastSquares::calcFaceCoeffs() const
{
    if (debug)
    {
        InfoInFunction << "start" << endl;
    }

    if (faceGradCoeffsPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    // Preliminaries
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const CompactListList<label>& cellStencils =  stencil().cellsStencil();
    const CompactListList<label>& faceStencils = faceGradStencil();
    const CompactListList<point>& faceQuadPoints =
        quadrature().faceQuadPoints();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    // Initialise faceGradCoeffsPtr_
    faceGradCoeffsPtr_.set
    (
        new List<CompactListList<vector>>(mesh.nFaces())
    );
    List<CompactListList<vector>>& faceGradCoeffs = *faceGradCoeffsPtr_;

    // Loop over internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        const label ownCellI = owner[faceI];
        const label neiCellI = neighbour[faceI];

        const labelUList& ownStencil = cellStencils[ownCellI];
        const labelUList& neiStencil = cellStencils[neiCellI];

        const labelUList& faceStencil = faceStencils[faceI];
        const UList<point>& quadPoints = faceQuadPoints[faceI];

        // Number of coefficients is equal to merged stencil size
        labelList rowSizes(quadPoints.size(), faceStencil.size());
        faceGradCoeffs[faceI] = CompactListList<vector>(rowSizes);

        // Map global cell IDs onto positions in the merged face stencil
        Map<label> faceCellToIndex(2*faceStencil.size());
        forAll(faceStencil, stencilI)
        {
            faceCellToIndex.insert(faceStencil[stencilI], stencilI);
        }

        // Get stencil list position for owner and neighbour
        const label ownIndex =
            faceCellToIndex[globalCells.toGlobal(ownCellI)];
        const label neiIndex =
            faceCellToIndex[globalCells.toGlobal(neiCellI)];

        vectorField ownCoeffs(ownStencil.size());
        vectorField neiCoeffs(neiStencil.size());

        forAll(quadPoints, qpI)
        {
            const point& quadPoint = quadPoints[qpI];
            vector ownSum = vector::zero;
            vector neiSum = vector::zero;
            UList<vector> coeffs = faceGradCoeffs[faceI][qpI];

            // Initialise to zero all coefficients of this quadrature point
            forAll(coeffs, stencilI)
            {
                coeffs[stencilI] = vector::zero;
            }

            // Get owner side coefficients
            cellGradCoeffsAtPoint(ownCellI, quadPoint, ownCoeffs);
            forAll(ownStencil, cellI)
            {
                const vector& coeff = ownCoeffs[cellI];
                const label faceStencilI =
                    faceCellToIndex[ownStencil[cellI]];

                coeffs[faceStencilI] += 0.5*coeff;
                ownSum += coeff;
            }
            coeffs[ownIndex] -= 0.5*ownSum;

            // Get neighbour side coefficients
            cellGradCoeffsAtPoint(neiCellI, quadPoint, neiCoeffs);
            forAll(neiStencil, cellI)
            {
                const vector& coeff = neiCoeffs[cellI];
                const label faceStencilI =
                    faceCellToIndex[neiStencil[cellI]];

                coeffs[faceStencilI] += 0.5*coeff;
                neiSum += coeff;
            }
            coeffs[neiIndex] -= 0.5*neiSum;
        }
    }

    if (debug)
    {
        InfoInFunction << "end" << endl;
    }
}


void kExactLeastSquares::makeStencils() const
{
    if (stencilPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    const label minStencilSize = minNn();
    // leastSquaresStencil expects Nc to include the central cell, whereas its
    // cellsStencil() contains neighbours only
    const label cellNn = minStencilSize + cellStencilExtraCells_ + 1;

    // We will use same size for face stencil,
    // but face stencil  will not be evaluated as it is not needed
    const label faceNn = cellNn;

    stencilPtr_.set
    (
        new leastSquaresStencil
        (
            mesh_,
            haloDepthScale_,
            faceNn,
            cellNn
        )
    );
}


void kExactLeastSquares::makeQuadrature() const
{
    if (quadraturePtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    quadraturePtr_.set
    (
        new fvMeshQuadrature
        (
            mesh_,
            polynomialOrder_,
            polynomialOrder_ - 1,
            true
        )
    );
}


void kExactLeastSquares::makeWeightFunction() const
{
    if (weightFuncPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    weightFuncPtr_.set
    (
        new weightFunction(weightFunctionCoeffs_)
    );
}


kExactLeastSquares::QRSolution kExactLeastSquares::QRSolve
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


void kExactLeastSquares::calcCellCoeffs() const
{
    if (debug)
    {
        InfoInFunction << "start" << endl;
    }

    if
    (
        cellGradCoeffsPtr_.valid()
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
    const globalIndex& globalCells = stencil().globalCells();
    const vectorField& CI = mesh.C().internalField();
    const CompactListList<label>& stencils = stencil().cellsStencil();

    // Allocate CompactListList, size is stencil
    labelList rowSizes(mesh.nCells());
    forAll(stencils, cI)
    {
        rowSizes[cI] = stencils[cI].size();
    }

    // Allocate CompactListList for gradient interpolation coefficients
    cellGradCoeffsPtr_.set(new CompactListList<vector>(rowSizes));
    CompactListList<vector>& cellGradCoeffs = *cellGradCoeffsPtr_;

    if (polynomialOrder() >= 2)
    {
        cellSecondGradCoeffsPtr_.set
        (
            new CompactListList<symmTensor>(rowSizes)
        );
    }
    if (polynomialOrder() >= 3)
    {
        cellThirdGradCoeffsPtr_.set
        (
            new CompactListList<symmTensor3rdOrder>(rowSizes)
        );
    }

    // Calculate Taylor series exponents, exponents differs for 2D and 3D case
    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(polynomialOrder_, exponents);

    // The constant coefficient is fixed by the owner-cell average and is not
    // included in exponents
    const label Np = exponents.size();

    // Compute and store derivatives rows position in matrix A
    const derivativeRows dRows = calcDerivativeRows(exponents, twoD);

    // Precompute factorials up to N
    List<scalar> factorials(polynomialOrder_ + 1, 1.0);
    for (label n = 1; n <= polynomialOrder_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    // Reference to remote centres
    const Map<vector>& remoteCI = stencil().remoteCentresMap();

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

    // Cell moments are exchanged only for the polynomial orders that use them
    const List<symmTensor>* secondMomentsPtr = nullptr;
    const List<symmTensor3rdOrder>* thirdMomentsPtr = nullptr;
    Map<symmTensor> remoteSecondMoments;
    Map<symmTensor3rdOrder> remoteThirdMoments;

    if (polynomialOrder() >= 2)
    {
        secondMomentsPtr = &secondOrderCellMoments();
        remoteSecondMoments =
            stencil().remoteFieldMap(*secondMomentsPtr);
    }
    if (polynomialOrder() >= 3)
    {
        thirdMomentsPtr = &thirdOrderCellMoments();
        remoteThirdMoments =
            stencil().remoteFieldMap(*thirdMomentsPtr);
    }

    // Return a volume-normalised central moment for a local or remote cell
    const auto centralMoment =
    [&]
    (
        const label globalCellID,
        const label i,
        const label j,
        const label k
    ) -> scalar
    {
        const label order = i + j + k;

        if (order == 0)
        {
            return 1.0;
        }
        else if (order == 1)
        {
            // Cell centres are volume centroids
            return 0.0;
        }
        else if (order == 2)
        {
            const symmTensor* momentPtr = nullptr;
            if (globalCells.isLocal(globalCellID))
            {
                momentPtr =
                    &(*secondMomentsPtr)[globalCells.toLocal(globalCellID)];
            }
            else
            {
                momentPtr = &remoteSecondMoments[globalCellID];
            }

            const symmTensor& moment = *momentPtr;

            if (i == 2) return moment.xx();
            if (j == 2) return moment.yy();
            if (k == 2) return moment.zz();
            if (i == 1 && j == 1) return moment.xy();
            if (i == 1 && k == 1) return moment.xz();
            if (j == 1 && k == 1) return moment.yz();
        }
        else if (order == 3)
        {
            const symmTensor3rdOrder* momentPtr = nullptr;
            if (globalCells.isLocal(globalCellID))
            {
                momentPtr =
                    &(*thirdMomentsPtr)[globalCells.toLocal(globalCellID)];
            }
            else
            {
                momentPtr = &remoteThirdMoments[globalCellID];
            }

            const symmTensor3rdOrder& moment = *momentPtr;

            if (i == 3) return moment.xxx();
            if (j == 3) return moment.yyy();
            if (k == 3) return moment.zzz();
            if (i == 2 && j == 1) return moment.xxy();
            if (i == 2 && k == 1) return moment.xxz();
            if (i == 1 && j == 2) return moment.xyy();
            if (i == 1 && j == 1 && k == 1) return moment.xyz();
            if (i == 1 && k == 2) return moment.xzz();
            if (j == 2 && k == 1) return moment.yyz();
            if (j == 1 && k == 2) return moment.yzz();
        }

        FatalErrorInFunction
            << "Unsupported central moment exponent ("
            << i << "," << j << "," << k << ")"
            << abort(FatalError);

        // Keeps compiler happy
        return 0.0;
    };

    // Return the cell average of a monomial centred about the owner cell
    const auto averageMonomial =
    [&]
    (
        const label globalCellID,
        const vector& d,
        const label i,
        const label j,
        const label k
    ) -> scalar
    {
        scalar average = 0.0;

        for (label a = 0; a <= i; ++a)
        {
            const scalar binomialI =
                factorials[i]/(factorials[a]*factorials[i - a]);

            for (label b = 0; b <= j; ++b)
            {
                const scalar binomialJ =
                    factorials[j]/(factorials[b]*factorials[j - b]);

                for (label c = 0; c <= k; ++c)
                {
                    const scalar binomialK =
                        factorials[k]/(factorials[c]*factorials[k - c]);

                    average +=
                        binomialI*binomialJ*binomialK
                       *pow(d.x(), i - a)
                       *pow(d.y(), j - b)
                       *pow(d.z(), k - c)
                       *centralMoment(globalCellID, a, b, c);
                }
            }
        }

        return average;
    };

    forAll(stencils, cellI)
    {
        const labelUList stencil = stencils[cellI];
        const vector& cellC = CI[cellI];
        const label ownerGlobalCellID = globalCells.toGlobal(cellI);

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

        // Number of neighbours in stencil
        const label Nn = stencil.size();

        // Initialise Q size, every entry in Q is set below
        Eigen::MatrixXd Q(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
                << "number of elements in Taylor order!" << nl
                << "Stencil size = " << Nn << ", Taylor elements = " << Np
                << abort(FatalError);
        }

        // Loop over stencil cells and compute their cell-average equations
        for (label cI = 0; cI < stencil.size(); ++cI)
        {
            const label neiGlobalCellID = stencil[cI];
            const vector dx = cellCentre(neiGlobalCellID) - cellC;

            forAll(exponents, p)
            {
                const FixedList<label, 3>& exponent = exponents[p];
                const label i = exponent[0];
                const label j = exponent[1];
                const label k = twoD ? 0 : exponent[2];
                const label order = i + j + k;

                // Compute factorial denominator
                const scalar factorialDenominator =
                    factorials[i]*factorials[j]*factorials[k];

                const scalar neighbourAverage = averageMonomial
                (
                    neiGlobalCellID,
                    dx,
                    i,
                    j,
                    k
                );
                const scalar ownerAverage = averageMonomial
                (
                    ownerGlobalCellID,
                    vector::zero,
                    i,
                    j,
                    k
                );

                Q(p, cI) =
                    (neighbourAverage - ownerAverage)
                   /(pow(h, order)*factorialDenominator);
            }
        }

        // Build W matrix
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(Nn);

        for (label cI = 0; cI < stencil.size(); cI++)
        {
            const label neiGlobalCellID = stencil[cI];
            scalar d = mag(cellCentre(neiGlobalCellID) - cellC);

            W.diagonal()[cI] = weightFunc().weight(d, maxDist);
        }

        QRSolution qrs = QRSolve(Q, W);

        const Eigen::MatrixXd& A = qrs.A;

        if (calcConditionNumber_)
        {
            cellConditionNumber()[cellI] = qrs.cond;
        }

        const label ix = rowOf(exponents, 1, 0, 0);
        const label iy = rowOf(exponents, 0, 1, 0);
        const label iz = twoD ? -1 : rowOf(exponents, 0, 0, 1);

        Eigen::RowVectorXd cxRow = A.row(ix)/h;
        Eigen::RowVectorXd cyRow = A.row(iy)/h;
        Eigen::RowVectorXd czRow =
            twoD ? Eigen::RowVectorXd::Zero(A.cols()) : (A.row(iz)/h).eval();

        for (label i = 0; i < A.cols(); ++i)
        {
            cellGradCoeffs[cellI][i] =
                vector(cxRow(i), cyRow(i), czRow(i));
        }

        if (polynomialOrder() >= 2)
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
        if (polynomialOrder() >= 3)
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
                (*cellThirdGradCoeffsPtr_)[cellI][i] =
                    symmTensor3rdOrder
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
            }
        }
    }

    if (calcConditionNumber_)
    {
        Info<< "Writing " << cellConditionNumber().name() << " to time = "
            << mesh.time().value() << endl;
        cellConditionNumber().write();
    }
    if (debug)
    {
        InfoInFunction << "end" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kExactLeastSquares::kExactLeastSquares
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
:
    leastSquaresScheme(mesh),
    stencilPtr_(),
    quadraturePtr_(),
    weightFuncPtr_(),
    polynomialOrder_(readInt(dict.lookup("polynomialOrder"))),
    cellStencilExtraCells_
    (
        dict.found("cellStencilExtraCells")
      ? readInt(dict.lookup("cellStencilExtraCells"))
      : readInt(dict.lookup("faceStencilExtraCells"))
    ),
    haloDepthScale_
    (
        dict.lookupOrDefault<scalar>("haloDepthScale", polynomialOrder_*2.5)
    ),
    weightFunctionCoeffs_(dict.subDict("weightFunctionCoeffs")),
    includePatchInStencils_(includePatchInStencils),
    calcConditionNumber_
    (
        dict.lookupOrDefault<Switch>("calcConditionNumber", false)
    ),
    cellGradCoeffsPtr_(),
    cellSecondGradCoeffsPtr_(),
    cellThirdGradCoeffsPtr_(),
    faceGradStencilPtr_(),
    faceGradCoeffsPtr_(),
    cellConditionNumberPtr_()
{
    if (polynomialOrder_ > 3 || polynomialOrder_ < 1)
    {
        FatalErrorInFunction
            << "Chosen polynomial order " << polynomialOrder_
            << " is not supported. Supported orders are: 1, 2 and 3."
            << abort(FatalError);
    }

    const label minStencilSize = minNn();

    // Number of cells in cell stencil
    const label cellNn = minStencilSize + cellStencilExtraCells_;

    const scalar cellNnRatio = scalar(cellNn)/scalar(minStencilSize);

    if (cellNnRatio < 1.2 || cellNnRatio > 15.0)
    {
        WarningInFunction
            << "Cell stencil size is outside valid range." << nl
            << "cellNn/minNn = " << cellNnRatio << nl
            << "Check 'cellStencilExtraCells' entry in the dictionary."
            << endl;
    }

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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

kExactLeastSquares::~kExactLeastSquares()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

const leastSquaresStencil& kExactLeastSquares::stencil() const
{
    if (stencilPtr_.empty())
    {
        makeStencils();
    }

    return *stencilPtr_;
}


const CompactListList<label>& kExactLeastSquares::faceGradStencil() const
{
    if (faceGradStencilPtr_.empty())
    {
        makeFaceGradStencil();
    }

    return autoPtrRef(faceGradStencilPtr_);
}


const List<CompactListList<vector>>&
kExactLeastSquares::faceGradCoeffs() const
{
    if (faceGradCoeffsPtr_.empty())
    {
        calcFaceCoeffs();
    }

    return autoPtrRef(faceGradCoeffsPtr_);
}


const List<symmTensor>& kExactLeastSquares::secondOrderCellMoments() const
{
    return quadrature().secondOrderCellMoments();
}


const List<symmTensor3rdOrder>&
kExactLeastSquares::thirdOrderCellMoments() const
{
    return quadrature().thirdOrderCellMoments();
}


void kExactLeastSquares::clear() const
{
    stencilPtr_.clear();
    quadraturePtr_.clear();
    weightFuncPtr_.clear();
    cellGradCoeffsPtr_.clear();
    cellSecondGradCoeffsPtr_.clear();
    cellThirdGradCoeffsPtr_.clear();
    faceGradStencilPtr_.clear();
    faceGradCoeffsPtr_.clear();
    cellConditionNumberPtr_.clear();
}


void kExactLeastSquares::fGrad
(
    const volScalarField& vf,
    CompactListList<vector>& result
) const
{
    this->fGrad<scalar>(vf, result);
}


void kExactLeastSquares::fGrad
(
    const volVectorField& vf,
    CompactListList<tensor>& result
) const
{
    this->fGrad<vector>(vf, result);
}


autoPtr<CompactListList<scalar>> kExactLeastSquares::patchFaceQuadValues
(
    const volScalarField&,
    const label
) const
{
    NotImplemented;
    return autoPtr<CompactListList<scalar>>();
}


autoPtr<CompactListList<vector>> kExactLeastSquares::patchFaceQuadValues
(
    const volVectorField&,
    const label
) const
{
    NotImplemented;
    return autoPtr<CompactListList<vector>>();
}


tmp<volVectorField> kExactLeastSquares::grad
(
    const volScalarField& vf
) const
{
    return this->grad<scalar>(vf);
}


tmp<volTensorField> kExactLeastSquares::grad
(
    const volVectorField& vf
) const
{
    return this->grad<vector>(vf);
}


tmp<volSymmTensorField> kExactLeastSquares::secondGrad
(
    const volScalarField& s
) const
{
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil().remoteCellLocation();

    const CompactListList<label>& stencils = stencil().cellsStencil();
    const CompactListList<symmTensor>& secondGradCoeffs =
        this->cellSecondGradCoeffs();

    const scalarField& sI = s.internalField();
    const List<Field<scalar>> remoteField = stencil().remoteFieldPerProc(sI);

    tmp<volSymmTensorField> tSecondGrad
    (
        new volSymmTensorField
        (
            IOobject
            (
                "secondGrad(" + s.name() + ")",
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
        const UList<label>& stencil = stencils[cellI];
        const UList<symmTensor>& coeffs = secondGradCoeffs[cellI];

        // Neighbour cell-average difference contribution
        forAll(stencil, cI)
        {
            secondGrad[cellI] +=
                coeffs[cI]
               *
                (
                    fieldValue
                    (
                        stencil[cI],
                        globalCells,
                        remoteLoc,
                        sI,
                        remoteField
                    )
                  - sI[cellI]
                );
        }
    }

    secondGrad.correctBoundaryConditions();

    return tSecondGrad;
}


autoPtr<List<symmTensor3rdOrder>> kExactLeastSquares::thirdGrad
(
    const volScalarField& s
) const
{
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil().remoteCellLocation();

    const CompactListList<label>& stencils = stencil().cellsStencil();
    const CompactListList<symmTensor3rdOrder>& thirdGradCoeffs =
        cellThirdGradCoeffs();

    const scalarField& sI = s.internalField();
    const List<Field<scalar>> remoteField = stencil().remoteFieldPerProc(sI);

    autoPtr<List<symmTensor3rdOrder>> tThirdGrad
    (
        new List<symmTensor3rdOrder>(mesh.nCells(), symmTensor3rdOrder::zero)
    );
    List<symmTensor3rdOrder>& thirdGrad = autoPtrRef(tThirdGrad);

    forAll(stencils, cellI)
    {
        const UList<label>& stencil = stencils[cellI];
        const UList<symmTensor3rdOrder>& coeffs = thirdGradCoeffs[cellI];

        // Neighbour cell-average difference contribution
        forAll(stencil, cI)
        {
            thirdGrad[cellI] +=
                coeffs[cI]
               *
                (
                    fieldValue
                    (
                        stencil[cI],
                        globalCells,
                        remoteLoc,
                        sI,
                        remoteField
                    )
                  - sI[cellI]
                );
        }
    }

    return tThirdGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
