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
#include "fixedDisplacementFvPatchVectorField.H"

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


scalar kExactLeastSquares::secondDerivativeMomentContraction
(
    const symmTensor& secondDerivative,
    const symmTensor& secondMoment
)
{
    return
        secondDerivative.xx()*secondMoment.xx()
      + 2.0*secondDerivative.xy()*secondMoment.xy()
      + 2.0*secondDerivative.xz()*secondMoment.xz()
      + secondDerivative.yy()*secondMoment.yy()
      + 2.0*secondDerivative.yz()*secondMoment.yz()
      + secondDerivative.zz()*secondMoment.zz();
}


scalar kExactLeastSquares::thirdDerivativeMomentContraction
(
    const symmTensor3rdOrder& thirdDerivative,
    const symmTensor3rdOrder& thirdMoment
)
{
    return
        thirdDerivative.xxx()*thirdMoment.xxx()
      + 3.0*thirdDerivative.xxy()*thirdMoment.xxy()
      + 3.0*thirdDerivative.xxz()*thirdMoment.xxz()
      + 3.0*thirdDerivative.xyy()*thirdMoment.xyy()
      + 6.0*thirdDerivative.xyz()*thirdMoment.xyz()
      + 3.0*thirdDerivative.xzz()*thirdMoment.xzz()
      + thirdDerivative.yyy()*thirdMoment.yyy()
      + 3.0*thirdDerivative.yyz()*thirdMoment.yyz()
      + 3.0*thirdDerivative.yzz()*thirdMoment.yzz()
      + thirdDerivative.zzz()*thirdMoment.zzz();
}


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


const List<List<FixedList<label, 2>>>&
kExactLeastSquares::faceBoundaryDataStencil() const
{
    if (faceBoundaryDataStencilPtr_.empty())
    {
        calcFaceCoeffs();
    }

    return autoPtrRef(faceBoundaryDataStencilPtr_);
}


const List<CompactListList<vector>>&
kExactLeastSquares::faceBoundaryDataCoeffs() const
{
    if (faceBoundaryDataCoeffsPtr_.empty())
    {
        calcFaceCoeffs();
    }

    return autoPtrRef(faceBoundaryDataCoeffsPtr_);
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


void kExactLeastSquares::cellValueCoeffsAtPoint
(
    const label cellI,
    const point& x,
    UList<scalar>& coeffs
) const
{
    const UList<label>& cellStencil = stencil().cellsStencil()[cellI];

    if (coeffs.size() != cellStencil.size() + 1)
    {
        FatalErrorInFunction
            << "Coefficient list size " << coeffs.size()
            << " does not match stencil plus owner size "
            << cellStencil.size() + 1 << " for cell " << cellI
            << abort(FatalError);
    }

    const vector r = x - mesh_.C()[cellI];
    const UList<vector>& gradCoeffs = cellGradCoeffs()[cellI];
    const CompactListList<symmTensor>* secondDerivativeCoeffsPtr = nullptr;
    const CompactListList<symmTensor3rdOrder>*
        thirdDerivativeCoeffsPtr = nullptr;
    const List<symmTensor>* secondMomentsPtr = nullptr;
    const List<symmTensor3rdOrder>* thirdMomentsPtr = nullptr;

    if (polynomialOrder() >= 2)
    {
        secondDerivativeCoeffsPtr = &cellSecondGradCoeffs();
        secondMomentsPtr = &secondOrderCellMoments();
    }
    if (polynomialOrder() >= 3)
    {
        thirdDerivativeCoeffsPtr = &cellThirdGradCoeffs();
        thirdMomentsPtr = &thirdOrderCellMoments();
    }

    scalar stencilCoeffSum = 0.0;

    forAll(cellStencil, stencilI)
    {
        scalar& coeff = coeffs[stencilI];
        coeff = r & gradCoeffs[stencilI];

        if (secondDerivativeCoeffsPtr)
        {
            const symmTensor& secondDerivativeCoeff =
                (*secondDerivativeCoeffsPtr)[cellI][stencilI];

            coeff +=
                0.5
               *(
                    (r & (secondDerivativeCoeff & r))
                  - secondDerivativeMomentContraction
                    (
                        secondDerivativeCoeff,
                        (*secondMomentsPtr)[cellI]
                    )
                );
        }
        if (thirdDerivativeCoeffsPtr)
        {
            const symmTensor3rdOrder& thirdDerivativeCoeff =
                (*thirdDerivativeCoeffsPtr)[cellI][stencilI];

            coeff +=
                (1.0/6.0)
               *(
                    cubicForm(thirdDerivativeCoeff, r)
                  - thirdDerivativeMomentContraction
                    (
                        thirdDerivativeCoeff,
                        (*thirdMomentsPtr)[cellI]
                    )
                );
        }

        stencilCoeffSum += coeff;
    }

    coeffs[cellStencil.size()] = 1.0 - stencilCoeffSum;
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

    // A fixed-value boundary face uses the owner-cell stencil and the owner
    // cell. Prescribed boundary values have separate addressing because they
    // are known data rather than global cell unknowns.
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (!includePatchInStencils_[patchI])
        {
            continue;
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;
            const label ownCellI = owner[faceI];
            const labelUList& ownStencil = cellStencils[ownCellI];

            labelHashSet boundaryStencil(ownStencil.size() + 1);
            boundaryStencil.insert(globalCells.toGlobal(ownCellI));

            forAll(ownStencil, stencilI)
            {
                boundaryStencil.insert(ownStencil[stencilI]);
            }

            labelList boundaryIDs(boundaryStencil.toc());
            sort(boundaryIDs);
            faceStencils[faceI].transfer(boundaryIDs);
        }
    }

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

    if
    (
        faceGradCoeffsPtr_.valid()
     || faceBoundaryDataStencilPtr_.valid()
     || faceBoundaryDataCoeffsPtr_.valid()
    )
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

    faceBoundaryDataStencilPtr_.set
    (
        new List<List<FixedList<label, 2>>>(mesh.nFaces())
    );
    List<List<FixedList<label, 2>>>& faceBoundaryDataStencil =
        *faceBoundaryDataStencilPtr_;

    faceBoundaryDataCoeffsPtr_.set
    (
        new List<CompactListList<vector>>(mesh.nFaces())
    );
    List<CompactListList<vector>>& faceBoundaryDataCoeffs =
        *faceBoundaryDataCoeffsPtr_;

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

    // Gather all fixed-value quadrature points belonging to each local cell.
    // One weighted reconstruction is used for all fixed-value faces of the
    // same cell.
    List<DynamicList<FixedList<label, 2>>> cellBoundaryData(mesh.nCells());

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (!includePatchInStencils_[patchI])
        {
            continue;
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;
            const label ownCellI = owner[faceI];

            forAll(faceQuadPoints[faceI], qpI)
            {
                FixedList<label, 2> address;
                address[0] = faceI;
                address[1] = qpI;
                cellBoundaryData[ownCellI].append(address);
            }
        }
    }

    const bool twoD = mesh.nGeometricD() == 2;
    const vectorField& CI = mesh.C().internalField();
    const Map<vector>& remoteCI = stencil().remoteCentresMap();

    DynamicList<FixedList<label, 3>> exponents;
    generateExponents(polynomialOrder_, exponents);
    const label Np = exponents.size();

    List<scalar> factorials(polynomialOrder_ + 1, 1.0);
    for (label n = 1; n <= polynomialOrder_; ++n)
    {
        factorials[n] = factorials[n - 1]*n;
    }

    const CompactListList<vector>& gradCoeffs = cellGradCoeffs();
    const CompactListList<symmTensor>* secondGradCoeffsPtr = nullptr;
    const CompactListList<symmTensor3rdOrder>* thirdGradCoeffsPtr = nullptr;
    const List<symmTensor>* secondMomentsPtr = nullptr;
    const List<symmTensor3rdOrder>* thirdMomentsPtr = nullptr;

    if (polynomialOrder() >= 2)
    {
        secondGradCoeffsPtr = &cellSecondGradCoeffs();
        secondMomentsPtr = &secondOrderCellMoments();
    }
    if (polynomialOrder() >= 3)
    {
        thirdGradCoeffsPtr = &cellThirdGradCoeffs();
        thirdMomentsPtr = &thirdOrderCellMoments();
    }

    const auto cellCentre =
    [&](const label globalCellID) -> vector
    {
        if (globalCells.isLocal(globalCellID))
        {
            return CI[globalCells.toLocal(globalCellID)];
        }

        return remoteCI[globalCellID];
    };

    const auto ownerCentralMoment =
    [&]
    (
        const label cellI,
        const label i,
        const label j,
        const label k
    ) -> scalar
    {
        const label order = i + j + k;

        if (order == 1)
        {
            return 0.0;
        }
        else if (order == 2)
        {
            const symmTensor& moment = (*secondMomentsPtr)[cellI];

            if (i == 2) return moment.xx();
            if (j == 2) return moment.yy();
            if (k == 2) return moment.zz();
            if (i == 1 && j == 1) return moment.xy();
            if (i == 1 && k == 1) return moment.xz();
            if (j == 1 && k == 1) return moment.yz();
        }
        else if (order == 3)
        {
            const symmTensor3rdOrder& moment = (*thirdMomentsPtr)[cellI];

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
            << "Unsupported owner central moment exponent ("
            << i << "," << j << "," << k << ")"
            << abort(FatalError);

        return 0.0;
    };

    const auto scaledDerivativeCoeff =
    [&]
    (
        const label cellI,
        const label stencilI,
        const FixedList<label, 3>& exponent,
        const scalar h
    ) -> scalar
    {
        const label i = exponent[0];
        const label j = exponent[1];
        const label k = twoD ? 0 : exponent[2];
        const label order = i + j + k;

        if (order == 1)
        {
            const vector& coeff = gradCoeffs[cellI][stencilI];

            if (i == 1) return h*coeff.x();
            if (j == 1) return h*coeff.y();
            if (k == 1) return h*coeff.z();
        }
        else if (order == 2)
        {
            const symmTensor& coeff =
                (*secondGradCoeffsPtr)[cellI][stencilI];
            const scalar h2 = h*h;

            if (i == 2) return h2*coeff.xx();
            if (j == 2) return h2*coeff.yy();
            if (k == 2) return h2*coeff.zz();
            if (i == 1 && j == 1) return h2*coeff.xy();
            if (i == 1 && k == 1) return h2*coeff.xz();
            if (j == 1 && k == 1) return h2*coeff.yz();
        }
        else if (order == 3)
        {
            const symmTensor3rdOrder& coeff =
                (*thirdGradCoeffsPtr)[cellI][stencilI];
            const scalar h3 = h*h*h;

            if (i == 3) return h3*coeff.xxx();
            if (j == 3) return h3*coeff.yyy();
            if (k == 3) return h3*coeff.zzz();
            if (i == 2 && j == 1) return h3*coeff.xxy();
            if (i == 2 && k == 1) return h3*coeff.xxz();
            if (i == 1 && j == 2) return h3*coeff.xyy();
            if (i == 1 && j == 1 && k == 1) return h3*coeff.xyz();
            if (i == 1 && k == 2) return h3*coeff.xzz();
            if (j == 2 && k == 1) return h3*coeff.yyz();
            if (j == 1 && k == 2) return h3*coeff.yzz();
        }

        FatalErrorInFunction
            << "Unsupported scaled derivative exponent ("
            << i << "," << j << "," << k << ")"
            << abort(FatalError);

        return 0.0;
    };

    // Calculate weighted coefficients for fixed-value boundary faces.
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (!includePatchInStencils_[patchI])
        {
            continue;
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;
            const label ownCellI = owner[faceI];
            const labelUList& cellStencil = cellStencils[ownCellI];
            const labelUList& faceStencil = faceStencils[faceI];
            const DynamicList<FixedList<label, 2>>& boundaryData =
                cellBoundaryData[ownCellI];

            faceBoundaryDataStencil[faceI] =
                List<FixedList<label, 2>>(boundaryData);

            const label Nn = cellStencil.size();
            const label Nc = boundaryData.size();

            const vector& cellC = CI[ownCellI];
            scalar maxDist = 0.0;

            forAll(cellStencil, stencilI)
            {
                maxDist = max
                (
                    maxDist,
                    mag(cellCentre(cellStencil[stencilI]) - cellC)
                );
            }

            const scalar h = 2.0*maxDist;
            Eigen::MatrixXd A(Np, Nn);
            Eigen::MatrixXd Winv = Eigen::MatrixXd::Zero(Nn, Nn);

            forAll(cellStencil, stencilI)
            {
                forAll(exponents, p)
                {
                    A(p, stencilI) = scaledDerivativeCoeff
                    (
                        ownCellI,
                        stencilI,
                        exponents[p],
                        h
                    );
                }

                const scalar d =
                    mag(cellCentre(cellStencil[stencilI]) - cellC);
                Winv(stencilI, stencilI) =
                    1.0/weightFunc().weight(d, maxDist);
            }

            Eigen::MatrixXd D(Nc, Np);
            Eigen::MatrixXd boundaryWinv =
                Eigen::MatrixXd::Zero(Nc, Nc);

            for (label boundaryI = 0; boundaryI < Nc; ++boundaryI)
            {
                const FixedList<label, 2>& address =
                    boundaryData[boundaryI];
                const point& x = faceQuadPoints[address[0]][address[1]];
                const vector r = x - cellC;

                boundaryWinv(boundaryI, boundaryI) =
                    1.0/weightFunc().weight(mag(r), maxDist);

                forAll(exponents, p)
                {
                    const FixedList<label, 3>& exponent = exponents[p];
                    const label i = exponent[0];
                    const label j = exponent[1];
                    const label k = twoD ? 0 : exponent[2];
                    const label order = i + j + k;
                    const scalar factorialDenominator =
                        factorials[i]*factorials[j]*factorials[k];
                    const scalar pointMonomial =
                        pow(r.x(), i)*pow(r.y(), j)*pow(r.z(), k);

                    D(boundaryI, p) =
                        (
                            pointMonomial
                          - ownerCentralMoment(ownCellI, i, j, k)
                        )
                       /(pow(h, order)*factorialDenominator);
                }
            }

            // A maps neighbour average differences to the cell-only scaled
            // derivative vector. Its weighted covariance gives the inverse
            // cell-only least-squares Hessian. The boundary system is the
            // Woodbury form of adding the weighted boundary rows to that
            // least-squares problem.
            const Eigen::MatrixXd Hinv = A*Winv*A.transpose();
            const Eigen::MatrixXd boundarySystem =
                boundaryWinv + D*Hinv*D.transpose();
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> boundaryQr
            (
                boundarySystem
            );

            const Eigen::MatrixXd correction =
                Hinv
               *D.transpose()
               *boundaryQr.solve(Eigen::MatrixXd::Identity(Nc, Nc));
            const Eigen::MatrixXd cellMap =
                (Eigen::MatrixXd::Identity(Np, Np) - correction*D)*A;
            const Eigen::MatrixXd boundaryMap = correction;

            labelList cellRowSizes
            (
                faceQuadPoints[faceI].size(),
                faceStencil.size()
            );
            faceGradCoeffs[faceI] =
                CompactListList<vector>(cellRowSizes);

            labelList boundaryRowSizes
            (
                faceQuadPoints[faceI].size(),
                Nc
            );
            faceBoundaryDataCoeffs[faceI] =
                CompactListList<vector>(boundaryRowSizes);

            Map<label> faceCellToIndex(2*faceStencil.size());
            forAll(faceStencil, stencilI)
            {
                faceCellToIndex.insert(faceStencil[stencilI], stencilI);
            }

            const label ownIndex =
                faceCellToIndex[globalCells.toGlobal(ownCellI)];

            forAll(faceQuadPoints[faceI], qpI)
            {
                const vector r = faceQuadPoints[faceI][qpI] - cellC;
                Eigen::MatrixXd L = Eigen::MatrixXd::Zero(3, Np);

                forAll(exponents, p)
                {
                    const FixedList<label, 3>& exponent = exponents[p];
                    const label i = exponent[0];
                    const label j = exponent[1];
                    const label k = twoD ? 0 : exponent[2];
                    const label order = i + j + k;
                    const scalar denominator =
                        pow(h, order)
                       *factorials[i]*factorials[j]*factorials[k];

                    if (i > 0)
                    {
                        L(0, p) =
                            i*pow(r.x(), i - 1)
                             *pow(r.y(), j)*pow(r.z(), k)/denominator;
                    }
                    if (j > 0)
                    {
                        L(1, p) =
                            j*pow(r.x(), i)
                             *pow(r.y(), j - 1)*pow(r.z(), k)/denominator;
                    }
                    if (k > 0)
                    {
                        L(2, p) =
                            k*pow(r.x(), i)
                             *pow(r.y(), j)*pow(r.z(), k - 1)/denominator;
                    }
                }

                const Eigen::MatrixXd cellGradientMap = L*cellMap;
                const Eigen::MatrixXd boundaryGradientMap = L*boundaryMap;
                UList<vector> coefficients = faceGradCoeffs[faceI][qpI];
                UList<vector> boundaryCoefficients =
                    faceBoundaryDataCoeffs[faceI][qpI];

                forAll(coefficients, stencilI)
                {
                    coefficients[stencilI] = vector::zero;
                }

                vector coefficientSum = vector::zero;

                forAll(cellStencil, stencilI)
                {
                    const vector coefficient
                    (
                        cellGradientMap(0, stencilI),
                        cellGradientMap(1, stencilI),
                        cellGradientMap(2, stencilI)
                    );
                    const label faceStencilI =
                        faceCellToIndex[cellStencil[stencilI]];

                    coefficients[faceStencilI] += coefficient;
                    coefficientSum += coefficient;
                }

                for (label boundaryI = 0; boundaryI < Nc; ++boundaryI)
                {
                    const vector coefficient
                    (
                        boundaryGradientMap(0, boundaryI),
                        boundaryGradientMap(1, boundaryI),
                        boundaryGradientMap(2, boundaryI)
                    );

                    boundaryCoefficients[boundaryI] = coefficient;
                    coefficientSum += coefficient;
                }

                coefficients[ownIndex] -= coefficientSum;
            }
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

            // Store second-derivative tensor
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
    faceBoundaryDataStencilPtr_(),
    faceBoundaryDataCoeffsPtr_(),
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
    clearFaceCentreValueCoeffs();
    stencilPtr_.clear();
    quadraturePtr_.clear();
    weightFuncPtr_.clear();
    cellGradCoeffsPtr_.clear();
    cellSecondGradCoeffsPtr_.clear();
    cellThirdGradCoeffsPtr_.clear();
    faceGradStencilPtr_.clear();
    faceGradCoeffsPtr_.clear();
    faceBoundaryDataStencilPtr_.clear();
    faceBoundaryDataCoeffsPtr_.clear();
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
    const volVectorField& vf,
    const label patchI
) const
{
    const fvPatchField<vector>& pf = vf.boundaryField()[patchI];

    if (isA<fixedDisplacementFvPatchVectorField>(pf))
    {
        const fixedDisplacementFvPatchVectorField& patchField =
            refCast<const fixedDisplacementFvPatchVectorField>(pf);

        return patchField.evaluateQuadrature();
    }

    FatalErrorInFunction
        << "Patch " << patchI << " fixes displacement, but quadrature-point "
        << "evaluation is only implemented for "
        << fixedDisplacementFvPatchVectorField::typeName
        << abort(FatalError);

    return autoPtr<CompactListList<vector>>();
}


scalar kExactLeastSquares::valueAtPoint
(
    const volScalarField& vf,
    const label cellID,
    const point& x
) const
{
    return this->evaluateAtPoint<scalar>(vf, cellID, x);
}


vector kExactLeastSquares::valueAtPoint
(
    const volVectorField& vf,
    const label cellID,
    const point& x
) const
{
    return this->evaluateAtPoint<vector>(vf, cellID, x);
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
