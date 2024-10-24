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

#include "coconutModel.H"
#include "addToRunTimeSelectionTable.H"
#include "QRMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coconutModel, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


tmp<scalarField> coconutModel::unroll
(
    const vectorField& x
) const
{
    // Prepare result field
    tmp<scalarField> tresult
    (
        new scalarField(x.size(), 0.0)
    );
    scalarField& result = tresult.ref();

    int i = 0;
    forAll(x, vecI)
    {
        result[i++] = x[vecI][vector::X];
        result[i++] = x[vecI][vector::Y];
        result[i++] = x[vecI][vector::Z];
    }

    return tresult;
}


tmp<vectorField> coconutModel::roll
(
    const scalarField& x
) const
{
    const int vecSize = x.size()/3;

    if (debug)
    {
        if (x.size() % 3 != 0)
        {
            FatalErrorInFunction
                << "x.size() % 3 != 0" << abort(FatalError);
        }
    }

    // Prepare result field
    tmp<vectorField> tresult
    (
        new vectorField(vecSize, vector::zero)
    );
    vectorField& result = tresult.ref();

    int i = 0;
    forAll(result, vecI)
    {
        result[vecI][vector::X] = x[i++];
        result[vecI][vector::Y] = x[i++];
        result[vecI][vector::Z] = x[i++];
    }

    return tresult;
}


void coconutModel::convertToMat
(
    scalarRectangularMatrix& vMat,
    const PtrList<vectorField>& v
) const
{
    if (v.size() == 0)
    {
        return;
    }

    // We assume all the vectorFields have the same length
    const int nVecRow = v[0].size();
    if (debug)
    {
        forAll(v, i)
        {
            if (v[i].size() != nVecRow)
            {
                FatalErrorInFunction
                    << "The vectorField have different lengths!"
                    << abort(FatalError);
            }
        }
    }

    // Place each vectorField in a column
    const int nRow = 3*nVecRow;
    const int nCol = v.size();
    vMat.resize(nRow, nCol);
    forAll(v, colI)
    {
        int rowI = 0;
        forAll(v[colI], vecRowI)
        {
            vMat(rowI++, colI) = v[colI][vecRowI][vector::X];
            vMat(rowI++, colI) = v[colI][vecRowI][vector::Y];
            vMat(rowI++, colI) = v[colI][vecRowI][vector::Z];
        }
    }
}


void coconutModel::combineAndLimit
(
    scalarRectangularMatrix& mat,
    const PtrList<vectorField>& aList,
    const PtrList<vectorField>& bList,
    const int modes
) const
{
    Info<< "aList.size = " << aList.size() << nl
        << "bList.size = " << bList.size() << nl
        << endl;

    // Info<< "aList[0].size = " << aList[0].size() << nl
    //     << "bList[0].size = " << bList[0].size() << nl
    //     << endl;

    // Convert PtrList<vectorField> to matrices
    scalarRectangularMatrix aMat;
    convertToMat(aMat, aList);
    scalarRectangularMatrix bMat;
    convertToMat(bMat, bList);

    Info<< "aMat.nRows() = " << aMat.nRows() << nl
        << "bMat.nRows() = " << bMat.nRows() << nl
        << endl;

    Info<< "aMat.nCols() = " << aMat.nCols() << nl
        << "bMat.nCols() = " << bMat.nCols() << nl
        << endl;

    // Limit vCurrS and vPrevS
    scalarRectangularMatrix a = limit(aMat, modes);
    scalarRectangularMatrix b = limit(bMat, modes);

    if (debug && a.nCols() > 0 && b.nCols() > 0)
    {
        if (a.nRows() != b.nRows())
        {
            FatalErrorInFunction
                << "a.nRows() != b.nRows()" << abort(FatalError);
        }
    }

    mat.resize
    (
        max(a.nRows(), b.nRows()), a.nCols() + b.nCols()
    );

    Info<< "mat.nRows() = " << mat.nRows() << nl
        << "mat.nCols() = " << mat.nCols() << nl
        << endl;

    for (int j = 0; j < a.nCols(); ++j)
    {
        for (int i = 0; i < a.nRows(); ++i)
        {
            mat(i, j) = a(i, j);
        }
    }
    for (int j = 0; j < b.nCols(); ++j)
    {
        for (int i = 0; i < b.nRows(); ++i)
        {
            mat(i, a.nCols() + j) = b(i, j);
        }
    }
}


void coconutModel::filter()
{
   notImplemented("To be implemented");

   // Filter vcurr, vprev
   // it may be easier to filter the matrix versions
}


scalarRectangularMatrix coconutModel::limit
(
    const scalarRectangularMatrix& matrix,
    const int modes
) const
{
    // If modes is negative, treat it as unspecified (similar to Python's None)
    if (modes < 0)
    {
        return matrix;
    }
    else if (modes == 0)
    {
        // Return a matrix with the same number of rows and 0 columns
        return scalarRectangularMatrix(matrix.nRows(), 0);
    }
    else
    {
        // Determine the number of columns to extract
        int colsToExtract = min(matrix.nCols(), modes);

        // Return the limited matrix (oldest columns)
        return matrix.subMatrix
        (
            0,
            matrix.nRows(),
            matrix.nCols() - colsToExtract,
            matrix.nCols()
        );

        // Alternatively, to return the most recent columns:
        // return matrix.subMatrix(0, matrix.rows(), 0, colsToExtract);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coconutModel::coconutModel(const int q)
:
    q_(q),
    added_(false),
    vCurr_(0),
    vPrev_(0),
    wCurr_(0),
    wPrev_(0),
    rRef_(0),
    xtRef_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coconutModel::initializeSolutionStep()
{
    rRef_ = vector::zero;
    xtRef_ = vector::zero;
    vCurr_.clear();
    wCurr_.clear();
    added_ = false;

    // Throw away the entries from the start until vPrev.size() <= max(q, 1)
    // Careful: vPrev/wPrev are stored in the opposite order in coconut
    for
    (
        int oldI = vPrev_.size() - max(q_, 1), newI = 0;
        oldI < vPrev_.size() && oldI >= 0;
        oldI++, newI++
    )
    {
        vPrev_[newI] = vPrev_[oldI];
    }

    if (vPrev_.size() > max(q_, 1))
    {
        vPrev_.resize(max(q_, 1));
    }
}


void coconutModel::finaliseSolutionStep()
{
    if (q_ > 0)
    {
        forAll(vCurr_, i)
        {
            vPrev_.append(new vectorField(vCurr_[i]));
            wPrev_.append(new vectorField(wCurr_[i]));
        }
    }
}


void coconutModel::add(const vectorField& r, const vectorField& xt)
{
    if (added_)
    {
        // Calculate the difference vectors
        const vectorField dr(r - rRef_);
        const vectorField dxt(xt - xtRef_);

        // Update vCurr_ and wCurr_
        vCurr_.append(new vectorField(dr));
        wCurr_.append(new vectorField(dxt));
    }
    else
    {
        added_ = true;
    }

    // Update reference values
    rRef_ = r;
    xtRef_ = xt;
}


bool coconutModel::isReady()
{
    if (vCurr_.size() + vPrev_.size() > 0)
    {
        return true;
    }

    return false;
}


tmp<vectorField> coconutModel::predict
(
    const vectorField& drVec, const int modes
)
{
    // Prepare result field
    tmp<vectorField> tresult
    (
        new vectorField(drVec.size(), vector::zero)
    );
    vectorField& result = tresult.ref();

    // Convert dr from a vectorField to a scalarField
    const scalarField dr(unroll(drVec));

    // Filter
    // Could we use pivoting instead of filtering??
    Warning
        << "filtering disabled" << endl;
    // filter();

    // Combine and limit vCurr and vPrev
    Info<< "combineAndLimit(v, ...)" << endl;
    scalarRectangularMatrix v;
    combineAndLimit(v, vCurr_, vPrev_, modes);

    // Combine and limit wCurr and wPrev
    Info<< "combineAndLimit(w, ...)" << endl;
    scalarRectangularMatrix w;
    combineAndLimit(w, wCurr_, wPrev_, modes);

    // Check if the matrix has any columns
    if (v.nCols() == 0)
    {
        Info<< "Least-squares model has no information to predict: "
            "zero returned" << endl;

        return tresult;
    }

    // Approximation for the inverse of the Jacobian from a least-squares model 

    // Perform QR decomposition on v
    QRMatrix<scalarRectangularMatrix> qr
    (
        v,
        QRMatrix<scalarRectangularMatrix>::modes::ECONOMY,
        //QRMatrix<scalarRectangularMatrix>::modes::FULL,
        QRMatrix<scalarRectangularMatrix>::outputs::BOTH_QR,
        false // pivoting
    );
    const scalarRectangularMatrix& qq = qr.Q();
    //const scalarRectangularMatrix& rr = qr.R();

    const scalarField b(qq.Tmul(dr));
    //const scalarField c(rr.solve(b));
    const scalarField c(qr.solve(b));

    // Multiply w by c
    const scalarField dxt(w*c);

    // Convert dxt from a scalarField into a vectorField
    result = roll(dxt);

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
