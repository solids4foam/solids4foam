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

#include "sparseMatrix.H"
#include <vector>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sparseMatrix, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sparseMatrix::sparseMatrix(const label size)
:
    refCount(),
    data_(size)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::sparseMatrix::nBlockRows() const
{
    label nBlockRows = 0;

    for (auto iter = data_.begin(); iter != data_.end(); ++iter)
    {
        const label rowI = iter.key()[0];

        nBlockRows = max(nBlockRows, rowI);
    }

    return nBlockRows + 1;
}


Foam::tensor& Foam::sparseMatrix::operator()
(
    const label rowI,
    const label colI
)
{
    // Create key for the entry
    FixedList<label, 2> key;
    key[0] = rowI;
    key[1] = colI;

    // Return a reference to the entry
    // If it does not exist then it will be initialised to zero first
    sparseMatrixData::iterator iter = data_.find(key);

    if (iter == data_.end())
    {
        data_.insert(key, tensor::zero);
        return *(data_.find(key));
    }
    else
    {
        return *iter;
    }
}


void Foam::sparseMatrix::print() const
{
    Info<< "void Foam::sparseMatrix::print() const" << endl;

    // Create a vector to store the matrix indices
    std::vector<FixedList<label, 2>> keys(data_.size());

    // Insert the matrix indices into the vector for all data
    int i = 0;
    for (auto iter = data_.begin(); iter != data_.end(); ++iter)
    {
        keys[i] = iter.key();
        i++;
    }

    // Define custom sorting criteria
    auto cmp = [](const FixedList<label, 2>& a, const FixedList<label, 2>& b)
    {
        if (a[0] < b[0])
        {
            return true;
        }
        else if (a[0] == b[0])
        {
            return a[1] < b[1];
        }
        else
        {
            return false;
        }
    };

    // Sort keys by row and column
    std::sort(keys.begin(), keys.end(), cmp);

    // Print out sorted values
    for (unsigned long k = 0; k < keys.size(); ++k)
    {
        const label rowI = keys[k][0];
        const label colI = keys[k][1];
        const tensor& coeff = data_[keys[k]];

        Info << "(" << rowI << ", " << colI << ") : " << coeff << endl;
    }
}


void Foam::sparseMatrix::operator+=
(
    const sparseMatrix& A
)
{
    // Copy the A data into the result one-by-one
    const sparseMatrixData& AData = A.data();
    for (auto AIter = AData.begin(); AIter != AData.end(); ++AIter)
    {
        // Check if the entry already exists
        sparseMatrixData::iterator iter =
            data_.find(AIter.key());

        if (iter == data_.end())
        {
            // Create a new entry
            data_.insert(AIter.key(), AIter.val());
        }
        else
        {
            // Add to the existing entry
            iter.val() += AIter.val();
        }
    }
}


void Foam::sparseMatrix::operator-=
(
    const sparseMatrix& A
)
{
    // Copy the A data into the result one-by-one
    const sparseMatrixData& AData = A.data();
    for (auto AIter = AData.begin(); AIter != AData.end(); ++AIter)
    {
        // Check if the entry already exists
        sparseMatrixData::iterator iter =
            data_.find(AIter.key());

        if (iter == data_.end())
        {
            // Create a new entry
            data_.insert(AIter.key(), -AIter.val());
        }
        else
        {
            // Add to the existing entry
            iter.val() -= AIter.val();
        }
    }
}


Foam::tmp<Foam::sparseMatrix> Foam::sparseMatrix::operator+
(
    const sparseMatrix& A
) const
{
    // Prepare the result
    tmp<sparseMatrix> tresult
    (
        new sparseMatrix(max(data_.size(), A.data().size()))
    );
    sparseMatrix& result = tresult.ref();

    // Copy data_ into the result
    result.data() = data_;

    // Copy the A data into the result one-by-one
    sparseMatrixData& resultData = result.data();
    const sparseMatrixData& AData = A.data();
    for (auto AIter = AData.begin(); AIter != AData.end(); ++AIter)
    {
        // Check if the entry already exists
        sparseMatrixData::iterator resultIter =
            resultData.find(AIter.key());

        if (resultIter == data_.end())
        {
            // Create a new entry
            resultData.insert(AIter.key(), AIter.val());
        }
        else
        {
            // Add to the existing entry
            resultIter.val() += AIter.val();
        }
    }

    return tresult;
}


Foam::tmp<Foam::sparseMatrix> Foam::sparseMatrix::operator*
(
    const scalarField& sf
) const
{
    // Prepare the result
    tmp<sparseMatrix> tresult(new sparseMatrix(data_.size()));
    sparseMatrix& result = tresult.ref();

    // Copy data_ into the result
    result.data() = data_;

    // Copy data_ into the result
    result.data() = data_;

    // Multiply each row by sf[rowI]
    sparseMatrixData& data = result.data();
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        const label rowI = iter.key()[0];
        iter.val() *= sf[rowI];
    }

    return tresult;
}


Foam::tmp<Foam::sparseMatrix> Foam::sparseMatrix::operator*
(
    const scalar s
) const
{
    // Prepare the result
    tmp<sparseMatrix> tresult(new sparseMatrix(data_.size()));
    sparseMatrix& result = tresult.ref();

    // Copy data_ into the result
    result.data() = data_;

    // Multiply each entry by s
    sparseMatrixData& data = result.data();
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        iter.val() *= s;
    }

    return tresult;
}


// ************************************************************************* //
