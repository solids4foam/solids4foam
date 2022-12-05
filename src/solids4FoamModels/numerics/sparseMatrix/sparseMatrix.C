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

#include "sparseMatrix.H"

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



// ************************************************************************* //
