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

#include "SparseMatrixTemplate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(SparseMatrixTemplate, 0);
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::SparseMatrixTemplate<Type>::SparseMatrixTemplate(const label size)
:
    refCount(),
    data_(size)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::SparseMatrixTemplate<Type>::nBlockRows() const
{
    label nBlockRows = 0;

    for (auto iter = data_.begin(); iter != data_.end(); ++iter)
    {
        const label rowI = iter.key()[0];

        nBlockRows = max(nBlockRows, rowI);
    }

    return nBlockRows + 1;
}


template<class Type>
void Foam::SparseMatrixTemplate<Type>::print() const
{
    Info<< "void Foam::SparseMatrixTemplate::print() const" << endl;

    for (auto iter = data_.begin(); iter != data_.end(); ++iter)
    {
        const label rowI = iter.key()[0];
        const label colI = iter.key()[1];
        const Type& val = *iter;

        Info<< "(" << rowI << ", " << colI << "): " << val << endl;
    }
}


template<class Type>
Type& Foam::SparseMatrixTemplate<Type>::operator()
(
    const label rowI,
    const label colI
)
{
    // Create key for the entry
    FixedList<label, 2> key;
    key[0] = rowI;
    key[1] = colI;

    FatalError
        << "stop" << abort(FatalError);
    // Return a reference to the entry
    // If it does not exist then it will be initialised to zero first
    typename SparseMatrixTemplateData::iterator iter = data_.find(key);

    if (iter == data_.end())
    {
        data_.insert(key, pTraits<Type>::zero);
        return *(data_.find(key));
    }
    else
    {
        return *iter;
    }
}



// ************************************************************************* //
