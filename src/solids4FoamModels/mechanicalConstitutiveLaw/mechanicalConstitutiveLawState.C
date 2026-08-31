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

#include "mechanicalConstitutiveLawState.H"

namespace Foam
{

// * * * * * * * * * * * * Private helper tables * * * * * * * * * * * * * //

template<>
HashTable<autoPtr<Field<scalar>>>&
mechanicalConstitutiveLawState::fields<scalar>()
{
    return scalarFields_;
}

template<>
const HashTable<autoPtr<Field<scalar>>>&
mechanicalConstitutiveLawState::fields<scalar>() const
{
    return scalarFields_;
}

template<>
HashTable<autoPtr<Field<scalar>>>&
mechanicalConstitutiveLawState::fields0<scalar>()
{
    return scalarFields0_;
}

template<>
const HashTable<autoPtr<Field<scalar>>>&
mechanicalConstitutiveLawState::fields0<scalar>() const
{
    return scalarFields0_;
}


template<>
HashTable<autoPtr<Field<vector>>>&
mechanicalConstitutiveLawState::fields<vector>()
{
    return vectorFields_;
}

template<>
const HashTable<autoPtr<Field<vector>>>&
mechanicalConstitutiveLawState::fields<vector>() const
{
    return vectorFields_;
}

template<>
HashTable<autoPtr<Field<vector>>>&
mechanicalConstitutiveLawState::fields0<vector>()
{
    return vectorFields0_;
}

template<>
const HashTable<autoPtr<Field<vector>>>&
mechanicalConstitutiveLawState::fields0<vector>() const
{
    return vectorFields0_;
}


template<>
HashTable<autoPtr<Field<tensor>>>&
mechanicalConstitutiveLawState::fields<tensor>()
{
    return tensorFields_;
}

template<>
const HashTable<autoPtr<Field<tensor>>>&
mechanicalConstitutiveLawState::fields<tensor>() const
{
    return tensorFields_;
}

template<>
HashTable<autoPtr<Field<tensor>>>&
mechanicalConstitutiveLawState::fields0<tensor>()
{
    return tensorFields0_;
}

template<>
const HashTable<autoPtr<Field<tensor>>>&
mechanicalConstitutiveLawState::fields0<tensor>() const
{
    return tensorFields0_;
}


template<>
HashTable<autoPtr<Field<symmTensor>>>&
mechanicalConstitutiveLawState::fields<symmTensor>()
{
    return symmTensorFields_;
}

template<>
const HashTable<autoPtr<Field<symmTensor>>>&
mechanicalConstitutiveLawState::fields<symmTensor>() const
{
    return symmTensorFields_;
}

template<>
HashTable<autoPtr<Field<symmTensor>>>&
mechanicalConstitutiveLawState::fields0<symmTensor>()
{
    return symmTensorFields0_;
}

template<>
const HashTable<autoPtr<Field<symmTensor>>>&
mechanicalConstitutiveLawState::fields0<symmTensor>() const
{
    return symmTensorFields0_;
}


// * * * * * * * * * * * * Field access helpers * * * * * * * * * * * * * //

template<class Type>
Field<Type>& mechanicalConstitutiveLawState::accessField
(
    HashTable<autoPtr<Field<Type>>>& table,
    const word& name
)
{
    if (!table.found(name))
    {
        table.insert
        (
            name,
            autoPtr<Field<Type>>(new Field<Type>(size_, Zero))
        );
    }

    return table[name]();
}


template<class Type>
const Field<Type>& mechanicalConstitutiveLawState::getField
(
    const HashTable<autoPtr<Field<Type>>>& table,
    const word& name
) const
{
    if (!table.found(name))
    {
        FatalErrorInFunction
            << "Requested state field '" << name
            << "' does not exist."
            << exit(FatalError);
    }

    return table[name]();
}


void mechanicalConstitutiveLawState::checkNotShadow(const word& what) const
{
    if (isShadow())
    {
        FatalErrorInFunction
            << "'" << what << "' would modify history through a shadow state."
            << nl
            << "A shadow aliases the old-time fields of its parent so that a "
            << "tangent query can evaluate a law without disturbing them. "
            << "Only current-time fields may be written through a shadow."
            << exit(FatalError);
    }
}


template<class Type>
const HashTable<autoPtr<Field<Type>>>&
mechanicalConstitutiveLawState::readableFields0() const
{
    // A shadow reads its parent's history, never its own
    if (isShadow())
    {
        return shadowedPtr_->readableFields0<Type>();
    }

    return fields0<Type>();
}


// * * * * * * * * * * * * Public interface * * * * * * * * * * * * * * * //

void mechanicalConstitutiveLawState::setSize(const label newSize)
{
    checkNotShadow("setSize");

    size_ = newSize;

    forAllIter(HashTable<autoPtr<Field<scalar>>>, scalarFields_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
    forAllIter(HashTable<autoPtr<Field<vector>>>, vectorFields_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
    forAllIter(HashTable<autoPtr<Field<tensor>>>, tensorFields_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
    forAllIter(HashTable<autoPtr<Field<symmTensor>>>, symmTensorFields_, iter)
    {
        iter()->setSize(newSize, Zero);
    }

    forAllIter(HashTable<autoPtr<Field<scalar>>>, scalarFields0_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
    forAllIter(HashTable<autoPtr<Field<vector>>>, vectorFields0_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
    forAllIter(HashTable<autoPtr<Field<tensor>>>, tensorFields0_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
    forAllIter(HashTable<autoPtr<Field<symmTensor>>>, symmTensorFields0_, iter)
    {
        iter()->setSize(newSize, Zero);
    }
}


void mechanicalConstitutiveLawState::storeOldTime()
{
    checkNotShadow("storeOldTime");

    // Scalars
#ifdef OPENFOAM_COM
    forAllConstIters(scalarFields0_, iter)
#else
    forAllConstIter(HashTable<autoPtr<Field<scalar>>>, scalarFields0_, iter)
#endif
    {
        const word& name = iter.key();
        if (!scalarFields_.found(name))
        {
            FatalErrorInFunction
                << "Old-time scalar state '" << name
                << "' has no corresponding current field."
                << exit(FatalError);
        }
        scalarFields0_[name]() = scalarFields_[name]();
    }

    // Vectors
#ifdef OPENFOAM_COM
    forAllConstIters(vectorFields0_, iter)
#else
    forAllConstIter(HashTable<autoPtr<Field<vector>>>, vectorFields0_, iter)
#endif
    {
        const word& name = iter.key();
        if (!vectorFields_.found(name))
        {
            FatalErrorInFunction
                << "Old-time vector state '" << name
                << "' has no corresponding current field."
                << exit(FatalError);
        }
        vectorFields0_[name]() = vectorFields_[name]();
    }

    // Tensors
#ifdef OPENFOAM_COM
    forAllConstIters(tensorFields0_, iter)
#else
    forAllConstIter(HashTable<autoPtr<Field<tensor>>>, tensorFields0_, iter)
#endif
    {
        const word& name = iter.key();
        if (!tensorFields_.found(name))
        {
            FatalErrorInFunction
                << "Old-time tensor state '" << name
                << "' has no corresponding current field."
                << exit(FatalError);
        }
        tensorFields0_[name]() = tensorFields_[name]();
    }

    // SymmTensors
#ifdef OPENFOAM_COM
    forAllConstIters(symmTensorFields0_, iter)
#else
    forAllConstIter
    (
        HashTable<autoPtr<Field<symmTensor>>>, symmTensorFields0_, iter
    )
#endif
    {
        const word& name = iter.key();
        if (!symmTensorFields_.found(name))
        {
            FatalErrorInFunction
                << "Old-time symmTensor state '" << name
                << "' has no corresponding current field."
                << exit(FatalError);
        }
        symmTensorFields0_[name]() = symmTensorFields_[name]();
    }
}


// * * * * * * * * * * * * Scalar accessors * * * * * * * * * * * * * //

Field<scalar>& mechanicalConstitutiveLawState::scalarField(const word& name)
{
    return accessField(fields<scalar>(), name);
}

const Field<scalar>& mechanicalConstitutiveLawState::scalarField
(
    const word& name
) const
{
    return getField(fields<scalar>(), name);
}

Field<scalar>& mechanicalConstitutiveLawState::scalarField0(const word& name)
{
    checkNotShadow("scalarField0");

    return accessField(fields0<scalar>(), name);
}

const Field<scalar>& mechanicalConstitutiveLawState::scalarField0
(
    const word& name
) const
{
    return getField(readableFields0<scalar>(), name);
}

const Field<scalar>& mechanicalConstitutiveLawState::getScalarField
(
    const word& name
) const
{
    return getField(fields<scalar>(), name);
}

const Field<scalar>& mechanicalConstitutiveLawState::getScalarField0
(
    const word& name
) const
{
    return getField(readableFields0<scalar>(), name);
}


// * * * * * * * * * * * * Vector accessors * * * * * * * * * * * * * //

Field<vector>& mechanicalConstitutiveLawState::vectorField(const word& name)
{
    return accessField(fields<vector>(), name);
}

const Field<vector>& mechanicalConstitutiveLawState::vectorField
(
    const word& name
) const
{
    return getField(fields<vector>(), name);
}

Field<vector>& mechanicalConstitutiveLawState::vectorField0(const word& name)
{
    checkNotShadow("vectorField0");

    return accessField(fields0<vector>(), name);
}

const Field<vector>& mechanicalConstitutiveLawState::vectorField0
(
    const word& name
) const
{
    return getField(readableFields0<vector>(), name);
}

const Field<vector>& mechanicalConstitutiveLawState::getVectorField
(
    const word& name
) const
{
    return getField(fields<vector>(), name);
}

const Field<vector>& mechanicalConstitutiveLawState::getVectorField0
(
    const word& name
) const
{
    return getField(readableFields0<vector>(), name);
}


// * * * * * * * * * * * * Tensor accessors * * * * * * * * * * * * * //

Field<tensor>& mechanicalConstitutiveLawState::tensorField(const word& name)
{
    return accessField(fields<tensor>(), name);
}

const Field<tensor>& mechanicalConstitutiveLawState::tensorField
(
    const word& name
) const
{
    return getField(fields<tensor>(), name);
}

Field<tensor>& mechanicalConstitutiveLawState::tensorField0(const word& name)
{
    checkNotShadow("tensorField0");

    return accessField(fields0<tensor>(), name);
}

const Field<tensor>& mechanicalConstitutiveLawState::tensorField0
(
    const word& name
) const
{
    return getField(readableFields0<tensor>(), name);
}

const Field<tensor>& mechanicalConstitutiveLawState::getTensorField
(
    const word& name
) const
{
    return getField(fields<tensor>(), name);
}

const Field<tensor>& mechanicalConstitutiveLawState::getTensorField0
(
    const word& name
) const
{
    return getField(readableFields0<tensor>(), name);
}


// * * * * * * * * * * * * SymmTensor accessors * * * * * * * * * * * * //

Field<symmTensor>& mechanicalConstitutiveLawState::symmTensorField
(
    const word& name
)
{
    return accessField(fields<symmTensor>(), name);
}

const Field<symmTensor>& mechanicalConstitutiveLawState::symmTensorField
(
    const word& name
) const
{
    return getField(fields<symmTensor>(), name);
}

Field<symmTensor>& mechanicalConstitutiveLawState::symmTensorField0
(
    const word& name
)
{
    checkNotShadow("symmTensorField0");

    return accessField(fields0<symmTensor>(), name);
}

const Field<symmTensor>& mechanicalConstitutiveLawState::symmTensorField0
(
    const word& name
) const
{
    return getField(readableFields0<symmTensor>(), name);
}

const Field<symmTensor>& mechanicalConstitutiveLawState::getSymmTensorField
(
    const word& name
) const
{
    return getField(fields<symmTensor>(), name);
}

const Field<symmTensor>& mechanicalConstitutiveLawState::getSymmTensorField0
(
    const word& name
) const
{
    return getField(readableFields0<symmTensor>(), name);
}

} // End namespace Foam

// ************************************************************************* //
