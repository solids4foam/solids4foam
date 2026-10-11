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
#include "compatibilityFunctions.H"

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
            autoPtr<Field<Type>>
            (
                new Field<Type>(size_, pTraits<Type>::zero)
            )
        );
    }

    return table[name]();
}


template<class Type>
const Field<Type>& mechanicalConstitutiveLawState::lookupField
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

mechanicalConstitutiveLawState& mechanicalConstitutiveLawState::child
(
    const word& name
) const
{
    HashTable<autoPtr<mechanicalConstitutiveLawState>>::iterator iter =
        children_.find(name);

    if (iter != children_.end())
    {
        return iter()();
    }

    if (isShadow())
    {
        // A shadow of this state must present shadows of the children, or a
        // sub-law evaluated through it would read and write the parent's own
        // history, which is exactly what shadowing exists to prevent
        children_.insert
        (
            name,
            autoPtr<mechanicalConstitutiveLawState>
            (
                new mechanicalConstitutiveLawState
                (
                    shadowedPtr_->child(name),
                    SHADOW
                )
            )
        );
    }
    else
    {
        children_.insert
        (
            name,
            autoPtr<mechanicalConstitutiveLawState>
            (
                new mechanicalConstitutiveLawState(size_)
            )
        );
    }

    return children_[name]();
}


bool mechanicalConstitutiveLawState::foundChild(const word& name) const
{
    return children_.found(name);
}


template<class Type>
void mechanicalConstitutiveLawState::resizeFields(const label newSize)
{
    forAllIters(fields<Type>(), iter)
    {
        iter()->setSize(newSize, pTraits<Type>::zero);
    }

    forAllIters(fields0<Type>(), iter)
    {
        iter()->setSize(newSize, pTraits<Type>::zero);
    }
}


template<class Type>
void mechanicalConstitutiveLawState::storeOldTimeFields()
{
    HashTable<autoPtr<Field<Type>>>& current = fields<Type>();
    HashTable<autoPtr<Field<Type>>>& old = fields0<Type>();

    forAllConstIters(old, iter)
    {
        const word& name = iter.key();

        if (!current.found(name))
        {
            FatalErrorInFunction
                << "Old-time " << pTraits<Type>::typeName << " state '"
                << name << "' has no corresponding current field."
                << exit(FatalError);
        }

        old[name]() = current[name]();
    }
}


void mechanicalConstitutiveLawState::setSize(const label newSize)
{
    checkNotShadow("setSize");

    size_ = newSize;

    forAllIters(children_, citer)
    {
        citer()->setSize(newSize);
    }

    resizeFields<scalar>(newSize);
    resizeFields<vector>(newSize);
    resizeFields<tensor>(newSize);
    resizeFields<symmTensor>(newSize);
}


void mechanicalConstitutiveLawState::storeOldTime()
{
    checkNotShadow("storeOldTime");

    // A composite's history lives in its children, so the rollover has to
    // reach them; otherwise a sub-law would read this step's values as though
    // they were last step's
    forAllIters(children_, citer)
    {
        citer()->storeOldTime();
    }

    storeOldTimeFields<scalar>();
    storeOldTimeFields<vector>();
    storeOldTimeFields<tensor>();
    storeOldTimeFields<symmTensor>();
}


// * * * * * * * * * * * * * * Typed accessors  * * * * * * * * * * * * * * //

template<class Type>
Field<Type>& mechanicalConstitutiveLawState::field(const word& name)
{
    return accessField(fields<Type>(), name);
}


template<class Type>
Field<Type>& mechanicalConstitutiveLawState::field0(const word& name)
{
    checkNotShadow(word(pTraits<Type>::typeName) + "Field0");

    return accessField(fields0<Type>(), name);
}


template<class Type>
const Field<Type>& mechanicalConstitutiveLawState::getField
(
    const word& name
) const
{
    return lookupField(fields<Type>(), name);
}


template<class Type>
const Field<Type>& mechanicalConstitutiveLawState::getField0
(
    const word& name
) const
{
    return lookupField(readableFields0<Type>(), name);
}


// The four types a state holds; the typed accessors in the header and
// anything written over the field type use these
#define makeStateFieldAccess(Type)                                             \
    template Field<Type>& mechanicalConstitutiveLawState::field<Type>          \
    (                                                                          \
        const word&                                                            \
    );                                                                         \
    template Field<Type>& mechanicalConstitutiveLawState::field0<Type>         \
    (                                                                          \
        const word&                                                            \
    );                                                                         \
    template const Field<Type>&                                                \
    mechanicalConstitutiveLawState::getField<Type>(const word&) const;         \
    template const Field<Type>&                                                \
    mechanicalConstitutiveLawState::getField0<Type>(const word&) const;

makeStateFieldAccess(scalar)
makeStateFieldAccess(vector)
makeStateFieldAccess(tensor)
makeStateFieldAccess(symmTensor)

#undef makeStateFieldAccess

} // End namespace Foam

// ************************************************************************* //
