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

Class
    Foam::mechanicalConstitutiveLawState

Description
    Container owning history-dependent mechanical constitutive state variables
    associated with a single mechanical constitutive law instance.

    The state stores per-integration-point variables (e.g. per cell, face,
    vertex, or higher-order integration point) such as plastic strain,
    hardening variables, or internal inelastic measures.

    State variables are stored as typed Field<T> objects indexed by name.
    Fields are created lazily on first access unless explicitly queried
    via get-only accessors.

    Both current-time and old-time storage are supported to enable implicit
    time integration and consistent return mapping algorithms.

    This class contains no constitutive logic and is agnostic to discretisation.
    It provides owned, typed storage only.

SourceFiles
    mechanicalConstitutiveLawState.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawState_H
#define mechanicalConstitutiveLawState_H

#include "mechanicalConstitutiveLawState.H"
#include "UIndirectList.H"
#include "Field.H"
#include "tensor.H"
#include "symmTensor.H"
#include "tmp.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class mechanicalConstitutiveLawState Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawState
{
    // Private data

        //- Current-time fields
        HashTable<autoPtr<Field<scalar>>>     scalarFields_;
        HashTable<autoPtr<Field<vector>>>     vectorFields_;
        HashTable<autoPtr<Field<tensor>>>     tensorFields_;
        HashTable<autoPtr<Field<symmTensor>>> symmTensorFields_;

        //- Old-time fields
        HashTable<autoPtr<Field<scalar>>> scalarFields0_;
        HashTable<autoPtr<Field<vector>>> vectorFields0_;
        HashTable<autoPtr<Field<tensor>>> tensorFields0_;
        HashTable<autoPtr<Field<symmTensor>>> symmTensorFields0_;

        label size_;

        //- Parent state shadowed by this one, or nullptr if this state owns
        //  its history.
        //  A shadow aliases the parent's old-time fields, which are read-only
        //  through it, and owns its own current-time fields
        const mechanicalConstitutiveLawState* shadowedPtr_;

        //- States owned by this one, keyed by name.
        //
        //  A composite law gives each sub-law a state of its own rather than
        //  sharing one and hoping the names do not collide. A sub-law then
        //  looks its history up by its own unqualified name, which is the
        //  point: prefixing a declaration does not prefix the lookups the
        //  sub-law already makes
        mutable HashTable<autoPtr<mechanicalConstitutiveLawState>> children_;


    // Private helper interface

        //- The current-time table for this type. Specialised for scalar,
        //  vector, tensor and symmTensor, so that the public accessors below
        //  are one line each rather than four near-copies
        template<class Type>
        HashTable<autoPtr<Field<Type>>>& fields();

        //- As above, const
        template<class Type>
        const HashTable<autoPtr<Field<Type>>>& fields() const;

        //- The old-time table for this type
        template<class Type>
        HashTable<autoPtr<Field<Type>>>& fields0();

        //- As above, const
        template<class Type>
        const HashTable<autoPtr<Field<Type>>>& fields0() const;

        //- Look a field up, creating it zeroed at this state's size if it is
        //  not there. This is the "declare or use" half of the interface: a
        //  law's first evaluation reaches for a history variable that does
        //  not exist yet and gets one
        template<class Type>
        Field<Type>& accessField
        (
            HashTable<autoPtr<Field<Type>>>& table,
            const word& name
        );

        //- Look a field up, and fail if it is not there. This is the "use
        //  what was declared" half: a law that misspells a name, or reaches
        //  for a variable it never declared, gets an error naming the field
        //  rather than a silent field of zeros that looks like an initial
        //  condition
        template<class Type>
        const Field<Type>& getField
        (
            const HashTable<autoPtr<Field<Type>>>& table,
            const word& name
        ) const;


        //- Fatal error if this state is a shadow, for operations that would
        //  write history
        void checkNotShadow(const word& what) const;

        //- The old-time table to read from: the parent's when shadowing
        template<class Type>
        const HashTable<autoPtr<Field<Type>>>& readableFields0() const;


    // No copy

        mechanicalConstitutiveLawState
        (
            const mechanicalConstitutiveLawState&
        ) = delete;
        void operator=(const mechanicalConstitutiveLawState&) = delete;


public:

    // Public Types

        //- Tag selecting the shadow constructor.
        //  Copy construction stays deleted: a shadow is not a copy, and the
        //  distinction matters enough to be spelled out at the call site
        enum shadowType
        {
            SHADOW
        };


    // Constructors

        explicit mechanicalConstitutiveLawState(const label size = 0)
        :
            size_(size),
            shadowedPtr_(nullptr)
        {}

        //- Construct a shadow of another state.
        //  The shadow aliases the parent's old-time fields, which are
        //  read-only through it, and owns its own current-time fields.
        //  A constitutive law is a function of the kinematics and the old-time
        //  state: it reads only old-time fields and writes only current-time
        //  fields. Evaluating a law into a shadow therefore produces the
        //  correct answer while leaving every one of the parent's fields
        //  untouched, which is exactly what a tangent query needs - including
        //  a finite-difference tangent, which evaluates repeatedly at
        //  perturbed kinematics and must not leave those perturbed values
        //  where endTimeStep() can read them or storeOldTime() can commit them.
        //  History is aliased rather than copied because it is only ever read,
        //  so the cost is one set of current-time fields and no copy of the
        //  parent's history.
        //  The parent must outlive the shadow
        mechanicalConstitutiveLawState
        (
            const mechanicalConstitutiveLawState& parent,
            const shadowType
        )
        :
            size_(parent.size()),
            shadowedPtr_(&parent)
        {}

        //- A shadow of a temporary would dangle as soon as the statement
        //  creating it ended, so make that a compile-time error rather than a
        //  comment asking the caller to be careful
        mechanicalConstitutiveLawState
        (
            const mechanicalConstitutiveLawState&&,
            const shadowType
        ) = delete;


    // Access

        label size() const
        {
            return size_;
        }

        //- Is this state a shadow of another?
        bool isShadow() const
        {
            return shadowedPtr_ != nullptr;
        }


    // Edit

        //- A state owned by this one, created on first use at this state's
        //  size. On a shadow it returns a shadow of the parent's child, so a
        //  tangent query evaluated into a shadow reaches shadowed history all
        //  the way down rather than only at the top level
        mechanicalConstitutiveLawState& child(const word& name) const;

        //- Whether a child of this name exists, without creating one
        bool foundChild(const word& name) const;

        //- Names of the children, in no particular order
        wordList childNames() const;

        void setSize(const label newSize);

        // Commit current state to old-time storage
        void storeOldTime();


    // Scalar fields
    //
    // Two families, and the difference between them is deliberate rather
    // than historical. The vector, tensor and symmTensor blocks below repeat
    // the same pattern, so this is the only place it is spelled out.
    //
    //   xField(name)    creates the field if it is not there
    //   getXField(name) fails if the field is not there
    //
    // A law wants the first when it is declaring what history it keeps -
    // typically on its first evaluation, where a plastic strain or a yield
    // stress has to come into existence - and the second everywhere after,
    // when it means to read back what it stored. Using the creating form to
    // read is the bug the two names exist to prevent: a misspelling, or a
    // field the law forgot to declare, then silently yields zeros that are
    // indistinguishable from a genuine initial condition, and the run
    // produces a plausible wrong answer instead of stopping.
    //
    // The trailing 0 marks old time. Current-time fields are what the law
    // writes; old-time fields are what storeOldTime() committed at the end
    // of the previous step, and are what a history-dependent law integrates
    // from.
    //
    // On a shadow state - see isShadow() - the old-time accessors read the
    // parent's fields, while the current-time ones are the shadow's own. That
    // is what lets a tangent query perturb a law without disturbing the
    // history it will be asked to integrate for real afterwards.

        //- Current-time scalar field, created zeroed if absent
        Field<scalar>& scalarField(const word& name);

        //- As above, const. Still creates: the state is mutable inside
        const Field<scalar>& scalarField(const word& name) const;

        //- Old-time scalar field, created zeroed if absent. On a shadow this
        //  is the parent's field
        Field<scalar>& scalarField0(const word& name);

        //- As above, const
        const Field<scalar>& scalarField0(const word& name) const;

        //- Current-time scalar field, fatal if it was never declared
        const Field<scalar>& getScalarField(const word& name) const;

        //- Old-time scalar field, fatal if it was never declared
        const Field<scalar>& getScalarField0(const word& name) const;


    // Vector fields
    //
    // As the scalar block above: the plain form creates, the get form fails,
    // and the trailing 0 is old time.

        //- Current-time vector field, created zeroed if absent
        Field<vector>& vectorField(const word& name);

        //- As above, const
        const Field<vector>& vectorField(const word& name) const;

        //- Old-time vector field, created zeroed if absent. On a shadow this is
        //  the parent's field
        Field<vector>& vectorField0(const word& name);

        //- As above, const
        const Field<vector>& vectorField0(const word& name) const;

        //- Current-time vector field, fatal if it was never declared
        const Field<vector>& getVectorField(const word& name) const;

        //- Old-time vector field, fatal if it was never declared
        const Field<vector>& getVectorField0(const word& name) const;


    // Tensor fields
    //
    // As the scalar block above: the plain form creates, the get form fails,
    // and the trailing 0 is old time.

        //- Current-time tensor field, created zeroed if absent
        Field<tensor>& tensorField(const word& name);

        //- As above, const
        const Field<tensor>& tensorField(const word& name) const;

        //- Old-time tensor field, created zeroed if absent. On a shadow this is
        //  the parent's field
        Field<tensor>& tensorField0(const word& name);

        //- As above, const
        const Field<tensor>& tensorField0(const word& name) const;

        //- Current-time tensor field, fatal if it was never declared
        const Field<tensor>& getTensorField(const word& name) const;

        //- Old-time tensor field, fatal if it was never declared
        const Field<tensor>& getTensorField0(const word& name) const;


    // Symmetric tensor fields
    //
    // As the scalar block above: the plain form creates, the get form fails,
    // and the trailing 0 is old time.

        //- Current-time symmTensor field, created zeroed if absent
        Field<symmTensor>& symmTensorField(const word& name);

        //- As above, const
        const Field<symmTensor>& symmTensorField(const word& name) const;

        //- Old-time symmTensor field, created zeroed if absent. On a shadow this is
        //  the parent's field
        Field<symmTensor>& symmTensorField0(const word& name);

        //- As above, const
        const Field<symmTensor>& symmTensorField0(const word& name) const;

        //- Current-time symmTensor field, fatal if it was never declared
        const Field<symmTensor>& getSymmTensorField(const word& name) const;

        //- Old-time symmTensor field, fatal if it was never declared
        const Field<symmTensor>& getSymmTensorField0(const word& name) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
