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
    Container class owning the mechanical constitutive state variables
    associated with a single mechanical constitutive law instance.

    The mechanicalConstitutiveLawState stores history-dependent variables
    required by mechanical constitutive laws, such as plastic strains,
    internal hardening variables, or inelastic deformation measures.
    These variables are owned by the state object and persist across
    time steps.

    State variables are stored as typed Field objects (e.g. scalar,
    vector, tensor, symmTensor), indexed by name. Fields are created
    on demand and resized consistently when the state size changes.
    This allows constitutive laws to declare and access only the state
    variables they require.

    The state object is independent of mesh topology and discretisation
    details. It stores data per integration point (e.g. per cell or
    per face or per vertex), but does not itself define how these points map
    to the computational mesh.

    This class contains no constitutive logic and performs no updates
    of state variables. Its sole responsibility is to provide owned,
    typed storage for mechanical constitutive history variables.

    Usage:
    A mechanicalConstitutiveLawState is created and owned by the
    mechanicalConstitutiveLawManager. During constitutive evaluation,
    non-owning views of the state are passed to a mechanical
    constitutive law via a mechanicalConstitutiveLawStateView.

    Notes:
    - The state object owns all stored fields and is responsible for
      their lifetime.
    - State fields are created lazily when first requested by a
      constitutive law.
    - The class is formulation-agnostic and may be used for linear,
      nonlinear, elastic, or inelastic mechanical laws.
    - Input/output and restart functionality may be layered on top of
      this class without changing its core responsibility.

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
    // Private Data

        // Storage for state fields, indexed by name.
        // Each entry owns a Field<T> of size nCells (per material).
        HashTable<autoPtr<Field<scalar>>> scalarFields_;
        HashTable<autoPtr<Field<vector>>> vectorFields_;
        HashTable<autoPtr<Field<tensor>>> tensorFields_;
        HashTable<autoPtr<Field<symmTensor>>> symmTensorFields_;

        // Storage for old-time state fields, indexed by name.
        // Each entry owns a Field<T> of size nCells (per material).
        HashTable<autoPtr<Field<scalar>>> scalarFields0_;
        HashTable<autoPtr<Field<vector>>> vectorFields0_;
        HashTable<autoPtr<Field<tensor>>> tensorFields0_;
        HashTable<autoPtr<Field<symmTensor>>> symmTensorFields0_;

        // Number of entries (e.g. cells) for this state
        label size_;

    // Private Member Functions

        //- No copy construct
        mechanicalConstitutiveLawState
        (
            const mechanicalConstitutiveLawState&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const mechanicalConstitutiveLawState&
        ) = delete;


public:


    // Constructors

        //- Construct empty state with given size
        explicit mechanicalConstitutiveLawState(const label size = 0)
        :
            size_(size)
        {}


    //- Destructor
    ~mechanicalConstitutiveLawState()
    {}


    // Member Functions

    // Access

        //- Number of entries
        label size() const
        {
            return size_;
        }

    // Edit

        //- Resize all stored state fields
        void setSize(const label newSize)
        {
            size_ = newSize;

            // Current-time fields

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

            // Old-time fields

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

        // Copy current state into old-time storage (commit)
        void storeOldTime()
        {
            // Check all old fields have a corresponding current field

            // Scalars
            forAllConstIters(scalarFields0_, iter)
            {
                const word& name = iter.key();

                if (!scalarFields_.found(name))
                {
                    FatalErrorInFunction
                        << "Old-time scalar state field '" << name
                        << "' exists but no corresponding current field is "
                        << "present. This indicates inconsistent state "
                        << "initialisation."
                        << exit(FatalError);
                }
            }

            // Vectors
            forAllConstIters(vectorFields0_, iter)
            {
                const word& name = iter.key();

                if (!vectorFields_.found(name))
                {
                    FatalErrorInFunction
                        << "Old-time vector state field '" << name
                        << "' exists but no corresponding current field is "
                        << "present. This indicates inconsistent state "
                        << "initialisation."
                        << exit(FatalError);
                }
            }

            // Tensors
            forAllConstIters(tensorFields0_, iter)
            {
                const word& name = iter.key();

                if (!tensorFields_.found(name))
                {
                    FatalErrorInFunction
                        << "Old-time tensor state field '" << name
                        << "' exists but no corresponding current field is "
                        << "present. This indicates inconsistent state "
                        << "initialisation."
                        << exit(FatalError);
                }
            }

            // SymmTensors
            forAllConstIters(symmTensorFields0_, iter)
            {
                const word& name = iter.key();

                if (!symmTensorFields_.found(name))
                {
                    FatalErrorInFunction
                        << "Old-time symmTensor state field '" << name
                        << "' exists but no corresponding current field is "
                        << "present. This indicates inconsistent state "
                        << "initialisation."
                        << exit(FatalError);
                }
            }

            // Update old time fields

            // Scalars
            forAllConstIters(scalarFields_, iter)
            {
                const word& name = iter.key();
                const Field<scalar>& cur = iter.val()();

                if (scalarFields0_.found(name))
                {
                    scalarFields0_[name]() = cur;
                }
            }

            // Vectors
            forAllConstIters(vectorFields_, iter)
            {
                const word& name = iter.key();
                const Field<vector>& cur = iter.val()();

                if (vectorFields0_.found(name))
                {
                    vectorFields0_[name]() = cur;
                }
            }

            // Tensors
            forAllConstIters(tensorFields_, iter)
            {
                const word& name = iter.key();
                const Field<tensor>& cur = iter.val()();

                if (tensorFields0_.found(name))
                {
                    tensorFields0_[name]() = cur;
                }
            }

            // SymmTensors
            forAllConstIters(symmTensorFields_, iter)
            {
                const word& name = iter.key();
                const Field<symmTensor>& cur = iter.val()();

                if (symmTensorFields0_.found(name))
                {
                    symmTensorFields0_[name]() = cur;
                }
            }
        }

        //- Access or create scalar state field
        Field<scalar>& scalarField(const word& name)
        {
            if (!scalarFields_.found(name))
            {
                scalarFields_.insert
                (
                    name,
                    autoPtr<Field<scalar>>(new Field<scalar>(size_, Zero))
                );
            }

            return scalarFields_[name]();
        }

        //- Access or create scalar old-time state field
        Field<scalar>& scalarField0(const word& name)
        {
            if (!scalarFields0_.found(name))
            {
                scalarFields0_.insert
                (
                    name,
                    autoPtr<Field<scalar>>(new Field<scalar>(size_, Zero))
                );
            }

            return scalarFields0_[name]();
        }

        //- Access or create vector state field
        Field<vector>& vectorField(const word& name)
        {
            if (!vectorFields_.found(name))
            {
                vectorFields_.insert
                (
                    name,
                    autoPtr<Field<vector>>(new Field<vector>(size_, Zero))
                );
            }

            return vectorFields_[name]();
        }

        //- Access or create vector old-time state field
        Field<vector>& vectorField0(const word& name)
        {
            if (!vectorFields0_.found(name))
            {
                vectorFields0_.insert
                (
                    name,
                    autoPtr<Field<vector>>(new Field<vector>(size_, Zero))
                );
            }

            return vectorFields0_[name]();
        }

        //- Access or create tensor state field
        Field<tensor>& tensorField(const word& name)
        {
            if (!tensorFields_.found(name))
            {
                tensorFields_.insert
                (
                    name,
                    autoPtr<Field<tensor>>(new Field<tensor>(size_, Zero))
                );
            }

            return tensorFields_[name]();
        }

        //- Access or create tensor old-time state field
        Field<tensor>& tensorField0(const word& name)
        {
            if (!tensorFields0_.found(name))
            {
                tensorFields0_.insert
                (
                    name,
                    autoPtr<Field<tensor>>(new Field<tensor>(size_, Zero))
                );
            }

            return tensorFields0_[name]();
        }

        //- Access or create symmTensor state field
        Field<symmTensor>& symmTensorField(const word& name)
        {
            if (!symmTensorFields_.found(name))
            {
                symmTensorFields_.insert
                (
                    name,
                    autoPtr<Field<symmTensor>>(new Field<symmTensor>(size_, Zero))
                );
            }

            return symmTensorFields_[name]();
        }

        //- Access or create scalar old-time state field
        Field<symmTensor>& symmTensorField0(const word& name)
        {
            if (!symmTensorFields0_.found(name))
            {
                symmTensorFields0_.insert
                (
                    name,
                    autoPtr<Field<symmTensor>>(new Field<symmTensor>(size_, Zero))
                );
            }

            return symmTensorFields0_[name]();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
