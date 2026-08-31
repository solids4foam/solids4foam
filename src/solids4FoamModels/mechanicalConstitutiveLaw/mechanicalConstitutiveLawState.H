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


    // Private helper interface

        template<class Type>
        HashTable<autoPtr<Field<Type>>>& fields();

        template<class Type>
        const HashTable<autoPtr<Field<Type>>>& fields() const;

        template<class Type>
        HashTable<autoPtr<Field<Type>>>& fields0();

        template<class Type>
        const HashTable<autoPtr<Field<Type>>>& fields0() const;


        template<class Type>
        Field<Type>& accessField
        (
            HashTable<autoPtr<Field<Type>>>& table,
            const word& name
        );

        template<class Type>
        const Field<Type>& getField
        (
            const HashTable<autoPtr<Field<Type>>>& table,
            const word& name
        ) const;


    // No copy

        mechanicalConstitutiveLawState
        (
            const mechanicalConstitutiveLawState&
        ) = delete;
        void operator=(const mechanicalConstitutiveLawState&) = delete;


public:

    // Constructors

        explicit mechanicalConstitutiveLawState(const label size = 0)
        :
            size_(size)
        {}


    // Access

        label size() const
        {
            return size_;
        }


    // Edit

        void setSize(const label newSize);

        // Commit current state to old-time storage
        void storeOldTime();


    // Scalar fields

        Field<scalar>& scalarField(const word& name);
        const Field<scalar>& scalarField(const word& name) const;

        Field<scalar>& scalarField0(const word& name);
        const Field<scalar>& scalarField0(const word& name) const;

        // Get-only (no creation)
        const Field<scalar>& getScalarField(const word& name) const;
        const Field<scalar>& getScalarField0(const word& name) const;


    // Vector fields

        Field<vector>& vectorField(const word& name);
        const Field<vector>& vectorField(const word& name) const;

        Field<vector>& vectorField0(const word& name);
        const Field<vector>& vectorField0(const word& name) const;

        const Field<vector>& getVectorField(const word& name) const;
        const Field<vector>& getVectorField0(const word& name) const;


    // Tensor fields

        Field<tensor>& tensorField(const word& name);
        const Field<tensor>& tensorField(const word& name) const;

        Field<tensor>& tensorField0(const word& name);
        const Field<tensor>& tensorField0(const word& name) const;

        const Field<tensor>& getTensorField(const word& name) const;
        const Field<tensor>& getTensorField0(const word& name) const;


    // Symmetric tensor fields

        Field<symmTensor>& symmTensorField(const word& name);
        const Field<symmTensor>& symmTensorField(const word& name) const;

        Field<symmTensor>& symmTensorField0(const word& name);
        const Field<symmTensor>& symmTensorField0(const word& name) const;

        const Field<symmTensor>& getSymmTensorField(const word& name) const;
        const Field<symmTensor>& getSymmTensorField0(const word& name) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
