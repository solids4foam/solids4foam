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

#ifdef USE_PETSC

#include "foamPetscSnesHelper.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class Type>
void foamPetscSnesHelper::extractFieldComponents
(
    const PetscScalar *x,
    Field<Type>& vf,
    const label xBlockSize,
    const label offset
) const
{
    // Obtain a pointer to the underlying scalar array in vf
    // We assume vf is stored as contiguous data, i.e.
    // (vx1, vy1, vz1, vx2, vy2, vz2, ..., vxN, vyN, vzN)
    scalar* vfPtr = reinterpret_cast<scalar*>(vf.begin());

    // Get the number of components per field element from the traits, e.g. 1
    // for a scalar, 3 for a vector, 6 for a symmTensor
    const label vfBlockSize = pTraits<Type>::nComponents;

    // Loop over each element in vf and copy the appropriate components from x
    forAll(vf, blockI)
    {
        for (label cmptI = 0; cmptI < vfBlockSize; ++cmptI)
        {
            vfPtr[blockI*vfBlockSize + cmptI] = x[offset + blockI*xBlockSize + cmptI];
        }
    }
}


template <class Type>
void foamPetscSnesHelper::insertFieldComponents
(
    const Field<Type>& vf,
    PetscScalar *x,
    const label xBlockSize,
    const label offset
) const
{
    // Obtain a pointer to the underlying scalar data in vf (assumes contiguous
    // storage)
    const scalar* vfPtr = reinterpret_cast<const scalar*>(vf.begin());

    // Get the number of components per field element from the traits (e.g., 1
    // for scalar, 3 for vector, etc.)
    const label vfBlockSize = pTraits<Type>::nComponents;

    // Loop over each element in vf and copy its components into the corresponding position in x
    forAll(vf, blockI)
    {
        for (label cmptI = 0; cmptI < vfBlockSize; ++cmptI)
        {
            x[offset + blockI*xBlockSize + cmptI] = vfPtr[blockI*vfBlockSize + cmptI];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
