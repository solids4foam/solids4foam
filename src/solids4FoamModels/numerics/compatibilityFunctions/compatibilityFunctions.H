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
    Foam::compatibilityFunctions

Description
      Collection of small helper functions that provide a uniform API across
      different OpenFOAM variants (foam-extend, OpenFOAM.org, OpenFOAM.com).

      These wrappers centralise differences in class and function interfaces,
      allowing the main codebase to remain free of scattered #ifdefs. Typical
      helpers include:

        - tmpRef()           : access the underlying object of a tmp<T>
        - autoPtrRef()       : access the underlying object of a autoPtr<T>
        - boundaryFieldRef() : obtain a modifiable boundary field reference
        - makeList<T>()      : construct a List<T> from an initializer_list

      Usage example:
          #include "compatibilityFunctions.H"

          auto& phi = Foam::tmpRef(tphi);
          boundaryRef(dualTraction)[patchI] = vector::zero;
          labelList ids = Foam::makeList<label>({0,1,2});

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef compatibilityFunctions_H
#define compatibilityFunctions_H

#include "tmp.H"
#include "autoPtr.H"
#include "GeometricFields.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // non-const
    template<class T>
    inline T& tmpRef(tmp<T>& t)
    {
#ifdef FOAMEXTEND
        return t();     // foam-extend
#else
        return t.ref(); // foundation/esi
#endif
    }

    // const
    template<class T>
    inline const T& tmpRef(const tmp<T>& t)
    {
#ifdef FOAMEXTEND
        return t();
#else
        return t.ref();
#endif
    }

    // non-const
    template<class T>
    inline T& autoPtrRef(autoPtr<T>& t)
    {
#ifdef OPENFOAM_COM
        return t.ref(); // com
#else
        return t();     // foam-extend/org
#endif
    }

    // const
    template<class T>
    inline const T& autoPtrRef(const autoPtr<T>& t)
    {
#ifdef OPENFOAM_COM
        return t.ref(); // com
#else
        return t();     // foam-extend/org
#endif
    }

    // Non-const field → modifiable primitiveField&
    template<class Type, template<class> class PatchField, class GeoMesh>
    inline Foam::Field<Type>& primitiveFieldRef
    (
        Foam::GeometricField<Type, PatchField, GeoMesh>& fld
    )
    {
#ifdef OPENFOAM_NOT_EXTEND
        return fld.primitiveFieldRef();
#else
        return fld.internalField();
#endif
    }

    // Non-const field → modifiable Boundary&
    template<class Type, template<class> class PatchField, class GeoMesh>
#ifdef FOAMEXTEND
    inline typename GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField&
#else
    inline typename GeometricField<Type, PatchField, GeoMesh>::Boundary&
#endif
    boundaryFieldRef(GeometricField<Type, PatchField, GeoMesh>& fld)
    {
#ifdef OPENFOAM_NOT_EXTEND
        return fld.boundaryFieldRef();
#else
        return fld.boundaryField();
#endif
    }

    // Provides a helper function Foam::makeList<T>() for constructing
    // Foam::List<T> objects from an initializer_list.
    // This function enables code such as
    //     labelList ids = Foam::makeList<label>({0,1,2});
    // which is supported in modern OpenFOAM via public initializer_list
    // constructors, but not in foam-extend. Using makeList() centralises
    // the difference and provides a single, portable syntax across all
    // OpenFOAM variants.
    template<class T>
    inline List<T> makeList(std::initializer_list<T> il)
    {
        List<T> L(il.size());
        label i = 0;

        for (const T& v : il)
        {
            L[i++] = v;
        }

        return L; // NRVO/move on modern forks; fine for small lists regardless
    }

    // fvm::ddt(const volVectorField&) which is missing in foam-extend
    inline Foam::tmp<Foam::fvMatrix<Foam::vector>>
    fvmDdtVectorCompat(const Foam::volVectorField& D)
    {
    #ifdef FOAMEXTEND
        return Foam::fvm::ddt(const_cast<Foam::volVectorField&>(D));
    #else
        return Foam::fvm::ddt(D);
    #endif
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
