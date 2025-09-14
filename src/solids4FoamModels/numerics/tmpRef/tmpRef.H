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
    Foam::tmpRef

Description
      Provides portable helper functions Foam::tmpRef() and Foam::autoRef()
      for accessing the contents of Foam::tmp<T> and Foam::autoPtr<T> objects
      across different OpenFOAM variants.

      - In OpenFOAM-foundation/ESI forks, tmp<T>::ref() is used.
      - In foam-extend, operator() is used instead.

      Using tmpRef() and autoPtrRef() removes the need for scattered #ifdefs in
      user code:
          tmp<volScalarField> tf = fvc::div(phi);
          volScalarField& f = Foam::tmpRef(tf);

          autoPtr<scalarField> apf(new scalarField(10, 0.0));
          volScalarField& f = Foam::autoPtrRef(apf);

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef tmpRef_H
#define tmpRef_H

#include "tmp.H"
#include "autoPtr.H"

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
#ifdef FOAMEXTEND
        return t();     // foam-extend
#else
        return t.ref(); // foundation/esi
#endif
    }

    // const
    template<class T>
    inline const T& autoPtrRef(const autoPtr<T>& t)
    {
#ifdef FOAMEXTEND
        return t();
#else
        return t.ref();
#endif
    }

    // Non-const field → modifiable Boundary&
    template<class Type, template<class> class PatchField, class GeoMesh>
#ifdef FOAMEXTEND
    inline typename GeometricField<Type, PatchField, GeoMesh>::GeometricBoundaryField&
#else
    inline typename GeometricField<Type, PatchField, GeoMesh>::Boundary&
#endif
    boundaryFieldRef
    (
        GeometricField<Type, PatchField, GeoMesh>& fld,
        const bool updateAccessTime = true
    )
    {
#ifdef FOAMEXTEND
        return fld.boundaryField();
#else
        return fld.boundaryFieldRef(updateAccessTime);
#endif
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
