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
    Foam::quadratureValuePatchField

Description
    Interface for fixed-value patch fields that evaluate their value at the
    face quadrature points used by the high-order least-squares schemes.

    A patch field inherits this interface alongside its usual base class.
    The schemes find it with a dynamic_cast (see patchQuadratureValues.H), so
    they do not need to know the patch field types.

SourceFiles
    patchQuadratureValues.H

Author
    Amirhossein Taran, UCD.

\*---------------------------------------------------------------------------*/

#ifndef quadratureValuePatchField_H
#define quadratureValuePatchField_H

#include "autoPtr.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshQuadrature;

/*---------------------------------------------------------------------------*\
                  Class quadratureValuePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class quadratureValuePatchField
{
public:

    //- Destructor
    virtual ~quadratureValuePatchField() = default;


    // Member Functions

        //- Return the patch values at the face quadrature points of quad
        virtual autoPtr<CompactListList<Type>> quadratureValues
        (
            const fvMeshQuadrature& quad
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
