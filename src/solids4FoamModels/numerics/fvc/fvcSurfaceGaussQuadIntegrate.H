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

InNamespace
    Foam::fvc

Description
    Surface integrate surfaceField creating a volField.
    Surface sum a surfaceField creating a volField.

SourceFiles
    fvcSurfaceIntegrate.C

Author
    Ivan Batistic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/


#ifndef fvcSurfaceGaussQuadIntegrate_H
#define fvcSurfaceGaussQuadIntegrate_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    void surfaceGaussQuadIntegrate
    (
        Field<Type>&,
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>>
    surfaceGaussQuadIntegrate
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>>
    surfaceGaussQuadIntegrate
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcSurfaceGaussQuadIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
