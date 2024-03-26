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

Function
    pointFieldFunctions

Description
    Functions for point fields

    The reason for this function is because the delta patch function returns the
    true delta vector in foam-extend, but it returns the normal vector in
    OpenFOAM.com and OpenFOAM.org.

SourceFile
    pointFieldFunctions.C

Author
    Federico Mazzanti, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef pointFieldFunctions_H
#define pointFieldFunctions_H

#include "fvPatch.H"
#include "pointFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Return the symmetric components of a pointTensorField
    tmp<pointSymmTensorField> symm(const pointTensorField& ptr);

    // Return the deviatoric components of a pointSymmTensorField
    tmp<pointSymmTensorField> dev(const pointSymmTensorField& ptr);

    // Return the trace of a pointSymmTensorField
    tmp<pointScalarField> tr(const pointSymmTensorField& ptr);
    
    // Return the cofactor of a pointTensorField
    tmp<pointTensorField> cof(const pointTensorField& ptr);
    
    // Return the square of a pointScalarField
    tmp<pointScalarField> sqr(const pointScalarField& ptr);
    
    // Return the power of a pointScalarField
    tmp<pointScalarField> pow(const pointScalarField& ptr, const scalar& power);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
