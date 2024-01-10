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
    patchCorrectionVectors

Description
    Function that returns the correction vectors on a patch given the patch
    object.

    The reason for this function is because the delta patch function returns the
    true delta vector in foam-extend, but it returns the normal vector in
    OpenFOAM.com and OpenFOAM.org.

SourceFile
    patchCorrectionVectors.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef patchCorrectionVectors_H
#define patchCorrectionVectors_H

#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Return the patch correction vectors
    tmp<vectorField> patchCorrectionVectors(const fvPatch& patch);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
