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
    newGgiInterpolation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Description
    Mass-conservative face interpolation: typedefs for polyPatch and faceZone

\*---------------------------------------------------------------------------*/

#ifndef newGgiInterpolation_H
#define newGgiInterpolation_H

#include "newGGIInterpolationTemplate.H"
#include "PrimitivePatchTemplate.H"
#include "face.H"
#include "SubList.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef newGGIInterpolation
    <
        PrimitivePatch<face, SubList, const pointField&>,
        PrimitivePatch<face, SubList, const pointField&>
    >   newGgiInterpolation;

    typedef newGGIInterpolation
    <
        PrimitivePatch<face, List, const pointField&>,
        PrimitivePatch<face, List, const pointField&>
    >   newGgiZoneInterpolation;

    typedef newGGIInterpolation
    <
        PrimitivePatch<face, List, pointField>,
        PrimitivePatch<face, List, pointField>
    >   newGgiStandAlonePatchInterpolation;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
