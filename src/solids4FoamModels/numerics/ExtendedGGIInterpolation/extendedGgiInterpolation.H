/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    extededGGIInterpolation

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Description
    Mass-conservative face interpolation: typedefs for polyPatch and faceZone

\*---------------------------------------------------------------------------*/

#ifndef extendedGGIInterpolation_H
#define extendedGGIInterpolation_H

#include "ExtendedGGIInterpolationTemplate.H"
#include "PrimitivePatchTemplate.H"
#include "face.H"
#include "SubList.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ExtendedGGIInterpolation
    <
        PrimitivePatch<face, SubList, const pointField&>,
        PrimitivePatch<face, SubList, const pointField&>
    >   extendedGgiInterpolation;

    typedef ExtendedGGIInterpolation
    <
        PrimitivePatch<face, List, const pointField&>,
        PrimitivePatch<face, List, const pointField&>
    >   extendedGgiZoneInterpolation;

    typedef ExtendedGGIInterpolation
    <
        PrimitivePatch<face, List, pointField>,
        PrimitivePatch<face, List, pointField>
    >   extendedGgiStandAlonePatchInterpolation;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
