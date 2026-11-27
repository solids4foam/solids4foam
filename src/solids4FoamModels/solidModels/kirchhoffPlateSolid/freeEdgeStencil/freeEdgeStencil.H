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

Namespace
    Foam::freeEdgeStencil

Description
    Connectivity of the edges along a finite area patch, used by the free
    edge conditions of the Kirchhoff plate model to differentiate along the
    plate edge and to classify its corners.

Author
    Ivan Batistic, FSB.

SourceFiles
    freeEdgeStencil.C

\*---------------------------------------------------------------------------*/

#ifndef freeEdgeStencil_H
#define freeEdgeStencil_H

#include "faPatch.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace freeEdgeStencil Declaration
\*---------------------------------------------------------------------------*/

namespace freeEdgeStencil
{
    //- For each edge of the patch, the other edges of the patch sharing one
    //  of its points: two for an edge inside the patch, one for an edge at
    //  an end of the patch
    labelListList edgeNeighbours(const faPatch& p);

    //- For each edge of the patch, its points that are not shared with
    //  another edge of the patch, i.e. the corners of the plate on the edge
    labelListList edgeEndPoints(const faPatch& p);

    //- True if the point lies on another patch where the field fixes its
    //  value, i.e. the corner joins a clamped or supported edge
    bool fixesValueAtPoint
    (
        const areaScalarField& fld,
        const faPatch& p,
        const label pointI
    );

    //- Unit tangent of each edge of the patch, oriented consistently along
    //  each chain of connected edges so that the edge values can be
    //  differentiated along the patch
    vectorField edgeTangents(const faPatch& p, const labelListList& nbrs);

} // End namespace freeEdgeStencil

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
