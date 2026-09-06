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

Description
    Which faces lie between two different materials.

    A least-squares gradient is computed from a cell and its neighbours, and a
    stencil that reaches across a material interface mixes two materials into
    one gradient. These functions say where those interfaces are so that a
    stencil can stop at them.

    A material is not the same thing as a cell zone. The materials are the
    entries of the "mechanical" list in mechanicalProperties, and each is named
    after the cell zone it occupies. A mesh may carry other zones that name no
    material - plateHole has zones "aluminium" and "steel" and a single
    material - and a face between two such zones is not a material interface.
    Keying on zones alone would cut a stencil in half for no reason and quietly
    degrade the gradient, so the material names decide.

SourceFiles
    cellZoneInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cellZoneInterface_H
#define cellZoneInterface_H

#include "boolList.H"
#include "compatibilityFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //- The material each cell belongs to, or -1 where a cell belongs to none.
    //  Every cell is -1 when the case has a single material, which is the
    //  common case and means there are no interfaces anywhere
    labelList cellMaterialID(const fvMesh& mesh);

    //- For each internal face, whether it lies between two different materials
    tmp<Field<bool>> cellZoneInterface
    (
        const fvMesh& mesh, const bool debug = true
    );

    //- For each face of a coupled patch, whether the cell on the other side
    //  is of a different material.
    //
    //  A material interface that runs along a processor boundary is not an
    //  internal face on either side, so without this a stencil would reach
    //  across it in parallel and not in serial, and the two would disagree
    List<boolList> cellZoneInterfaceCoupled(const fvMesh& mesh);

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
