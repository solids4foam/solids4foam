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
    dualMeshToMeshMap

Description
    Class to provide maps between a dual mesh and the corresponding primary
    mesh, for example, dualFaceToCellMap, etc.

    The class is constructed from a meshDualiser object, which already
    contains many of the required maps. Additional maps are constructed using
    these base maps.

    All maps are created on construction.

    Take care as entries will be -1 if the map does not exist for the entity.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    dualMeshToMeshMap.C

\*---------------------------------------------------------------------------*/

#ifndef dualMeshToMeshMap_H
#define dualMeshToMeshMap_H

#include "meshDualiser.H"
#include "point.H"
//#include "tetPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                     Class dualMeshToMeshMap Declaration
\*---------------------------------------------------------------------------*/

class dualMeshToMeshMap
{
    // Private data

        //- Const reference to primary mesh
        const polyMesh& mesh_;

        //- Const reference to dual mesh
        const polyMesh& dualMesh_;

        //- From point to dual cell. Either single entry or one entry per
        //  pointCells
        const labelListList pointToDualCells_;

        //- From point to dual cell, assuming a single entry
        labelList pointToDualCell_;

        //- From point to dual faces
        //  Points may map to multiple dual faces, particularly at boundaries
        //  Note that this map is approximate and returns all faces in the dual
        // cell that the point maps to
        labelListList pointToDualFaces_;

        //- From point to dual point (or -1 if not feature point)
        const labelList pointToDualPoint_;

        //- From cell to dual point (at cell centre). All cells become points
        const labelList cellToDualPoint_;

        //- From face to dual point
        const labelList faceToDualPoint_;

        //- From edge to dual point
        const labelList edgeToDualPoint_;

        //- From dual face to cell
        //  This may not be uniquely defined for some boundary faces; in that
        //  case, the value is set to -1
        labelList dualFaceToCell_;

        //- From dual cell to point
        labelList dualCellToPoint_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dualMeshToMeshMap(const dualMeshToMeshMap&);

        //- Disallow default bitwise assignment
        void operator=(const dualMeshToMeshMap&);

public:

    //- Runtime type information
    TypeName("dualMeshToMeshMap");


    // Constructors

        //- Construct from primary mesh, dual mesh and meshDualiser
        dualMeshToMeshMap
        (
            const polyMesh& mesh,
            const polyMesh& dualMesh,
            const meshDualiser& dualiser
        );

    // Destructor

        virtual ~dualMeshToMeshMap();


    // Member Functions

        // Access

            //- Const reference to primary mesh
            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Const reference to dual mesh
            const polyMesh& dualMesh() const
            {
                return dualMesh_;
            }

            //- Return point to dual cell
            const labelListList& pointToDualCells() const
            {
                return pointToDualCells_;
            }

            //- Return point to dual cell
            const labelList& pointToDualCell() const
            {
                return pointToDualCell_;
            }

            //- From point to dual face(s)
            const labelListList& pointToDualFaces() const
            {
                return pointToDualFaces_;
            }

            //- From point to dual point
            const labelList& pointToDualPoint() const
            {
                return pointToDualPoint_;
            }

            //- From cell to dual point
            const labelList& cellToDualPoint() const
            {
                return cellToDualPoint_;
            }

            //- From face to dual point
            const labelList& faceToDualPoint() const
            {
                return faceToDualPoint_;
            }

            //- From edge to dual point
            const labelList& edgeToDualPoint() const
            {
                return edgeToDualPoint_;
            }

            //- From dual face to primary cell
            const labelList& dualFaceToCell() const
            {
                return dualFaceToCell_;
            }

            //- From dual cell to primary point
            const labelList& dualCellToPoint() const
            {
                return dualCellToPoint_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
