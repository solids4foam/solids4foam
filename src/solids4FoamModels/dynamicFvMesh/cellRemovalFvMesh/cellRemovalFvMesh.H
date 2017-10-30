/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cellRemovalFvMesh

Description
    Mesh allowing cell removal based on a run-time selectable cell removal law.

SourceFiles
    cellRemovalFvMesh.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cellRemovalFvMesh_H
#define cellRemovalFvMesh_H

#include "dynamicFvMesh.H"
#include "directTopoChange.H"
#include "removeCells.H"
#include "cellRemovalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellRemovalFvMesh Declaration
\*---------------------------------------------------------------------------*/

class cellRemovalFvMesh
:
    public dynamicFvMesh
{
    // Private Data Members

        //- Settings dictionary
        dictionary dict_;

        //- Write mesh when there is a topo change
        const Switch writeMesh_;

        //- Remove dead cells: these are small islands of cells which have
        //  become disconnected from the main mesh. We define an island as a
        //  bunch of cells disconnected from the main mesh having less than 1%
        //  of the total number of cells in the mesh
        const Switch removeDeadCells_;

        //- Topo changer for changing the mesh
        mutable directTopoChange* meshModPtr_;

        //- Cell remover
        removeCells cellRemover_;

        //- Law for governing which cells are removed
        autoPtr<cellRemovalLaw> lawPtr_;


    // Private Member Functions

        //- Calculate mesh modifier
        void calcMeshMod() const;

        //- Return the mesh modifier
        directTopoChange& meshModifier();

        //- Update volField values on the newly exposed faces
        //  The values on the new faces are extrapolated from the internal
        //  field; currently we use zero gradient extrapolation
        template<class Type>
        void updateVolFieldsExposedFaces
        (
            const mapPolyMesh& map,
            const labelList& exposedFaces
        ) const;

        // Delete demand driven data
        void clearOut() const;

        //- Disallow default bitwise copy construct
        cellRemovalFvMesh(const cellRemovalFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const cellRemovalFvMesh&);


public:

    //- Runtime type information
    TypeName("cellRemovalFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit cellRemovalFvMesh(const IOobject& io);


    // Destructor

        virtual ~cellRemovalFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cellRemovalFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
