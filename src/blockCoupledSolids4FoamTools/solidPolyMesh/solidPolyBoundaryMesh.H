/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    solidPolyBoundaryMesh

Description
    polyBoundaryMesh for the solidPolyMesh class.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    solidPolyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef solidPolyBoundaryMesh_H
#define solidPolyBoundaryMesh_H

#include "fvPatchList.H"

//#include "faceListFwd.H"
#include "lduInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class solidPolyMesh;
class polyBoundaryMesh;
//class globalTetPolyPatch;

/*---------------------------------------------------------------------------*\
                 Class solidPolyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class solidPolyBoundaryMesh
// :
//     public fvPatchList
{
    // Private data

        //- Const reference to mesh
        const solidPolyMesh& mesh_;

        //- Const reference to polyBoundaryMesh
        const polyBoundaryMesh& polyBdry_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidPolyBoundaryMesh(const solidPolyBoundaryMesh&);

        //- Disallow default bitwise assignment
        void operator=(const solidPolyBoundaryMesh&);

        //- Add fvPatches corresponding to the given solidPolyBoundaryMesh
        void addPatches(const polyBoundaryMesh&);

public:

    // Constructors

        //- Construct from polyBoundaryMesh
        solidPolyBoundaryMesh
        (
            const solidPolyMesh&,
            const polyBoundaryMesh&
        );


    // Member functions

        // Access

        //- Return the mesh reference
            const solidPolyMesh& mesh() const
            {
                return mesh_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;

            //- Return const reference to polyBoundaryMesh
            const polyBoundaryMesh& boundaryMesh() const
            {
                return polyBdry_;
            }

            //- Return reference to shared processor point patch
            //const globalTetPolyPatch& globalPatch() const;

        // Edit

            // Update topology
            // void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
