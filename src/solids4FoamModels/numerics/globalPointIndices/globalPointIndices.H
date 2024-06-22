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
    globalPointIndices

Description
    Class for creating a unique global indices for the points as well as a list
    of which points are owned by each processor. The heavy lifting is performed
    by the globalPoints and globalMeshData classes.

Author
    Philip Cardiff, UCD.

SourceFiles
    globalPointIndices.C

\*---------------------------------------------------------------------------*/

#ifndef globalPointIndices_H
#define globalPointIndices_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef FixedList<label, 2> procPoint;
    typedef List<procPoint> procPointList;

/*---------------------------------------------------------------------------*\
                       Class globalPointIndices Declaration
\*---------------------------------------------------------------------------*/

class globalPointIndices
{
    // Private data

        //- Const reference to the polyMesh
        const polyMesh& mesh_;

        //- List indicating which points are owned by the current processor
        boolList ownedByThisProc_;

        //- Local-to-global point map
        //  For each point on the current processor, return the global point
        //  index
        //  -1 for points not owned by this proc
        //  For larger meshes, label should be long int or long long int
        labelList localToGlobalPointMap_;

        //- Stencil size, which only counts points owned by this proc
        labelList stencilSizeOwned_;

        //- Stencil size, which only counts points not-owned by this proc. In
        //  general, some of these points will be on this processor and some
        //  will not
        labelList stencilSizeNotOwned_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        globalPointIndices(const globalPointIndices&);

        //- Disallow default bitwise assignment
        void operator=(const globalPointIndices&);

public:

    //- Runtime type information
    TypeName("globalPointIndices");


    // Constructors

        //- Construct given the a mesh
        globalPointIndices(const polyMesh& mesh);


    // Destructor

        virtual ~globalPointIndices()
        {}


    // Member Functions

        // Access

            //- Const access to ownedByThisProc list
            const boolList& ownedByThisProc() const
            {
                return ownedByThisProc_;
            }

            //- Const access to local-to-global point map
            const labelList& localToGlobalPointMap() const
            {
                return localToGlobalPointMap_;
            }

            //- Const access stencilSizeOwned
            const labelList& stencilSizeOwned() const
            {
                return stencilSizeOwned_;
            }

            //- Const access stencilSizeNotOwned
            const labelList& stencilSizeNotOwned() const
            {
                return stencilSizeNotOwned_;
            }

            //- Return sum of stencilSizeOwned and stencilSizeNotOwned
            labelList stencilSize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
