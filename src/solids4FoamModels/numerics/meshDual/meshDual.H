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
    meshDual

Description
    Class for constructing a dual (fvMesh) mesh of the provided primary
    (fvMesh) mesh.

    Maps between the primary and dual mesh are constructed and stored.

Author
    Philip Cardiff, UCD.

SourceFiles
    meshDual.C

\*---------------------------------------------------------------------------*/

#ifndef meshDual_H
#define meshDual_H

#include "fvMesh.H"
#include "dualMeshToMeshMap.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class meshDual Declaration
\*---------------------------------------------------------------------------*/

class meshDual
:
    public fvMesh
{
    // Private data

        //- Const-access to the primary mesh
        const fvMesh& mesh_;

        //- Const-access to the dictionary
        const dictionary& dict_;

        //- Dual mesh to primary mesh map
        mutable autoPtr<dualMeshToMeshMap> dualMeshToMeshMapPtr_;


    // Private Member Functions

        //- Make the dual mesh
        //  This function assumes the primary mesh has already been copied to
        //  the dual mesh
        void makeDualMesh();

        //- Mark features for dual mesh
        //  Copied from polyDualMeshApp.C
        void simpleMarkFeatures
        (
            const polyMesh& mesh,
            const PackedBoolList& isBoundaryEdge,
            const scalar featureAngle,
            const bool concaveMultiCells,
            const bool doNotPreserveFaceZones,
            labelList& featureFaces,
            labelList& featureEdges,
            labelList& singleCellFeaturePoints,
            labelList& multiCellFeaturePoints
        ) const;

        //- Disallow default bitwise copy construct
        meshDual(const meshDual&);

        //- Disallow default bitwise assignment
        void operator=(const meshDual&);

public:

    //- Runtime type information
    TypeName("meshDual");


    // Constructors

        //- Construct from the primary mesh and a dictionary
        meshDual(const fvMesh& mesh, const dictionary& dict);


    // Destructor

        virtual ~meshDual()
        {}


    // Member Functions

        // Access

            //- Const-access to the primary mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Const-access to the time
            const Time& runTime() const
            {
                return mesh_.time();
            }

            //- Const-access to the dual mesh to mesh map
            const dualMeshToMeshMap& dualMeshMap() const
            {
                if (dualMeshToMeshMapPtr_.empty())
                {
                    FatalErrorIn
                    (
                        "const dualMeshToMeshMap& dualMeshToMeshMap() const"
                    )   << "Pointer not set!" << abort(FatalError);
                }

                return dualMeshToMeshMapPtr_();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
