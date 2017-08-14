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
    faMeshDecomposition

Description
    Automatic faMesh decomposition class

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    faMeshDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef faMeshDecomposition_H
#define faMeshDecomposition_H

#include "fvMesh.H"
#include "faMesh.H"
#include "labelList.H"
#include "SLList.H"
#include "PtrList.H"
#include "point.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

/*---------------------------------------------------------------------------*\
                           Class faMeshDecomposition Declaration
\*---------------------------------------------------------------------------*/

class faMeshDecomposition
:
    public faMesh
{
    // Private data

        //- Mesh decomposition control dictionary
        IOdictionary decompositionDict_;

        //- Number of processors in decomposition
        label nProcs_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;

        //- Processor label for each cell
        labelList faceToProc_;

        //- Face labels for each processor mesh
        labelListList procFaceLabels_;

        //-
        List<Map<label> > procMeshEdgesMap_;

        //- Number of internal edges for each processor mesh
        labelList procNInternalEdges_;

        //- Edge labels for patches of processor meshes
        List<List<List<label> > > procPatchEdgeLabels_;

        //- Labels of points for each processor
        labelListList procPatchPointAddressing_;

        //- Labels of edges for each processor
        labelListList procPatchEdgeAddressing_;

        //- Labels of edges for each processor
        labelListList procEdgeAddressing_;

        //- Labels of faces for each processor
        labelListList procFaceAddressing_;

        //- Original patch index for every processor patch
        labelListList procBoundaryAddressing_;

        //- Sizes for processor mesh patches
        // Excludes inter-processor boundaries
        labelListList procPatchSize_;

        //- Start indices for processor patches
        // Excludes inter-processor boundaries
        labelListList procPatchStartIndex_;

        //- Neighbour processor ID for inter-processor boundaries
        labelListList procNeighbourProcessors_;

        //- Sizes for inter-processor patches
        labelListList procProcessorPatchSize_;

        //- Start indices for inter-processor patches
        labelListList procProcessorPatchStartIndex_;

        //- List of globally shared point labels
        labelList globallySharedPoints_;

        //- Are there cyclic-parallel faces
        bool cyclicParallel_;


    // Private Member Functions

        void distributeFaces();

public:

    // Constructors

        //- Construct from fvMesh
        faMeshDecomposition(const fvMesh& mesh);


    // Destructor

        ~faMeshDecomposition();


    // Member Functions

        //- Number of processor in decomposition
        label nProcs() const
        {
            return nProcs_;
        }

        //- Is the decomposition data to be distributed for each processor
        bool distributed() const
        {
            return distributed_;
        }

        //- Decompose mesh. Optionally remove zero-sized patches.
        void decomposeMesh(const bool filterEmptyPatches);

        //- Write decomposition
        bool writeDecomposition();

        //- Cell-processor decomposition labels
        const labelList& faceToProc() const
        {
            return faceToProc_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
