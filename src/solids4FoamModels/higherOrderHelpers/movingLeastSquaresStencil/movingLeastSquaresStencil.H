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
    movingLeastSquaresStencil

Description
    Constructs cell- and face-centred stencils for Moving Least Squares (MLS)
    reconstruction using a nearest-neighbour strategy.

    For local stencils, neighbouring cells are identified using layered
    topological searches on the fvMesh. For non-local (remote) contributions
    in parallel runs, candidate cells are identified using spatial bounding
    boxes and exchanged between processors via MPI. Remote stencil entries
    are stored using global indices and geometric information to ensure
    consistency across processor boundaries.

    Face-centred stencils are constructed first, and corresponding
    cell-centred stencils are subsequently assembled by combining the
    surrounding face stencils. Although constructing cell stencils directly
    would reduce the total number of stencil builds, the current approach
    was chosen for improved robustness, particularly for faces near
    processor boundaries and other geometrically challenging regions.

    The resulting stencils are intended for high-order MLS-based gradient
    and reconstruction operators in solids4foam, with explicit support for
    parallel execution.

SourceFile
    movingLeastSquaresStencil.C

Author
    Ivan Batistic, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef movingLeastSquaresStencil_H
#define movingLeastSquaresStencil_H

#include "Map.H"
#include "Field.H"
#include "tmp.H"
#include "globalIndex.H"
#include "treeBoundBox.H"
#include "CompactListList.H"
#include "boolList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;


class movingLeastSquaresStencil
{

    // Private data

        //- Store reference to the mesh
        const fvMesh& mesh_;

        //- Global cell addressing object for converting between local and
        //  global indices
        const globalIndex globalCells_;

        //- First halo depth is multiplied with this number to get estimation
        //  of multi-layer thickness
        const scalar haloDepthScale_;

        //- Stencil size for faces
        const label N_;

         //- Stencil size for cells
        const label Nc_;

        //- Marked cells at processor boundary
        mutable autoPtr<boolList> procCellsPtr_;

        //- Face-centred MLS stencils (global cell IDs)
        mutable autoPtr<CompactListList<label>> facesStencilPtr_;

        //- Cell-centred MLS stencils assembled from face stencils (global cell IDs)
        mutable autoPtr<CompactListList<label>> cellsStencilPtr_;

        //- Remote (off-processor) cells referenced by stencils, grouped per processor
        mutable autoPtr<List<labelList>> remoteCellsPerProcPtr_;

        //- Geometric centres corresponding to remoteCellsPerProc
        mutable autoPtr<List<vectorField>> remoteCentresPerProcPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingLeastSquaresStencil(const movingLeastSquaresStencil&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const movingLeastSquaresStencil&) = delete;

        //- Calculate axis-aligned bounding box for owned cells
        treeBoundBox calcOwnedCellsBox() const;

        //- Calculate axis-aligned bounding box for owned faces
        treeBoundBox calcOwnedFacesBox() const;

        //- Get average depth of first-halo cells for each processors
        List<scalar> calcFirstHaloDepth() const;

        //- Mark cells at processor boundary
        void calcProcessorCells() const;

        //- Check processor overlap between bounding box for remote processors
        /// cells and bounding box for local processor faces
        labelList checkProcessorOverlap
        (
            const List<treeBoundBox>& allOwnedCellsBox,
            const treeBoundBox& ownedFacesBox
        ) const;

        //- Get remote cells candidates for stencil construction
        List<labelList> remoteCandidates
        (
            const treeBoundBox& ownedFacesBox,
            const labelList& procToQuery
        ) const;

        //- Get cell centre for remote candidates
        List<vectorField> remoteCandidatesCellCentres
        (
            const List<labelList>& remoteCandidates,
            const labelList& procToQuery
        ) const;

        //- Build face stencil
        labelList buildFacesStencil
        (
             const label faceI,
             const List<labelList>& remoteCells,
             const List<vectorField>& remoteCellsCentres,
             const scalar relTol = 1e-8,
             const scalar overSampleFactor = 1.5,
             const scalar searchExpansionFactor = 1.3
        ) const;

        //- Filter unused cells and store cells IDs and their cell centres
        void filterUnusedCandidates
        (
            const List<labelList>& remoteCellsPerProc,
            const List<vectorField>& remoteCellsCentresPerProc
        ) const;

        //- Calculate faces stencil
        void calcFacesStencil() const;

        //- Calculate cells stencil
        void calcCellsStencil(const scalar relTol = 1e-8) const;

public:

    TypeName("movingLeastSquaresStencil");

    // Constructors
    movingLeastSquaresStencil
    (
        const fvMesh& mesh,
        const scalar haloDepthScale,
        const label N,
        const label Nc
    );

    // Desctructor
    virtual ~movingLeastSquaresStencil() = default;

    // Member Functions

        //- Marked cells at processor boundary
        const boolList& procCells() const;

        //- Returns face stencils
        const CompactListList<label>& facesStencil() const;

        //- Returns cell stencils
        const CompactListList<label>& cellsStencil() const;

        //- Return remote cells per processor written in global indexing
        const List<labelList>& remoteCellsPerProc() const;

        //- Return remote cell centres per processor, same order as remoteCells
        const List<vectorField>& remoteCentresPerProc() const;

        //- Clear all data, call after mesh topology change
        void clear() const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
