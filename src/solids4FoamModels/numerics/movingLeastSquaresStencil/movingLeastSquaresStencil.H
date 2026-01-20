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
    Build cell and face centered stencils using neirest N neighbbours
    approach. Remote cells are stored using global coordinates. Local cells
    are collected using layer technique while linear search is conducted for
    remote cells.

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
#include "fvMesh.H"
#include "globalIndex.H"
#include "Pstream.H"
#include "PstreamBuffers.H"
#include "treeBoundBox.H"
#include "processorPolyPatch.H"
#include "surfaceFields.H"
#include "volFields.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

        //- Stencil size
        const label N_;

        //- Marked cells at processor boundary
        boolList procPatchesCells_;

        //- Face stencils (global IDs)
        CompactListList<label> facesStencil_;

         //- Cell centre stencils (global IDs)
        CompactListList<label> cellsStencil_;

        //- Remote cells used by stencils (global IDs), grouped by proc
        List<labelList> remoteCellsPerProc_;

        //- Matching cell centres for remoteCellsPerProc
        List<vectorField> remoteCentresPerProc_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingLeastSquaresStencil(const movingLeastSquaresStencil&);

        //- Disallow default bitwise assignment
        void operator=(const movingLeastSquaresStencil&);

        //- Calculate axis-aligned bounding box for owned cells
        treeBoundBox calcOwnedCellsBox() const;

        //- Calculate axis-aligned bounding box for owned faces
        treeBoundBox calcOwnedFacesBox() const;

        //- Get average depth of first-halo cells for each processors
        List<scalar> calcFirstHaloDepth() const;

        //- Check processor overlap between bounding box for remote processors
        /// cells and bounding box for local processor faces
        labelList checkProcessorOverlap
        (
            const List<treeBoundBox>& allOwnedCellsBox,
            const treeBoundBox& ownedFacesBox
        );

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
             const scalar relTol = 1e-10,
             const scalar overSampleFactor = 1.5
        ) const;

        //- Filter unused cells and store cells IDs and their cell centres
        void filterUnusedCandidates
        (
            const List<labelList>& remoteCellsPerProc,
            const List<vectorField>& remoteCellsCentresPerProc
         );

        //- Calculate faces stencil
        void calcFacesStencil();


public:

    TypeName("movingLeastSquaresStencil");

    // Constructors
    movingLeastSquaresStencil
    (
        const fvMesh& mesh,
        const scalar& haloDepthScale,
        const label& N
    );

    // Desctructor
    virtual ~movingLeastSquaresStencil() = default;

    // Member Functions

        //- Returns face stencils
        inline const CompactListList<label>& facesStencil() const
        {
            return facesStencil_;
        }

        //- Return stencil for faceI
        inline const labelUList faceStencil(const label faceI) const
        {
            return facesStencil_[faceI];
        }


};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
