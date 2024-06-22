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
    Foam::advancingFrontAMIS4F

Description
    Base class for Arbitrary Mesh Interface (AMI) methods

SourceFiles
    advancingFrontAMIS4F.C

\*---------------------------------------------------------------------------*/

#ifndef advancingFrontAMIS4F_H
#define advancingFrontAMIS4F_H

#include "className.H"
#include "DynamicList.H"
#include "faceAreaIntersect.H"
#include "pointList.H"
#include "AMIInterpolationS4F.H"
#include "standAlonePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class advancingFrontAMIS4F Declaration
\*---------------------------------------------------------------------------*/

class advancingFrontAMIS4F
:
    public AMIInterpolationS4F
{

private:

    // Private Member Functions

        //- No copy assignment
        void operator=(const advancingFrontAMIS4F&) = delete;


        // Parallel operations

            label calcOverlappingProcs
            (
                const List<treeBoundBoxList>& procBb,
                const treeBoundBox& bb,
                boolList& overlaps
            ) const;

            void distributePatches
            (
                const mapDistribute& map,
                const standAlonePatch& pp,
                const globalIndex& gi,
                List<faceList>& faces,
                List<pointField>& points,
                List<labelList>& tgtFaceIDs
            ) const;

            void distributeAndMergePatches
            (
                const mapDistribute& map,
                const standAlonePatch& tgtPatch,
                const globalIndex& gi,
                faceList& tgtFaces,
                pointField& tgtPoints,
                labelList& tgtFaceIDs
            ) const;

            autoPtr<mapDistribute> calcProcMap
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch
            ) const;


protected:

    // Protected data

        //- Storage for src-side triangle decomposition
        List<DynamicList<face>> srcTris_;

        //- Storage for tgt-side triangle decomposition
        List<DynamicList<face>> tgtTris_;

        //- Demand-driven extended target mesh (distributed parallel usage)
        autoPtr<standAlonePatch> extendedTgtPatchPtr_;

        //- Extended patch faces
        faceList extendedTgtFaces_;

        //- Extended patch points
        pointField extendedTgtPoints_;

        //- Extended patch face IDs
        labelList extendedTgtFaceIDs_;

        //- Extended patch map
        autoPtr<mapDistribute> extendedTgtMapPtr_;

        //- Labels of faces that are not overlapped by any target faces
        //- (should be empty for correct functioning for fully covered AMIs)
        labelList srcNonOverlap_;

        //- Octree used to find face seeds
        autoPtr<indexedOctree<treeType>> treePtr_;

        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;


    // Protected Member Functions

        // Helper functions

            //- Create a map that extends tgtPatch so that it covers srcPatch
            void createExtendedTgtPatch();

            //- Check AMI patch coupling
            void checkPatches() const;

            virtual bool calculate
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch,
                const autoPtr<searchableSurface>& surfPtr = nullptr
            );

            //- Initialise walk and return true if all ok
            bool initialiseWalk
            (
                label& srcFacei,
                label& tgtFacei
            );

            //- Write triangle intersection to OBJ file
            void writeIntersectionOBJ
            (
                const scalar area,
                const face& f1,
                const face& f2,
                const pointField& f1Points,
                const pointField& f2Points
            ) const;


        // Common AMI method functions

            label findTargetFace
            (
                const label srcFacei,
                const UList<label>& excludeFaces = UList<label>::null(),
                const label srcFacePti = -1
            ) const;

            //- Add faces neighbouring facei to the ID list
            void appendNbrFaces
            (
                const label facei,
                const standAlonePatch& patch,
                const DynamicList<label>& visitedFaces,
                DynamicList<label>& faceIDs
            ) const;

            //- Helper function to decompose a patch
            void triangulatePatch
            (
                const standAlonePatch& patch,
                List<DynamicList<face>>& tris,
                List<scalar>& magSf
            ) const;


public:

    //- Runtime type information
    TypeName("advancingFrontAMIS4F");

    // Constructors

        //- Construct from components
        advancingFrontAMIS4F
        (
            const dictionary& dict,
            const bool reverseTarget
        );

        //- Construct from components
        advancingFrontAMIS4F
        (
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const faceAreaIntersect::triangulationMode triMode =
                faceAreaIntersect::tmMesh
        );

        //- Construct as copy
        advancingFrontAMIS4F(const advancingFrontAMIS4F& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolationS4F> clone() const
        {
            return autoPtr<AMIInterpolationS4F>(new advancingFrontAMIS4F(*this));
        }


    //- Destructor
    virtual ~advancingFrontAMIS4F() = default;


    // Member Functions

        //- Return const access to the source patch
        inline const standAlonePatch& srcPatch() const;

        //- Return const access to the target patch
        inline const standAlonePatch& tgtPatch() const;

        //- Labels of faces that are not overlapped by any target faces
        //  Note: this should be empty for correct functioning
        inline const labelList& srcNonOverlap() const;

        //- Flag to indicate that interpolation patches are conformal, i.e.
        //- should fully cover each other
        virtual bool conformal() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "advancingFrontAMIS4FI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
