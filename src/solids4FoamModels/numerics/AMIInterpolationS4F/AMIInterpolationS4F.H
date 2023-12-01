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
    Foam::AMIInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitive patches with an arbitrary mesh interface (AMI).

    Based on the algorithm given in:

        Conservative interpolation between volume meshes by local Galerkin
        projection, Farrell PE and Maddison JR, 2011, Comput. Methods Appl.
        Mech Engrg, Volume 200, Issues 1-4, pp 89-100

    Interpolation requires that the two patches should have opposite
    orientations (opposite normals).  The 'reverseTarget' flag can be used to
    reverse the orientation of the target patch.

\*---------------------------------------------------------------------------*/

#ifndef AMIInterpolationS4F_H
#define AMIInterpolationS4F_H

#include "className.H"
#include "DynamicList.H"
#include "searchableSurface.H"
#include "treeBoundBoxList.H"
#include "boolList.H"
#include "primitivePatch.H"
#include "faceAreaIntersect.H"
#include "globalIndex.H"
#include "ops.H"
#include "refPtr.H"
#include "Enum.H"
#include "pointList.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"
#include "standAlonePatch.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class AMIInterpolationS4F Declaration
\*---------------------------------------------------------------------------*/

class AMIInterpolationS4F
{
public:

    // Public data types

        static bool cacheIntersections_;


protected:

    //- Local typedef to octree tree-type
    typedef treeDataPrimitivePatch<standAlonePatch> treeType;

    // Protected data

        //- Flag to indicate that the two patches must be matched/an overlap
        //- exists between them
        bool requireMatch_;

        //- Flag to indicate that the two patches are co-directional and
        //- that the orientation of the target patch should be reversed
        const bool reverseTarget_;

        //- Threshold weight below which interpolation is deactivated
        scalar lowWeightCorrection_;

        //- Index of processor that holds all of both sides. -1 in all other
        //- cases
        label singlePatchProc_;

            //***********************************************************//
            //  JN: from extendedggiinterpolation - possibly add here?
            //  JN: don't forget to add in .C file
            //***********************************************************//

            //- Master face into which each point of slave patch is projected
            //mutable List<labelPair>* slavePointAddressingPtr_;

            //- Distance to intersection for master patch points
            //mutable scalarField* slavePointDistancePtr_;

            //***********************************************************//

        // Source patch

            //- Source face areas
            scalarList srcMagSf_;

            //- Addresses of target faces per source face
            labelListList srcAddress_;

            //- Weights of target faces per source face
            scalarListList srcWeights_;

            //- Sum of weights of target faces per source face
            scalarField srcWeightsSum_;

            //- Centroid of target faces per source face
            pointListList srcCentroids_;

            //- Source patch points if input points are manipulated, e.g.
            //- projected
            pointField srcPatchPts_;

            //- Source patch using manipulated input points
            //refPtr<primitivePatch> tsrcPatch0_;
            refPtr<standAlonePatch> tsrcPatch0_;

            //- Source map pointer - parallel running only
            autoPtr<mapDistribute> srcMapPtr_;


        // Target patch

            //- Target face areas
            scalarList tgtMagSf_;

            //- Addresses of source faces per target face
            labelListList tgtAddress_;

            //- Weights of source faces per target face
            scalarListList tgtWeights_;

            //- Sum of weights of source faces per target face
            scalarField tgtWeightsSum_;

            //- Centroid of source faces per target face
            pointListList tgtCentroids_;

            //- Target patch points if input points are manipulated, e.g.
            //- projected
            pointField tgtPatchPts_;

            //- Target patch using manipulated input points
            //refPtr<primitivePatch> ttgtPatch0_;
            refPtr<standAlonePatch> ttgtPatch0_;

            //- Target map pointer - parallel running only
            autoPtr<mapDistribute> tgtMapPtr_;

            //- Up-to-date flag
            bool upToDate_;

            //- Patches are globalPolyPatches (same copy of global patch on all
            //  processors
            bool useGlobalPolyPatch_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const AMIInterpolationS4F&) = delete;


        // Initialisation

            //- Reset the octree for the patch face search
            autoPtr<indexedOctree<treeType>> createTree
            (
                const standAlonePatch& patch
            ) const;

            //- Calculate if patches are on multiple processors
            label calcDistribution
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch
            ) const;

            //- Project points to surface
            void projectPointsToSurface
            (
                const searchableSurface& surf,
                pointField& pts
            ) const;


        // Access

            //- Return the orginal src patch with optionally updated points
            inline const standAlonePatch& srcPatch0() const;

            //- Return the orginal tgt patch with optionally updated points
            inline const standAlonePatch& tgtPatch0() const;


        // Evaluation

            //- Normalise the (area) weights - suppresses numerical error in
            //- weights calculation
            //  NOTE: if area weights are incorrect by 'a significant amount'
            //     normalisation may stabilise the solution, but will introduce
            //     numerical error!
            static void normaliseWeights
            (
                const scalarList& patchAreas,
                const word& patchName,
                const labelListList& addr,
                scalarListList& wght,
                scalarField& wghtSum,
                const bool conformal,
                const bool output,
                const scalar lowWeightTol
            );


        // Constructor helpers

            static void agglomerate
            (
                const autoPtr<mapDistribute>& targetMap,
                const scalarList& fineSrcMagSf,
                const labelListList& fineSrcAddress,
                const scalarListList& fineSrcWeights,

                const labelList& sourceRestrictAddressing,
                const labelList& targetRestrictAddressing,

                scalarList& srcMagSf,
                labelListList& srcAddress,
                scalarListList& srcWeights,
                scalarField& srcWeightsSum,
                autoPtr<mapDistribute>& tgtMap
            );


public:

    //- Runtime type information
    TypeName("AMIInterpolationS4F");

    // Selection tables

        //- Selection table for dictionary construction
        declareRunTimeSelectionTable
        (
            autoPtr,
            AMIInterpolationS4F,
            dict,
            (
                const dictionary& dict,
                const bool reverseTarget
            ),
            (
                dict,
                reverseTarget
            )
        );

        //- Selection table for component-wise construction
        declareRunTimeSelectionTable
        (
            autoPtr,
            AMIInterpolationS4F,
            component,
            (
                const bool requireMatch,
                const bool reverseTarget,
                const scalar lowWeightCorrection
            ),
            (
                requireMatch,
                reverseTarget,
                lowWeightCorrection
            )
        );

        //- Selector for dictionary
        static autoPtr<AMIInterpolationS4F> New
        (
            const word& modelName,
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Selector for components
        static autoPtr<AMIInterpolationS4F> New
        (
            const word& modelName,
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1
        );


    // Constructors

        //- Construct from dictionary
        AMIInterpolationS4F
        (
            const dictionary& dict,
            const bool reverseTarget = false,
            const bool useGlobalPolyPatch = false
        );

        //- Construct from components
        AMIInterpolationS4F
        (
            const standAlonePatch& srcPatch,
            const standAlonePatch& tgtPatch,
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const bool useGlobalPolyPatch = false
        );

        //- Construct from components with patches; added by JN
        AMIInterpolationS4F
        (
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const bool useGlobalPolyPatch = false
        );

        //- Construct from agglomeration of AMIInterpolationS4F. Agglomeration
        //- passed in as new coarse size and addressing from fine from coarse
        AMIInterpolationS4F
        (
            const AMIInterpolationS4F& fineAMI,
            const labelList& sourceRestrictAddressing,
            const labelList& neighbourRestrictAddressing,
            const bool useGlobalPolyPatch = false
        );

        //- Construct as copy
        AMIInterpolationS4F
        (
            const AMIInterpolationS4F& ami,
            const bool useGlobalPolyPatch = false
        );

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolationS4F> clone() const
        {
            return autoPtr<AMIInterpolationS4F>::New(*this);
        }


    //- Destructor
    virtual ~AMIInterpolationS4F() = default;


    // Member Functions

        // Access

            //- Access to the up-to-date flag
            inline bool upToDate() const;

            //- Access to the up-to-date flag
            inline bool& upToDate();

            //- Access to the distributed flag
            inline bool distributed() const;

            //- Access to the requireMatch flag
            inline bool requireMatch() const;

            //- Access to the requireMatch flag
            inline bool setRequireMatch(const bool flag);

            //- Access to the reverseTarget flag
            inline bool reverseTarget() const;

            //- Threshold weight below which interpolation is deactivated
            inline scalar lowWeightCorrection() const;

            //- Return true if employing a 'lowWeightCorrection'
            inline bool applyLowWeightCorrection() const;

            //- Set to -1, or the processor holding all faces (both sides) of
            //- the AMI
            inline label singlePatchProc() const;

            //- PC: access to the useGlobalPolyPatch switch to allow the use
            //  of a globalPolyPatch. This should be called before the calculate
            //  function
            bool& useGlobalPolyPatch()
            {
                return useGlobalPolyPatch_;
            }

            // Source patch

                //- Return const access to source patch face areas
                inline const List<scalar>& srcMagSf() const;

                //- Return access to source patch face areas
                inline List<scalar>& srcMagSf();

                //- Return const access to source patch addressing
                inline const labelListList& srcAddress() const;

                //- Return access to source patch addressing
                inline labelListList& srcAddress();

                //- Return const access to source patch weights
                inline const scalarListList& srcWeights() const;

                //- Return access to source patch weights
                inline scalarListList& srcWeights();

                //- Return const access to normalisation factor of source
                //- patch weights (i.e. the sum before normalisation)
                inline const scalarField& srcWeightsSum() const;

                //- Return access to normalisation factor of source
                //- patch weights (i.e. the sum before normalisation)
                inline scalarField& srcWeightsSum();

                //- Return const access to source patch face centroids
                inline const pointListList& srcCentroids() const;

                //- Return access to source patch face centroids
                inline pointListList& srcCentroids();

                //- Source map pointer - valid only if singlePatchProc = -1
                //- This gets source data into a form to be consumed by
                //- tgtAddress, tgtWeights
                inline const mapDistribute& srcMap() const;


            // Target patch

                //- Return const access to target patch face areas
                inline const List<scalar>& tgtMagSf() const;

                //- Return access to target patch face areas
                inline List<scalar>& tgtMagSf();

                //- Return const access to target patch addressing
                inline const labelListList& tgtAddress() const;

                //- Return access to target patch addressing
                inline labelListList& tgtAddress();

                //- Return const access to target patch weights
                inline const scalarListList& tgtWeights() const;

                //- Return access to target patch weights
                inline scalarListList& tgtWeights();

                //- Return const access to normalisation factor of target
                //- patch weights (i.e. the sum before normalisation)
                inline const scalarField& tgtWeightsSum() const;

                //- Return access to normalisation factor of target
                //- patch weights (i.e. the sum before normalisation)
                inline scalarField& tgtWeightsSum();

                //- Target map pointer -  valid only if singlePatchProc=-1.
                //- This gets target data into a form to be consumed by
                //- srcAddress, srcWeights
                inline const mapDistribute& tgtMap() const;

            //***********************************************************//
            //  JN: from extendedggiinterpolation - possibly add here?
            //***********************************************************//

            //- Return distance to intersection for patch points
            //const scalarField& slavePointDistanceToIntersection() const;

            //- JN: this should be a private member function
            //- JN: OK to put it here?
            //- Calculate point weights
            //void calcSlavePointAddressing() const;

            //***********************************************************//

        // Manipulation

            //- Update addressing, weights and (optional) centroids
            virtual bool calculate
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch,
                const autoPtr<searchableSurface>& surfPtr = nullptr
            );

            //- Set the maps, addresses and weights from an external source
            void reset
            (
                autoPtr<mapDistribute>&& srcToTgtMap,
                autoPtr<mapDistribute>&& tgtToSrcMap,
                labelListList&& srcAddress,
                scalarListList&& srcWeights,
                labelListList&& tgtAddress,
                scalarListList&& tgtWeights
            );

            //- Append additional addressing and weights
            void append
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch
            );

            //- Normalise the weights
            void normaliseWeights(const bool conformal, const bool output);


        // Evaluation

            // Low-level

                //- Interpolate from target to source with supplied op
                //- to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToSource
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;

                //- Interpolate from source to target with supplied op
                //- to combine existing value with remote value and weight
                template<class Type, class CombineOp>
                void interpolateToTarget
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>::null()
                ) const;


            //- Interpolate from target to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field to source with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field to target with supplied op
            template<class Type, class CombineOp>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const CombineOp& cop,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target to source
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToSource
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source tmp field
            template<class Type>
            tmp<Field<Type>> interpolateToTarget
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;


        // Point intersections

            //- Return source patch face index of point on target patch face
            label srcPointFace
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch,
                const vector& n,
                const label tgtFacei,
                point& tgtPoint
            )
            const;

            //- Return target patch face index of point on source patch face
            label tgtPointFace
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch,
                const vector& n,
                const label srcFacei,
                point& srcPoint
            )
            const;


        // Checks

            //- Check if src addresses are present in tgt addresses and
            //- viceversa
            bool checkSymmetricWeights(const bool log) const;

            //- Write face connectivity as OBJ file
            void writeFaceConnectivity
            (
                const standAlonePatch& srcPatch,
                const standAlonePatch& tgtPatch,
                const labelListList& srcAddress
            ) const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIInterpolationS4FI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMIInterpolationS4FTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
