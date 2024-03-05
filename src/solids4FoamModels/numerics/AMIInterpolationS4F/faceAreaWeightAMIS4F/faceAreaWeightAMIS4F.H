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
    Foam::faceAreaWeightAMIS4F

Description
    Face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    faceAreaWeightAMIS4F.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightAMIS4F_H
#define faceAreaWeightAMIS4F_H

#include "advancingFrontAMIS4F.H"
#include "standAlonePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faceAreaWeightAMIS4F Declaration
\*---------------------------------------------------------------------------*/

class faceAreaWeightAMIS4F
:
    public advancingFrontAMIS4F
{
private:

    // Private data

        //- Flag to restart uncovered source faces
        const bool restartUncoveredSourceFace_;



protected:

    // Protected Member Functions

        //- No copy assignment
        void operator=(const faceAreaWeightAMIS4F&) = delete;

        //- Initialise the geometry
        void initGeom
        (
            const standAlonePatch& srcPatch,
            const standAlonePatch& tgtPatch,
            const globalIndex& globalTgtFaces,
            labelList& extendedTgtFaceIDs
        );


        // Marching front

            //- Calculate addressing, weights and centroids using temporary
            //- storage
            virtual void calcAddressing
            (
                List<DynamicList<label>>& srcAddress,
                List<DynamicList<scalar>>& srcWeights,
                List<DynamicList<point>>& srcCentroids,
                List<DynamicList<label>>& tgtAddress,
                List<DynamicList<scalar>>& tgtWeights,
                label srcFacei,
                label tgtFacei
            );

            //- Determine overlap contributions for source face srcFacei
            virtual bool processSourceFace
            (
                const label srcFacei,
                const label tgtStartFacei,
                DynamicList<label>& nbrFaces,
                DynamicList<label>& visitedFaces,
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght
            );

            //- Attempt to re-evaluate source faces that have not been included
            virtual void restartUncoveredSourceFace
            (
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght
            );

            //- Set the source and target seed faces
            virtual bool setNextFaces
            (
                label& startSeedi,
                label& srcFacei,
                label& tgtFacei,
                const bitSet& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces,
                const bool errorOnNotFound = true
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            virtual void calcInterArea
            (
                const label srcFacei,
                const label tgtFacei,
                scalar& area,
                vector& centroid
            ) const;

            //- Return true if faces overlap
            virtual bool overlaps
            (
                const label srcFacei,
                const label tgtFacei,
                const scalar threshold
            ) const;


public:

    //- Runtime type information
    TypeName("faceAreaWeightAMIS4F");


    // Constructors

        //- Construct from dictionary
        faceAreaWeightAMIS4F
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        faceAreaWeightAMIS4F
        (
            const bool requireMatch,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const faceAreaIntersect::triangulationMode triMode =
                faceAreaIntersect::tmMesh,
            const bool restartUncoveredSourceFace = true
        );

        //- Construct as copy
        faceAreaWeightAMIS4F(const faceAreaWeightAMIS4F& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolationS4F> clone() const
        {
            return autoPtr<AMIInterpolationS4F>(new faceAreaWeightAMIS4F(*this));
        }


    //- Destructor
    virtual ~faceAreaWeightAMIS4F() = default;


    // Member Functions

        //- Update addressing, weights and (optional) centroids
        virtual bool calculate
        (
            const standAlonePatch& srcPatch,
            const standAlonePatch& tgtPatch,
            const autoPtr<searchableSurface>& surfPtr = nullptr
        );

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
