/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::globalPoints

Description
    Calculates points shared by more than two processor patches or cyclic
    patches.

    Is used in globalMeshData. (this info is needed for point-edge
    communication where you do all but these shared points with patch to
    patch communication but need to do a reduce on these shared points)

    Works purely topological and using local communication only.
    Needs:
      - domain to be one single domain (i.e. all faces can be reached through
        face-cell walk).
      - patch face ordering to be ok
      - f[0] ordering on patch faces to be ok.

    Works by constructing equivalence lists for all the points on processor
    patches. These list are procPointList and give processor and meshPoint
    label on that processor.
    E.g.
    @verbatim
          ((7 93)(4 731)(3 114))
    @endverbatim

    means point 93 on proc7 is connected to point 731 on proc4 and 114 on proc3.
    It then gets the lowest numbered processor (the 'master') to request a
    sharedPoint label from processor0 and it redistributes this label back to
    the other processors in the equivalence list.

    Algorithm:
        - get meshPoints of all my points on processor patches and initialize
          equivalence lists to this.
     loop
        - send to all neighbours in relative form:
            - patchFace
            - index in face
        - receive and convert into meshPoints. Add to to my equivalence lists.
        - mark meshPoints for which information changed.
        - send data for these meshPoints again
     endloop until nothing changes

    At this point one will have complete point-point connectivity for all
    points on processor patches. Now

        - remove point equivalences of size 2. These are just normal points
          shared between two neighbouring procPatches.
        - collect on each processor points for which it is the master
        - request number of sharedPointLabels from the Pstream::master.

    This information gets redistributed to all processors in a similar way
    as that in which the equivalence lists were collected:

        - initialize the indices of shared points I am the master for
     loop
        - send my known sharedPoints + meshPoints to all neighbours
        - receive from all neighbour. Find which meshPoint on my processor
          the sharedpoint is connected to
        - mark indices for which information has changed
     endloop until nothing changes.


SourceFiles
    globalPoints.C

\*---------------------------------------------------------------------------*/

#ifndef globalPoints_H
#define globalPoints_H

#include "DynamicList.H"
#include "Map.H"
#include "labelList.H"
#include "FixedList.H"
#include "primitivePatch.H"
#include "className.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyBoundaryMesh;
class cyclicPolyPatch;

/*---------------------------------------------------------------------------*\
                           Class globalPoints Declaration
\*---------------------------------------------------------------------------*/

class globalPoints
{
    // Private classes

        //- Define procPointList as holding a list of meshPoint/processor labels
        typedef FixedList<label, 2> procPoint;
        typedef List<procPoint> procPointList;

    // Private data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Sum of points on processor patches (unfiltered, point on 2 patches
        //  counts as 2)
        const label nPatchPoints_;

        //- All points on boundaries and their corresponding connected points
        //  on other processors.
        DynamicList<procPointList> procPoints_;

        //- Map from mesh point to index in procPoints
        Map<label> meshToProcPoint_;

        //- Shared points used by this processor (= global point number)
        labelList sharedPointAddr_;

        //- My meshpoints corresponding to the shared points
        labelList sharedPointLabels_;

        //- Total number of shared points.
        label nGlobalPoints_;


    // Private Member Functions

        //- Count all points on processorPatches. Is all points for which
        //  information is collected.
        static label countPatchPoints(const polyBoundaryMesh&);

        //- Add information about patchPointI in relative indices to send
        //  buffers (patchFaces, indexInFace etc.)
        static void addToSend
        (
            const primitivePatch&,
            const label patchPointI,
            const procPointList&,
            DynamicList<label>& patchFaces,
            DynamicList<label>& indexInFace,
            DynamicList<procPointList>& allInfo
        );

        //- Merge info from neighbour into my data
        static bool mergeInfo
        (
            const procPointList& nbrInfo,
            procPointList& myInfo
        );

        //- Store (and merge) info for meshPointI
        bool storeInfo(const procPointList& nbrInfo, const label meshPointI);

        //- Initialize procPoints_ to my patch points. allPoints = true:
        //  seed with all patch points, = false: only boundaryPoints().
        void initOwnPoints(const bool allPoints, labelHashSet& changedPoints);

        //- Send subset of procPoints to neighbours
        void sendPatchPoints(const labelHashSet& changedPoints) const;

        //- Receive neighbour points and merge into my procPoints.
        void receivePatchPoints(labelHashSet& changedPoints);

        //- Remove entries of size 2 where meshPoint is in provided Map.
        //  Used to remove normal face-face connected points.
        void remove(const Map<label>&);

        //- Get indices of point for which I am master (lowest numbered proc)
        labelList getMasterPoints() const;

        //- Send subset of shared points to neighbours
        void sendSharedPoints(const labelList& changedIndices) const;

        //- Receive shared points and update subset.
        void receiveSharedPoints(labelList& changedIndices);


        //- Should move into cyclicPolyPatch but some ordering problem
        //  keeps on giving problems.
        static edgeList coupledPoints(const cyclicPolyPatch&);

        //- Disallow default bitwise copy construct
        globalPoints(const globalPoints&);

        //- Disallow default bitwise assignment
        void operator=(const globalPoints&);


public:

        //- Declare name of the class and its debug switch
        ClassName("globalPoints");


    // Constructors

        //- Construct from mesh
        globalPoints(const polyMesh& mesh);


    // Member Functions

        // Access

            label nPatchPoints() const
            {
                return nPatchPoints_;
            }

            const DynamicList<procPointList>& procPoints() const
            {
                return procPoints_;
            }

            const Map<label>& meshToProcPoint() const
            {
                return meshToProcPoint_;
            }

            //- shared points used by this processor (= global point number)
            const labelList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- my meshpoints corresponding to the shared points
            const labelList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
