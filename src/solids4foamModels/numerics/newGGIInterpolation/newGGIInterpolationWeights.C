/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    Mass-conservative face interpolation of face data between two
    primitivePatches

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Modification by:
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "newGGIInterpolationTemplate.H"
#include "objectHit.H"
#include "boolList.H"
#include "DynamicList.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::debug::tolerancesSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::areaErrorTol_
(
    "GGIAreaErrorTol",
    1.0e-8,
    "Minimum GGI face to face intersection area. "
    "The smallest accepted GGI weighting factor."
);

template<class MasterPatch, class SlavePatch>
const Foam::debug::tolerancesSwitch
newGGIInterpolation<MasterPatch, SlavePatch>::featureCosTol_
(
    "GGIFeatureCosTol",
    0.8,
    "Minimum cosine value between 2 GGI patch neighbouring facet normals."
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::calcAddressing() const
{
    if
    (
        masterAddrPtr_
     || masterWeightsPtr_
     || slaveAddrPtr_
     || slaveWeightsPtr_
     || uncoveredMasterAddrPtr_
     || uncoveredSlaveAddrPtr_
    )
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "calcAddressing() const"
        )   << "Evaluation of GGI weighting factors:" << endl;
    }

    // Create the dynamic lists to hold the addressing

    // The final master/slave list, after filtering out the "false" neighbours
    List<DynamicList<label, 8> > masterNeighbors(masterPatch_.size());
    List<DynamicList<label, 8> > slaveNeighbors(slavePatch_.size());

    // The weights
    List<DynamicList<scalar, 8> > masterNeighborsWeights(masterPatch_.size());
    List<DynamicList<scalar, 8> > slaveNeighborsWeights(slavePatch_.size());

    // First, find a rough estimate of each slave and master facet
    // neighborhood by filtering out all the faces located outside of
    // an Axis-Aligned Bounding Box (AABB).  Warning: This algorithm
    // is based on the evaluation of AABB boxes, which is pretty fast;
    // but still the complexity of the algorithm is n^2, wich is
    // pretty bad for GGI patches composed of 100,000 of facets...  So
    // here is the place where we could certainly gain major speedup
    // for larger meshes.

    // The candidates master neighbours
    // Choice of algorithm:
    // 1) Axis-aligned bounding box
    // 2) Octree search with bounding box
    // 3) 3-D vector distance
    // 4) n-Squared search
    labelListList candidateMasterNeighbors;

    if (reject_ == AABB)
    {
         findNeighboursAABB(candidateMasterNeighbors);
    }
    else if (reject_ == BB_OCTREE)
    {
         findNeighboursBBOctree(candidateMasterNeighbors);
    }
    else if (reject_ == THREE_D_DISTANCE)
    {
         findNeighbours3D(candidateMasterNeighbors);
    }
    else if (reject_ == N_SQUARED)
    {
        candidateMasterNeighbors.setSize(masterPatch_.size());

        // Mark N-squared search
        labelList nSquaredList(slavePatch_.size());
        forAll (nSquaredList, i)
        {
            nSquaredList[i] = i;
        }

        forAll (candidateMasterNeighbors, j)
        {
            candidateMasterNeighbors[j] = nSquaredList;
        }
    }

    // Next, we move to the 2D world.  We project each slave and
    // master face onto a local plane defined by the master face
    // normal.  We filter out a few false neighbors using the
    // Separating Axes Theorem

    // It is in this local plane that we will refine our list of
    // neighbors.  So for a given a neighbor face, we need as many
    // projections as there are neighbors closeby.

    const pointField& masterPatchPoints = masterPatch_.points();
    const vectorField masterPatchNormals = masterPatch_.faceNormals();

    // Store the polygon made by projecting the face points onto the
    // face normal
    // The master faces polygons
    List<pointField>  masterFace2DPolygon(masterPatch_.size());

    // Tolerance factor for the Separation of Axes Theorem == distErrorTol_

    forAll (masterPatch_, faceMi)
    {
        // First, we make sure that all the master faces points are
        // recomputed onto the 2D plane defined by the master faces
        // normals.
        // For triangles, this is useless, but for N-gons
        // with more than 3 points, this is essential.
        // The intersection between the master and slave faces will be
        // done in these 2D reference frames

        // A few basic information to keep close-by
        vector currentMasterFaceNormal = masterPatchNormals[faceMi];
        vector currentMasterFaceCentre =
            masterPatch_[faceMi].centre(masterPatchPoints);

        scalarField facePolygonErrorProjection;

        // Project the master faces points onto the normal face plane to
        // form a flattened polygon
        masterFace2DPolygon[faceMi] =
            projectPointsOnPlane
            (
                masterPatch_[faceMi].points(masterPatchPoints),
                currentMasterFaceCentre,
                currentMasterFaceNormal,
                facePolygonErrorProjection
            );

        // Next we compute an orthonormal basis (u, v, w) aligned with
        // the face normal for doing the 3D to 2D projection.
        //
        // "w" is aligned on the face normal.  We need to select a "u"
        // direction, it can be anything as long as it lays on the
        // projection plane.  We chose to use the direction from the
        // master face center to the most distant projected master face
        // point on the plane.  Finally, we get "v" by evaluating the
        // cross-product w^u = v.  And we make sure that u, v, and w are
        // normalized.
        //
        //
        //                                                       u  =  vector from face center to most distant projected master face point.
        //                                               /       .
        //           ^y                                / |       .      .w = normal to master face
        //           |                               /   |       .    .
        //           |                             /     |       .  .
        //           |                            |      |       .
        //           |                            |      /        .
        //           |                            |    /           .
        //           |                            |  /              .
        //           ---------> x                 |/                 .
        //          /                                                 v = w^u
        //         /
        //        /
        //       z
        //
        //

        orthoNormalBasis uvw =
            computeOrthonormalBasis
            (
                currentMasterFaceCentre,
                currentMasterFaceNormal,
                masterFace2DPolygon[faceMi]
            );

        // Recompute the master polygon into this orthoNormalBasis
        // We should only see a rotation along the normal of the face here
        List<point2D> masterPointsInUV;
        scalarField masterErrorProjectionAlongW;

        masterPointsInUV =
            projectPoints3Dto2D
            (
                uvw,
                currentMasterFaceCentre,
                masterFace2DPolygon[faceMi],
                masterErrorProjectionAlongW   // Should be at zero all the way
            );

        // Compute the surface area of the polygon;
        // We need this for computing the weighting factors
        scalar surfaceAreaMasterPointsInUV = area2D(masterPointsInUV);

        // Check if polygon is CW.. Should not, it should be CCW; but
        // better and cheaper to check here
        if (surfaceAreaMasterPointsInUV < 0.0)
        {
            reverse(masterPointsInUV);
            surfaceAreaMasterPointsInUV = -surfaceAreaMasterPointsInUV;

            // Just generate a warning until we can verify this is a non issue
            InfoIn
            (
                "void newGGIInterpolation<MasterPatch, SlavePatch>::"
                "calcAddressing()"
            )   << "The master projected polygon was CW instead of CCW.  "
                << "This is strange..."  << endl;
        }

        // Next, project the candidate master neighbours faces points
        // onto the same plane using the new orthonormal basis
        const labelList& curCMN = candidateMasterNeighbors[faceMi];

        forAll (curCMN, neighbI)
        {
            // For each points, compute the dot product with u,v,w.  The
            // [u,v] component will gives us the 2D cordinates we are
            // looking for for doing the 2D intersection The w component
            // is basically the projection error normal to the projection
            // plane

            // NB: this polygon is most certainly CW w/r to the uvw
            // axis because of the way the normals are oriented on
            // each side of the GGI interface... We will switch the
            // polygon to CCW in due time...
            List<point2D> neighbPointsInUV;
            scalarField neighbErrorProjectionAlongW;

            // We use the xyz points directly, with a possible transformation
            pointField curSlaveFacePoints =
                slavePatch_[curCMN[neighbI]].points(slavePatch_.points());

            if (doTransform())
            {
                // Transform points to master plane
                if (forwardT_.size() == 1)
                {
                    transform
                    (
                        curSlaveFacePoints,
                        forwardT_[0],
                        curSlaveFacePoints
                    );
                }
                else
                {
                    transform
                    (
                        curSlaveFacePoints,
                        forwardT_[curCMN[neighbI]],
                        curSlaveFacePoints
                    );
                }
            }

            // Apply the translation offset in order to keep the
            // neighbErrorProjectionAlongW values to a minimum
            if (doSeparation())
            {
                if (forwardSep_.size() == 1)
                {
                    curSlaveFacePoints += forwardSep_[0];
                }
                else
                {
                    curSlaveFacePoints += forwardSep_[curCMN[neighbI]];
                }
            }

            neighbPointsInUV =
                projectPoints3Dto2D
                (
                    uvw,
                    currentMasterFaceCentre,
                    curSlaveFacePoints,
                    neighbErrorProjectionAlongW
                );

            // We are now ready to filter out the "bad" neighbours.
            // For this, we will apply the Separating Axes Theorem
            // http://en.wikipedia.org/wiki/Separating_axis_theorem.

            // This will be the second and last quick reject test.
            // We will use the 2D projected points for both the master
            // patch and its neighbour candidates
            if
            (
                detect2dPolygonsOverlap
                (
                    masterPointsInUV,
                    neighbPointsInUV,
                    sqrt(areaErrorTol_()) // distErrorTol
                )
            )
            {
                // We have an overlap between the master face and this
                // neighbor face.
                label faceMaster = faceMi;
                label faceSlave  = curCMN[neighbI];

                // Compute the surface area of the neighbour polygon;
                // We need this for computing the weighting factors
                scalar surfaceAreaNeighbPointsInUV = area2D(neighbPointsInUV);

                // Check for CW polygons. It most certainly is, and
                // the polygon intersection algorithms are expecting
                // to work with CCW point ordering for the polygons
                if (surfaceAreaNeighbPointsInUV < 0.0)
                {
                    reverse(neighbPointsInUV);
                    surfaceAreaNeighbPointsInUV = -surfaceAreaNeighbPointsInUV;
                }


                // We compute the intersection area using the
                // Sutherland-Hodgman algorithm.  Of course, if the
                // intersection area is 0, that would constitute the last and
                // final reject test, but it would also be an indication that
                // our 2 previous rejection tests are a bit laxed...  or that
                // maybe we are in presence of concave polygons....
                scalar intersectionArea =
                    polygonIntersection
                    (
                        masterPointsInUV,
                        neighbPointsInUV
                    );

                if (intersectionArea > VSMALL) // Or > areaErrorTol_ ???
                {
                    // We compute the GGI weights based on this
                    // intersection area, and on the individual face
                    // area on each side of the GGI.

                    // Since all the intersection have been computed
                    // in the projected UV space we need to compute
                    // the weights using the surface area from the
                    // faces projection as well. That way, we make
                    // sure all our factors will sum up to 1.0.

                    masterNeighbors[faceMaster].append(faceSlave);
                    slaveNeighbors[faceSlave].append(faceMaster);

                    masterNeighborsWeights[faceMaster].append
                    (
                        intersectionArea/surfaceAreaMasterPointsInUV
                    );

                    slaveNeighborsWeights[faceSlave].append
                    (
                        intersectionArea/surfaceAreaNeighbPointsInUV
                    );
                }
                else
                {
                    WarningIn
                    (
                        "newGGIInterpolation<MasterPatch, SlavePatch>::"
                        "calcAddressing()"
                    )   << "polygonIntersection is returning a "
                        << "zero surface area between " << nl
                        << "     Master face: " << faceMi
                        << " and Neighbour face: " << curCMN[neighbI]
                        << " intersection area = " << intersectionArea << nl
                        << "Please check the two quick-check algorithms for "
                        << "newGGIInterpolation.  Something is  missing."
                        << endl;
                }
            }
        }

        // We went through all the possible neighbors for this face.
    }


    // Allocate the member attributes and pack addressing
    masterAddrPtr_ = new labelListList(masterPatch_.size());
    labelListList& ma  = *masterAddrPtr_;

    masterWeightsPtr_ = new scalarListList(masterPatch_.size());
    scalarListList& maW = *masterWeightsPtr_;

    forAll (ma, mfI)
    {
        ma[mfI].transfer(masterNeighbors[mfI].shrink());
        maW[mfI].transfer(masterNeighborsWeights[mfI].shrink());
    }

    slaveAddrPtr_ = new labelListList(slavePatch_.size());
    labelListList& sa = *slaveAddrPtr_;

    slaveWeightsPtr_ = new scalarListList(slavePatch_.size());
    scalarListList& saW = *slaveWeightsPtr_;

    forAll (sa, sfI)
    {
        sa[sfI].transfer(slaveNeighbors[sfI].shrink());
        saW[sfI].transfer(slaveNeighborsWeights[sfI].shrink());
    }

    // Now that the neighbourhood is known, let's go hunting for
    // non-overlapping faces
    uncoveredMasterAddrPtr_ =
        new labelList
        (
            findNonOverlappingFaces(maW, masterNonOverlapFaceTol_)
        );

    uncoveredSlaveAddrPtr_ =
        new labelList
        (
            findNonOverlappingFaces(saW, slaveNonOverlapFaceTol_)
        );

    // Rescaling the weighting factors so they will sum up to 1.0
    // See the comment for the method ::rescaleWeightingFactors() for
    // more information.  By default, we always rescale.  But for some
    // special kind of GGI interpolation, like the mixingPlaneGGI,
    // then we need the brute values, so no rescaling in that
    // case. Hence the little flag rescaleGGIWeightingFactors_

    if (rescaleGGIWeightingFactors_)
    {
        rescaleWeightingFactors();
    }
}


// Rescaling the weighting factors so they will sum up to 1.0 This is
// necessary for the slave weighting factors because intersection with
// master neighbours are usually computed from different projection
// plane, so the weighting factor don't quite sum up to 1.0 For the
// slave weighting factors, we are usually talking of delta of the
// order of 10e-6 here.
//
// For the master weighting factor, this is another story. We truly
// expect that the master weighting factors will exactly sum up to 1.0
// if all the neighbours are properly identified.
//
// However, for concentric circular geometry, if the circumferantial
// resolution is too coarse, we will end up with some part of the face
// surface that are not taken into account because they do not
// physically overlap any neighbours.  For example, think of 2
// concentric circular patches, slightly rotated one relatively to the
// other.  A good case: the ercoftac conical diffuser, Case0...  GGI
// located right between the cylindrical and conical parts, rotate the
// cylindrical by 15 degrees.  For this case, we will need to devise a
// decent strategy in order to intelligently take care of these
// "missing weights"
//
// The purpose of the ::rescaleWeightingFactors() method is mainly for
// this.
template<class MasterPatch, class SlavePatch>
void newGGIInterpolation<MasterPatch, SlavePatch>::
rescaleWeightingFactors() const
{
    scalarListList& maW = *masterWeightsPtr_;
    scalarListList& saW = *slaveWeightsPtr_;

    // Memorize the largest correction needed in order to provide some
    // basic info to the user
    scalar largestSWC = 0;
    scalar sumSWC = 0;
    scalar curSWC = 0;

    scalar largestMWC = 0;
    scalar sumMWC = 0;
    scalar curMWC = 0;

    // Rescaling the slave weights
    if (debug)
    {
        if
        (
            uncoveredMasterFaces().size() > 0
         || uncoveredSlaveFaces().size() > 0
        )
        {
            InfoIn
            (
                "void newGGIInterpolation<MasterPatch, SlavePatch>::"
                "rescaleWeightingFactors() const"
            )   << "Uncovered faces found.  On master: "
                << uncoveredMasterFaces().size()
                << " on slave: " << uncoveredSlaveFaces().size() << endl;
        }
      }

    forAll (saW, saWi)
    {
        scalar slaveWeightSum = Foam::sum(saW[saWi]);

        if (saW[saWi].size() > 0)
        {
            saW[saWi] = saW[saWi]/slaveWeightSum;

            // Some book-keeping
            curSWC = mag(1.0 - slaveWeightSum);
            largestSWC = Foam::max(largestSWC, curSWC);

            sumSWC += curSWC;
        }
    }

    // Rescaling the master weights
    forAll (maW, maWi)
    {
        scalar masterWeightSum = Foam::sum(maW[maWi]);

        if (maW[maWi].size() > 0)
        {
            maW[maWi] = maW[maWi]/masterWeightSum;

            // Some book-keeping
            curMWC = mag(1.0 - masterWeightSum);
            largestMWC = Foam::max(largestMWC, curMWC);

            sumMWC += curMWC;
        }
    }

    if (debug)
    {
        if (saW.size() > 0 && maW.size() > 0)
        {
            Info<< "  Largest slave weighting factor correction : "
                << largestSWC
                << " average: " << sumSWC/saW.size() << nl
                << "  Largest master weighting factor correction: "
                << largestMWC
                << " average: " << sumMWC/maW.size() << endl;
        }
    }
}


// Find non-overlapping faces from both master and slave patches
// The default non-overlapping criteria is total absence of neighbours.
// Later on, ths criteria will be based on minimum surface intersection, or
// minimum weight factor
template<class MasterPatch, class SlavePatch>
tmp<labelField>
newGGIInterpolation<MasterPatch, SlavePatch>::findNonOverlappingFaces
(
    const scalarListList& patchWeights,
    const scalar& nonOverlapFaceTol   //  = min sum of the neighbour weights
) const
{
    tmp<labelField> tpatchFaceNonOverlapAddr(new labelField());
    labelField& patchFaceNonOverlapAddr = tpatchFaceNonOverlapAddr();

    DynamicList<label, 64> patchFaceNonOverlap(patchWeights.size());

    // Scan the list of patch weights, looking for empty lists
    forAll (patchWeights, paWi)
    {
        scalar sumWeightsFace = sum(patchWeights[paWi]);

        if (sumWeightsFace <= nonOverlapFaceTol)
        {
            // Store local index.
            patchFaceNonOverlap.append(paWi);
        }
    }

    if (patchFaceNonOverlap.size() > 0)
    {
        patchFaceNonOverlapAddr.transfer(patchFaceNonOverlap.shrink());
    }

    if (debug)
    {
        InfoIn("newGGIInterpolation::findNonOverlappingFaces")
            << "   : Found " << patchFaceNonOverlapAddr.size()
            << " non-overlapping faces for this GGI patch" << endl;
    }

    return tpatchFaceNonOverlapAddr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
