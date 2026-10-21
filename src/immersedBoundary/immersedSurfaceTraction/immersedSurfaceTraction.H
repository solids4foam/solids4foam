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
    Foam::immersedSurfaceTraction

Description
    Traction exerted by the fluid on the surface of an immersed body,
    evaluated at quadrature points that move with the surface.

    The quadrature points are placed on the surface triangles at a spacing of
    about the local cell width h: in three dimensions at the centroids of a
    uniform subdivision of each triangle, and in two dimensions (one empty
    direction) along the intersection of the surface with the mid-plane of
    the mesh, where each point represents the thickness of the mesh. Points
    outside the mesh are not used. Each point is defined by its triangle and
    barycentric coordinates, so that it follows the motion or deformation of
    the surface. The candidate points are used with one of three
    quadratures:

    - clustered (the default): the candidates are grouped into clusters,
      each represented by one point (see below), so that the quadrature
      does not depend on the order of the triangles;
    - none: every candidate point;
    - greedy (legacy): where the triangles are smaller than the cells, the
      candidates are thinned to a spacing of about 0.7 mean cell widths,
      kept in the order of the surface triangles, with the area of each
      removed point given to the nearest kept point. The kept points depend
      on the order of the triangles, and the area that each takes from its
      neighbours is not centred on it, so that the quadrature of a smooth
      traction has a first-order error (e.g. a lift whose sign depends on
      the order of the triangles); it is kept to reproduce earlier results.

    A body driven by a fluid-solid interface always uses none: the traction
    is returned to the solid face by face (faceValues), and each face then
    needs its own points.

    With the quadrature clustered, the candidate points are made in the
    reference configuration of the surface (immersedBody::points0) and
    grouped by the cell of a uniform grid of spacing s, the dominant signed
    axis of their normal (so that the two sides of a thin body are never
    merged) and whether they are inside the mesh, and each such group is
    split into the connected components of the triangles of its candidates
    (so that the layers of a folded surface or the faces of nearby bodies
    are never merged). The spacing s is a fraction tractionPointSpacing of
    the smallest width of the cells containing the candidates in the
    reference configuration, and the grid has its origin at the origin of
    the coordinates. The clusters are therefore a function of the reference
    surface, the mesh and tractionPointSpacing only: they do not depend on
    the order of the triangles, the decomposition or the time at which they
    are made (e.g. a restart), provided that the mesh is the same (a mesh
    whose cells change, e.g. by adaptive refinement, can give different
    clusters). The spacing is not changed as the body moves: a body moving
    into finer cells keeps the clusters of its reference configuration.

    At each time each cluster is represented by one point, its current area
    centroid projected on the current plane of the triangle of its candidate
    nearest to the reference centroid, with the area A = sum(A_i) and the
    vector area a = sum(A_i n_i) of its candidates on their current
    triangles, and the normal a/|a|. The force integrates the tangential
    traction with A and the pressure with the vector area a: the total area
    is kept, the pressure force of a closed surface is invariant to a
    constant pressure (the vector areas sum to that of the candidates), and
    the first moment of the area of each cluster about its point vanishes in
    the tangent plane. The traction at the points (traction) is the physical
    traction, -p n + tau, with the unit normal n = a/|a|: the vector area is
    used in the force only. faceValues gives the value of a cluster to each
    of its candidates, with its area.

    At each point \f$ x_w \f$, with the outward unit normal n, the pressure
    and the velocity relative to the rigid body velocity \f$ U_b \f$ are
    fitted by weighted quadratic least squares, centred on the point,
    over the cells within 3 h of it. The weights decrease with the distance
    to the point, and are zero for the cells behind the tangent plane of the
    surface, rising smoothly to their full value at h/2 in front of it, so
    that the fits use fluid cells only and vary continuously as the body
    moves.
    The velocity fit has no constant term, so that the velocity equals the
    body velocity at the point. The traction exerted by the fluid on the
    surface is

    \f[
        t = \rho \left[ -p_w n + \nu (I - n n)
            \cdot \left( \nabla (U - U_b) \cdot n \right) \right]
    \f]

    with the wall pressure \f$ p_w \f$ and the velocity gradient from the
    fits, where p is the kinematic pressure: on a no-slip wall the tangential
    derivatives of \f$ U - U_b \f$ vanish, and the rigid body velocity has
    no strain. The fits are solved by Cholesky decomposition, or with the
    pseudo-inverse of their matrices for rank-deficient stencils.

    The sums of the fits over the processors are reduced, so that the cells
    of a fit may be on any processor, and each fit is then solved on one
    processor. Points outside the mesh at a given time are not used. The
    stencils and widths are found again when the body or the mesh moves.

Author
    Philip Cardiff, UCD.

SourceFiles
    immersedSurfaceTraction.C

\*---------------------------------------------------------------------------*/

#ifndef immersedSurfaceTraction_H
#define immersedSurfaceTraction_H

#include "immersedBody.H"
#include "FixedList.H"
#include "Enum.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedSurfaceTraction Declaration
\*---------------------------------------------------------------------------*/

class immersedSurfaceTraction
{
public:

    // Public Data Types

        //- Quadratures of the surface traction
        enum quadratureType
        {
            greedy,     //!< Thinning in the order of the triangles (legacy)
            none,       //!< Every candidate point
            clustered   //!< One point per cluster of candidate points
        };

        //- Names of the quadratures
        static const Enum<quadratureType> quadratureTypeNames;


private:

    // Private Data

        //- Reference to the fluid mesh
        const fvMesh& mesh_;

        //- The body
        const immersedBody& body_;

        //- Solved directions
        const Vector<label> solD_;

        //- Quadrature
        const quadratureType quadrature_;

        //- Size of the clusters relative to the smallest width of the cells
        //  containing the surface (clustered)
        const scalar spacingFactor_;

        //- Parent surface face of each point
        labelList faces_;

        //- Barycentric coordinates of each point in its face
        List<FixedList<scalar, 3>> barycentric_;

        //- Surface area represented by each point (including the mesh
        //  thickness in two dimensions)
        scalarField areas_;

        //- Initial area of each point, from which the area of a point of a
        //  deforming body is scaled (greedy and none)
        scalarField areas0_;

        //- Initial area of the parent face of each point (greedy and none)
        scalarField faceAreas0_;

        //- Local cell width at each point
        scalarField h_;

        //- Current position of each point
        pointField points_;

        //- Current outward unit normal at each point
        vectorField normals_;

        //- Local cells of the least squares stencil of each point
        List<labelList> stencils_;

        //- Parent surface face of each candidate point (clustered)
        labelList candidateFaces_;

        //- Barycentric coordinates of each candidate point in its face
        //  (clustered)
        List<FixedList<scalar, 3>> candidateBarycentric_;

        //- Area of each candidate point in the reference configuration
        //  (clustered)
        scalarField candidateAreas0_;

        //- Current area of each candidate point (clustered)
        scalarField candidateAreas_;

        //- Area of the parent face of each candidate point in the reference
        //  configuration (clustered)
        scalarField candidateFaceAreas0_;

        //- Candidate points of each cluster, i.e. of each point (clustered)
        List<labelList> members_;

        //- Ratio |a|/A of the magnitude of the vector area to the area of
        //  each point, with which the force integrates its pressure: 1
        //  unless clustered
        scalarField pressureScale_;


    // Private Member Functions

        //- Width of the cell containing each point, the largest over the
        //  processors, and 0 if no processor contains it
        scalarField cellWidths(const pointField& pts) const;

        //- Create the quadrature points
        void createPoints();

        //- Update the positions and normals of the points
        void updatePoints();

        //- Wall pressure and viscous traction at each point, from the least
        //  squares fits (zero at the points outside the mesh or without a
        //  fit), the same on every processor
        void wallValues
        (
            const volVectorField& U,
            const volScalarField& p,
            const scalar nu,
            scalarField& pw,
            vectorField& tau
        );


public:

    // Constructors

        //- Construct from the body, the mesh, the quadrature and the size of
        //  the clusters relative to the smallest width of the cells
        //  containing the surface (clustered only)
        immersedSurfaceTraction
        (
            const immersedBody& body,
            const fvMesh& mesh,
            const quadratureType quadrature = clustered,
            const scalar spacingFactor = 0.5
        );

        //- No copy construct
        immersedSurfaceTraction(const immersedSurfaceTraction&) = delete;

        //- No copy assignment
        void operator=(const immersedSurfaceTraction&) = delete;


    //- Destructor
    ~immersedSurfaceTraction() = default;


    // Member Functions

        //- Number of quadrature points
        label size() const
        {
            return faces_.size();
        }

        //- Current quadrature points
        const pointField& points() const
        {
            return points_;
        }

        //- Kinematic traction (per unit density) exerted by the fluid on the
        //  surface at each quadrature point: -p n + tau, with the wall
        //  pressure p, the unit normal n and the viscous traction tau at the
        //  point, and zero at the points outside the mesh
        tmp<vectorField> traction
        (
            const volVectorField& U,
            const volScalarField& p,
            const scalar nu
        );

        //- Area-weighted mean of the values at the quadrature points over
        //  each surface face, where the points outside the mesh have no
        //  value but keep their area, and the area of the quadrature points
        //  of each face: zero for a face without quadrature points. With the
        //  clustered quadrature, the value of a cluster is given to each of
        //  its candidates, with its area: the face values of the traction
        //  then integrate the pressure with the area, not the vector area,
        //  of the cluster (see force)
        void faceValues
        (
            const vectorField& pointValues,
            vectorField& faceValues,
            scalarField& faceAreas
        ) const;

        //- Force and torque (about CofR) exerted by the fluid on the body:
        //  the sum over the points of rho*(-p |a| n + A tau), with the area
        //  A and the vector area a of each point (|a| = A unless clustered)
        void force
        (
            const volVectorField& U,
            const volScalarField& p,
            const scalar nu,
            const scalar rho,
            const point& CofR,
            vector& F,
            vector& T
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
