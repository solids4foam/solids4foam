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
    Foam::immersedMarkers

Description
    Lagrangian markers on the surface of an immersed body, with the
    interpolation of cell fields to the markers and the spreading of marker
    forces to the cells, for the marker (multi-direct forcing) immersed
    boundary method.

    The markers are placed on the surface triangles at a spacing of about
    spacing times the local cell width: in three dimensions at the centroids
    of a uniform subdivision of each triangle, and in two dimensions (one
    empty direction) along the intersection of the surface with the mid-plane
    of the mesh, where each marker represents the thickness of the mesh. Each
    marker is defined by its triangle and barycentric coordinates, so that it
    follows the motion or deformation of the surface.

    The interpolation weights of a marker are the Roma et al. (1999)
    three-point kernel, with a support of 1.5 times the local cell width h in
    each solved direction, weighted by the cell volumes and corrected to
    reproduce linear fields (a reproducing kernel), so that they sum to one
    and have zero first moments on general meshes:

    \f[
        w_{mc} = \tilde{w}_{mc} \, p(x_c - X_m) \cdot M^{-1} e_0, \qquad
        M = \sum_c \tilde{w}_{mc} \, p \, p^T
    \f]

    with \f$ p = (1, (x - X_m)/h) \f$. Spreading uses the same weights, so
    that it is the adjoint of the interpolation and conserves the force:
    \f$ f_c = \sum_m w_{mc} F_m \Delta V_m / V_c \f$, with the marker volume
    \f$ \Delta V_m = A_m h_m \f$.

    The cells of the kernel support may be on any processor: the moments
    and interpolated values of each marker are summed over the processors.

    References:
    \verbatim
        Roma, A.M., Peskin, C.S., Berger, M.J. (1999). An adaptive version
        of the immersed boundary method. Journal of Computational Physics,
        153, 509-534.
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    immersedMarkers.C

\*---------------------------------------------------------------------------*/

#ifndef immersedMarkers_H
#define immersedMarkers_H

#include "immersedBody.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class immersedMarkers Declaration
\*---------------------------------------------------------------------------*/

class immersedMarkers
{
    // Private Data

        //- Reference to the fluid mesh
        const fvMesh& mesh_;

        //- The body
        const immersedBody& body_;

        //- Marker spacing, relative to the local cell width
        const scalar spacing_;

        //- Retraction of the forcing points into the body, relative to the
        //  local cell width
        const scalar retraction_;

        //- Solved directions
        const Vector<label> solD_;

        //- Number of solved directions
        const label nD_;

        //- Parent surface face of each marker
        labelList faces_;

        //- Barycentric coordinates of each marker in its face
        List<FixedList<scalar, 3>> barycentric_;

        //- Surface area represented by each marker (including the mesh
        //  thickness in two dimensions)
        scalarField areas_;

        //- Current forcing points
        pointField points_;

        //- Current surface points
        pointField surfacePoints_;

        //- Current outward unit normals
        vectorField normals_;

        //- Current body velocity at the surface points
        vectorField velocities_;

        //- Kernel width of each marker
        scalarField h_;

        //- Marker volumes, area*h
        scalarField volumes_;

        //- Local kernel stencil of each marker (compressed rows)
        labelList offsets_;
        labelList cells_;
        scalarField weights_;


    // Private Member Functions

        //- Width of the cell containing each point, the largest over the
        //  processors, and 0 if no processor contains it
        scalarField cellWidths(const pointField& pts) const;

        //- Width of a cell: its smallest bounding box width in the solved
        //  directions
        scalar cellWidth(const label celli) const;

        //- Create the markers
        void createMarkers();

        //- Roma et al. three-point kernel
        static scalar phi(const scalar r);


public:

    // Constructors

        //- Construct from the body and the mesh
        immersedMarkers
        (
            const immersedBody& body,
            const fvMesh& mesh,
            const scalar spacing,
            const scalar retraction
        );

        //- No copy construct
        immersedMarkers(const immersedMarkers&) = delete;

        //- No copy assignment
        void operator=(const immersedMarkers&) = delete;


    //- Destructor
    ~immersedMarkers() = default;


    // Member Functions

        // Access

            //- Number of markers
            label size() const
            {
                return faces_.size();
            }

            //- Current forcing points
            const pointField& points() const
            {
                return points_;
            }

            //- Current body velocity at the markers
            const vectorField& velocities() const
            {
                return velocities_;
            }

            //- Marker volumes
            const scalarField& volumes() const
            {
                return volumes_;
            }

            //- Marker areas
            const scalarField& areas() const
            {
                return areas_;
            }

            //- Parent surface face of each marker
            const labelList& faces() const
            {
                return faces_;
            }


        // Edit

            //- Update the markers and their stencils for the current
            //  configuration of the body and the mesh
            void update();


        // Evaluation

            //- Interpolate a cell field to the markers
            tmp<vectorField> interpolate(const UList<vector>& vf) const;

            //- Interpolate a cell field to the markers
            tmp<scalarField> interpolate(const UList<scalar>& sf) const;

            //- Add the spreading of the marker forces F (accelerations) to
            //  the cell field f
            void spread(const UList<vector>& F, UList<vector>& f) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
