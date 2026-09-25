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
    the surface.

    At each point \f$ x_w \f$, with the outward unit normal n, the velocity
    and pressure and their gradients are fitted by weighted linear least
    squares at the image point \f$ x_I = x_w + d n \f$ (d = 1.5 h), over the
    cells within 2 h of it. The weights are a smooth function of the distance
    to the image point, times a smooth function of the distance \f$ \psi \f$
    of the cell centre from the tangent plane of the surface, which is zero
    for the cells within h/4 of the surface or inside the body, so that the
    fit uses fluid cells only and varies continuously as the body moves.

    The wall pressure is extrapolated with the fitted pressure gradient, and
    the tangential viscous traction is

    \f[
        p_w = p_I - d \, \nabla p_I \cdot n
    \f]

    \f[
        \tau_w = \nu \left[ \frac{(I - n n)(U_I - U_b)}{d}
      - \frac{d}{2} \frac{(I - n n)(a_b + \nabla p_I)}{\nu} \right]
    \f]

    from the tangential momentum balance at the wall, where \f$ U_b \f$
    and \f$ a_b \f$ are the velocity and acceleration of the body and p is
    the kinematic pressure. The second term makes the tangential traction
    second order accurate. The traction on the body is
    \f$ \rho (-p_w n + \tau_w) \f$.

    The sums of the fits over the processors are reduced, so that the cells
    of a fit may be on any processor.

Author
    Philip Cardiff, UCD.

SourceFiles
    immersedSurfaceTraction.C

\*---------------------------------------------------------------------------*/

#ifndef immersedSurfaceTraction_H
#define immersedSurfaceTraction_H

#include "immersedBody.H"
#include "FixedList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class immersedSurfaceTraction Declaration
\*---------------------------------------------------------------------------*/

class immersedSurfaceTraction
{
    // Private Data

        //- Reference to the fluid mesh
        const fvMesh& mesh_;

        //- The body
        const immersedBody& body_;

        //- Solved directions
        const Vector<label> solD_;

        //- Parent surface face of each point
        labelList faces_;

        //- Barycentric coordinates of each point in its face
        List<FixedList<scalar, 3>> barycentric_;

        //- Surface area represented by each point (including the mesh
        //  thickness in two dimensions)
        scalarField areas_;

        //- Local cell width at each point
        scalarField h_;

        //- Current positions and outward unit normals
        pointField points_;
        vectorField normals_;


    // Private Member Functions

        //- Width of a cell: its smallest bounding box width in the solved
        //  directions
        scalar cellWidth(const label celli) const;

        //- Width of the cell containing each point, the largest over the
        //  processors, and 0 if no processor contains it
        scalarField cellWidths(const pointField& pts) const;

        //- Create the quadrature points
        void createPoints();

        //- Update the positions and normals of the points
        void updatePoints();


public:

    // Constructors

        //- Construct from the body and the mesh
        immersedSurfaceTraction
        (
            const immersedBody& body,
            const fvMesh& mesh
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
        //  surface at each quadrature point
        tmp<vectorField> traction
        (
            const volVectorField& U,
            const volScalarField& p,
            const scalar nu
        );

        //- Force and torque (about CofR) exerted by the fluid on the body
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
