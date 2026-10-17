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
    Foam::fsiImmersedInterface

Description
    An immersed fluid-solid interface: the solid interface patch drives an
    immersed body of the fluid (see fsiImmersedBoundary), and there is no
    fluid patch.

    The surface of the body is built from the solid interface zone (the
    global patch of the solid interface patch) and, optionally, closure
    zones (other patches of the solid, e.g. a clamped root), whose union
    must be a closed surface with the normals pointing out of the solid;
    the points that the zones share are merged within a small tolerance.
    The faces are triangulated, and each triangle keeps its parent interface
    face. In a two-dimensional fluid mesh (one empty direction), the
    patches of the empty direction are not needed: the points on either
    side of the mid-plane of the fluid mesh are moved beyond its bounds in
    the empty direction, and the boundary loops of the surface are capped
    by triangulated polygons, so that the body extends through the mesh.

    The points of the surface are the interface zone points, which follow
    the interface displacement of the coupling loop (relaxed or accelerated
    by the coupling scheme), and the closure zone points, which stay at
    their positions at the start of the time step; a fixed offset in the
    empty direction is added in two dimensions. The point velocities of the
    interface are the displacement increment over the time step.

    The traction on each interface face is the area-weighted mean of the
    tractions on its triangles from the immersed boundary, where the areas
    are those represented by the quadrature points of the triangles, so
    that the force on a face partly outside the mesh is that of its part in
    the mesh; a face without quadrature points has zero traction.

Author
    Philip Cardiff, UCD.

SourceFiles
    fsiImmersedInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fsiImmersedInterface_H
#define fsiImmersedInterface_H

#include "fsiImmersedBoundary.H"
#include "standAlonePatch.H"
#include "PtrList.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fsiImmersedInterface Declaration
\*---------------------------------------------------------------------------*/

class fsiImmersedInterface
{
    // Private Data

        //- Immersed boundary option of the fluid
        fsiImmersedBoundary& ib_;

        //- Index of the body in the option
        const label bodyi_;

        //- Name of the body
        const word bodyName_;

        //- Zone of each surface point: -1 for the interface zone, or the
        //  index of the closure zone
        labelList pointZones_;

        //- Index of each surface point in its zone
        labelList pointIndices_;

        //- Offset of each surface point in the empty direction of a
        //  two-dimensional fluid mesh
        vectorField pointOffsets_;

        //- Parent interface zone face of each triangle of the surface, or
        //  -1 for the triangles of the closure zones and the caps
        labelList triangleFaces_;

        //- Number of interface zone faces
        const label nZoneFaces_;

        //- Interface zone points at the start of the time step
        pointField zonePoints0_;

        //- Closure zone points at the start of the time step
        List<pointField> closurePoints0_;

        //- Number of interface faces without quadrature points at the last
        //  traction evaluation
        mutable label nFacesWithoutPoints_;


    // Private Member Functions

        //- Build the surface from the zones and give it to the body
        void buildSurface
        (
            const standAlonePatch& zone,
            const pointField& zonePoints,
            const List<const standAlonePatch*>& closureZones,
            const List<pointField>& closurePoints,
            const Vector<label>& fluidSolutionD,
            const boundBox& fluidBounds
        );

        //- No copy construct
        fsiImmersedInterface(const fsiImmersedInterface&) = delete;

        //- No copy assignment
        void operator=(const fsiImmersedInterface&) = delete;


public:

    // Constructors

        //- Construct from the option, the body, the interface zone and the
        //  closure zones with their reference points, and the solution
        //  directions and bounds of the fluid mesh
        fsiImmersedInterface
        (
            fsiImmersedBoundary& ib,
            const label bodyi,
            const word& bodyName,
            const standAlonePatch& zone,
            const pointField& zonePoints,
            const List<const standAlonePatch*>& closureZones,
            const List<pointField>& closurePoints,
            const Vector<label>& fluidSolutionD,
            const boundBox& fluidBounds
        );


    //- Destructor
    ~fsiImmersedInterface() = default;


    // Member Functions

        //- Name of the body
        const word& bodyName() const
        {
            return bodyName_;
        }

        //- Interface zone points at the start of the time step
        const pointField& zonePoints0() const
        {
            return zonePoints0_;
        }

        //- Number of interface faces without quadrature points at the last
        //  traction evaluation
        label nFacesWithoutPoints() const
        {
            return nFacesWithoutPoints_;
        }

        //- Set the interface and closure zone points at the start of the
        //  time step
        void setPoints0
        (
            const pointField& zonePoints,
            const List<pointField>& closurePoints
        );

        //- Move the surface: the interface zone points are displaced from
        //  their positions at the start of the time step by the given
        //  displacement, with the given velocity
        void move
        (
            const vectorField& zoneDisplacement,
            const vectorField& zoneVelocity
        );

        //- Traction exerted by the fluid on each interface zone face, from
        //  the latest fluid solution
        tmp<vectorField> zoneTraction() const;

        //- Force on the whole body from the momentum exchange with the
        //  fluid, from the latest fluid solution
        vector momentumExchangeForce() const
        {
            return ib_.momentumExchangeForce(bodyi_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
