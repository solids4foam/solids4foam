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
    Foam::fsiImmersedBoundary

Description
    Abstract interface of an immersed boundary option of the fluid model
    whose bodies may be driven by a fluid-solid interaction: the fluid-solid
    interface gives a body its surface, moves the surface in each coupling
    iteration, and takes the traction exerted by the fluid on the surface.

    The interface is header-only, so that the immersed boundary library,
    which implements it (fv::immersedBoundaryForce), does not depend on
    solids4FoamModels. The fluid-solid interface finds the option among the
    finite volume options of the fluid mesh by a dynamic cast to this class.

    The surface is a closed triangulated surface with outward normals,
    replicated on every processor; its points and their velocities are given
    in the same order on every processor.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef fsiImmersedBoundary_H
#define fsiImmersedBoundary_H

#include "pointField.H"
#include "vectorField.H"
#include "scalarField.H"
#include "faceList.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fsiImmersedBoundary Declaration
\*---------------------------------------------------------------------------*/

class fsiImmersedBoundary
{
public:

    // Constructors

        //- Null constructor
        fsiImmersedBoundary() = default;


    //- Destructor
    virtual ~fsiImmersedBoundary() = default;


    // Member Functions

        //- Index of the body with the given name, or -1 if there is none
        virtual label findBody(const word& bodyName) const = 0;

        //- Can the surface of the body be driven by the fluid-solid
        //  interface (its motion is fsiDriven)?
        virtual bool fsiDriven(const label bodyi) const = 0;

        //- Set the reference surface of the body: the points, and the
        //  triangles as faces of three points, with the normals pointing
        //  out of the body
        virtual void setReferenceSurface
        (
            const label bodyi,
            const pointField& points,
            const faceList& triangles
        ) = 0;

        //- Move the surface of the body to the given point positions, with
        //  the given point velocities. The bodies are re-positioned in the
        //  fluid mesh at the next assembly of the momentum equation
        virtual void moveSurface
        (
            const label bodyi,
            const pointField& points,
            const vectorField& velocities
        ) = 0;

        //- Traction (force per unit area, including the fluid density)
        //  exerted by the fluid on each triangle of the surface of the body,
        //  from the latest fluid solution, and the area of each triangle
        //  represented by its quadrature points in the mesh (zero for a
        //  triangle without quadrature points, whose traction is zero)
        virtual tmp<vectorField> surfaceTraction
        (
            const label bodyi,
            scalarField& areas
        ) = 0;

        //- Force on the body from the momentum exchange with the fluid
        //  (including the fluid density), from the latest fluid solution,
        //  for checking the conservation of the surface traction
        virtual vector momentumExchangeForce(const label bodyi) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
