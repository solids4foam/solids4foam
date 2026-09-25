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
    Foam::immersedBody

Description
    A body immersed in the fluid mesh, given by a closed surface (e.g. an STL
    file) and a prescribed motion.

    The body gives the solid volume fraction (occupancy) lambda of the cells
    it covers:

    - signedDistance (the default): lambda = 0.5 + d/w, clipped to [0, 1],
      where d is the distance of the cell centre to the surface, positive
      inside the body, and w is the width of the cell normal to the surface
      (the width of its bounding box in the direction of the normal of the
      nearest surface face). lambda varies continuously as the body moves.
    - vertexFraction: half the fraction of the cell vertices inside the
      surface, plus a half if the cell centre is inside the surface, as in
      the non-convex body of openHFDIB-DEM.

    Cells with an occupancy greater than 1 - threshold are internal cells, and
    cells with an occupancy greater than threshold are surface cells.

    The surface is replicated on every processor, so every processor moves it
    identically and locates the cells of its own part of the mesh.

    Usage, as an entry of the bodies dictionary of immersedBoundaryForce:
    \verbatim
    cylinder
    {
        // Surface file, relative to constant/triSurface
        surface     "cylinder.stl";

        // Prescribed motion, optional, default static
        motion
        {
            type        sinusoidalTranslation;
            amplitude   0.25;
            period      4;
        }

        // Centre of rotation for the torque, optional, default the mean of
        // the reference surface points
        CofR        (0 0 0);

        // Occupancy, optional, default signedDistance
        occupancy   signedDistance;
    }
    \endverbatim

Author
    Derived from openHFDIB-DEM (https://github.com/techMathGroup/openHFDIB-DEM,
    commit b4ea321), itself derived from openHFDIB
    (https://github.com/fmuni/openHFDIB). openHFDIB-DEM contributors:
    Federico Municchi (2016), Martin Isoz (2019-*), Martin Kotouč Šourek
    (2019-*), Ondřej Studeník (2020-*). openHFDIB and openHFDIB-DEM are
    distributed under the GNU Lesser General Public License version 3, whose
    terms allow this derived code to be distributed under the GNU General
    Public License version 3.

    Rewritten for prescribed motion in solids4foam by Philip Cardiff, UCD,
    from the version adapted by Sairam Pamulaparthi Venkata, UCD, in
    cardiacFoam.

SourceFiles
    immersedBody.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBody_H
#define immersedBody_H

#include "fvMesh.H"
#include "volFields.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class immersedBody Declaration
\*---------------------------------------------------------------------------*/

class immersedBody
{
    // Private Data

        //- Name of the body
        const word name_;

        //- Reference to the fluid mesh
        const fvMesh& mesh_;

        //- Surface in the current configuration
        triSurface surface_;

        //- Surface points in the reference configuration
        const pointField points0_;

        //- Search engine for the current surface
        autoPtr<triSurfaceSearch> searchPtr_;

        //- Prescribed motion
        autoPtr<immersedBodyMotion> motionPtr_;

        //- Centre of rotation for the torque, in the reference configuration
        const point CofR0_;

        //- Occupancy from the signed distance of the cell centres to the
        //  surface, rather than from the fraction of cell vertices inside
        const bool signedDistance_;

        //- Time of the current configuration
        scalar time_;

        //- Local cells fully covered by the body
        labelList internalCells_;

        //- Local cells partially covered by the body
        labelList surfaceCells_;


    // Private Member Functions

        //- Read the surface
        static triSurface readSurface(const fvMesh& mesh, const fileName& f);


public:

    // Constructors

        //- Construct from name, dictionary and mesh
        immersedBody
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        immersedBody(const immersedBody&) = delete;

        //- No copy assignment
        void operator=(const immersedBody&) = delete;


    //- Destructor
    ~immersedBody() = default;


    // Member Functions

        // Access

            //- Name of the body
            const word& name() const
            {
                return name_;
            }

            //- Does the body move?
            bool moving() const
            {
                return motionPtr_->moving();
            }

            //- Surface in the current configuration
            const triSurface& surface() const
            {
                return surface_;
            }

            //- Local cells fully covered by the body
            const labelList& internalCells() const
            {
                return internalCells_;
            }

            //- Local cells partially covered by the body
            const labelList& surfaceCells() const
            {
                return surfaceCells_;
            }

            //- Centre of rotation for the torque, in the current
            //  configuration
            point CofR() const;


        // Edit

            //- Move the surface to its configuration at time t
            void move(const scalar t);

            //- Locate the covered cells and add the occupancy of the body to
            //  lambda, which is clipped to [0, 1]
            void addOccupancy
            (
                volScalarField& lambda,
                const scalar threshold
            );

            //- Set the velocity of the body in the covered cells of Ui
            void setVelocity(volVectorField& Ui) const;


        // Evaluation

            //- Force on the fluid by the body, given by the immersed boundary
            //  forcing f (an acceleration): -rho*sum(f*V) over the covered
            //  cells. This is also the force on the body by the fluid,
            //  excluding the inertia of the fluid inside the body
            vector force(const volVectorField& f, const scalar rho) const;

            //- Torque about CofR corresponding to force
            vector torque(const volVectorField& f, const scalar rho) const;

            //- Momentum of the fluid inside the body, per unit density:
            //  sum(lambda*U*V) over the covered cells
            vector momentum
            (
                const volVectorField& U,
                const volScalarField& lambda
            ) const;


        // Write

            //- Write the current surface to the given file
            void writeSurface(const fileName& file) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
