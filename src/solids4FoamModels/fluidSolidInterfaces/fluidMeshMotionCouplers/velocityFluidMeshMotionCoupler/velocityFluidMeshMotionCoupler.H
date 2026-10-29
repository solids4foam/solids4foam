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
    velocityFluidMeshMotionCoupler

Description
    Fluid mesh motion coupler for velocity-based mesh motion solvers, e.g.
    velocityLaplacian, which solve for the point velocity field pointMotionU.

    The pointMotionU value on each fluid interface patch is set to the change
    in the patch point displacement divided by the time step.

Author
    Philip Cardiff, UCD.

SourceFiles
    velocityFluidMeshMotionCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef velocityFluidMeshMotionCoupler_H
#define velocityFluidMeshMotionCoupler_H

#include "fluidMeshMotionCoupler.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class velocityFluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class velocityFluidMeshMotionCoupler
:
    public fluidMeshMotionCoupler
{
    // Private data

        //- Point velocity field of the mesh motion solver
        pointVectorField& pointMotionU_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        velocityFluidMeshMotionCoupler(const velocityFluidMeshMotionCoupler&);

        //- Disallow default bitwise assignment
        void operator=(const velocityFluidMeshMotionCoupler&);


protected:

    // Protected Constructors

        //- Construct from components, where the pointMotionU field is
        //  looked up on motionMesh, e.g. a sub-mesh of the fluid mesh
        velocityFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const fvMesh& motionMesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


public:

    //- Runtime type information
    TypeName("velocity");


    // Constructors

        //- Construct from components
        velocityFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~velocityFluidMeshMotionCoupler()
        {}


    // Member Functions

        //- Set the motion of the given fluid interface patch
        virtual void setPatchMotion
        (
            const label patchID,
            const vectorField& iterDisplacement
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
