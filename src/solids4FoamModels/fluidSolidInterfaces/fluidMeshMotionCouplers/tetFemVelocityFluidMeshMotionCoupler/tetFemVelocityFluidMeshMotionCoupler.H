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
    tetFemVelocityFluidMeshMotionCoupler

Description
    Fluid mesh motion coupler for the foam-extend tetrahedral finite element
    mesh motion solvers, which solve for the point velocity field motionU.

    The motionU value on each fluid interface patch is set to the change in
    the patch point displacement divided by the time step.

Author
    Philip Cardiff, UCD.

SourceFiles
    tetFemVelocityFluidMeshMotionCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef OPENFOAM_NOT_EXTEND

#ifndef tetFemVelocityFluidMeshMotionCoupler_H
#define tetFemVelocityFluidMeshMotionCoupler_H

#include "fluidMeshMotionCoupler.H"
#include "tetPointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class tetFemVelocityFluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class tetFemVelocityFluidMeshMotionCoupler
:
    public fluidMeshMotionCoupler
{
    // Private data

        //- Point velocity field of the mesh motion solver
        tetPointVectorField& motionU_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetFemVelocityFluidMeshMotionCoupler
        (
            const tetFemVelocityFluidMeshMotionCoupler&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetFemVelocityFluidMeshMotionCoupler&);


public:

    //- Runtime type information
    TypeName("tetFemVelocity");


    // Constructors

        //- Construct from components
        tetFemVelocityFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~tetFemVelocityFluidMeshMotionCoupler()
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

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
