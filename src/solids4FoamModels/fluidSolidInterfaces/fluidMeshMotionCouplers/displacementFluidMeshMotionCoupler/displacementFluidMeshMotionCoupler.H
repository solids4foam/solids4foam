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
    displacementFluidMeshMotionCoupler

Description
    Fluid mesh motion coupler for displacement-based mesh motion solvers,
    e.g. displacementLaplacian and displacementSBRStress, which solve for the
    point displacement field pointDisplacement relative to the reference
    points, points0, of the solver.

    The pointDisplacement value on each fluid interface patch is set to the
    total displacement of the patch points from their reference points:

        pointDisplacement = points + iterDisplacement - points0

    where points are the current patch points and iterDisplacement is the
    change in their displacement since the previous mesh update. The value
    therefore depends only on the current geometry: it is the same whatever
    the number of outer iterations, it includes any direct moves of the
    interface points (interfaceDeformationLimit), and it does not depend on
    the history of the boundary values, e.g. on restart.

    Unlike a velocity-based solver, the mesh for a given interface position
    does not depend on how the interface got there, so the mesh does not
    drift under periodic interface motion.

Author
    Philip Cardiff, UCD.

SourceFiles
    displacementFluidMeshMotionCoupler.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_NOT_EXTEND

#ifndef displacementFluidMeshMotionCoupler_H
#define displacementFluidMeshMotionCoupler_H

#include "fluidMeshMotionCoupler.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class displacementFluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class displacementFluidMeshMotionCoupler
:
    public fluidMeshMotionCoupler
{
    // Private data

        //- Displacement-based mesh motion solver
        displacementMotionSolver& motionSolver_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        displacementFluidMeshMotionCoupler
        (
            const displacementFluidMeshMotionCoupler&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementFluidMeshMotionCoupler&);


public:

    //- Runtime type information
    TypeName("displacement");


    // Constructors

        //- Construct from components
        displacementFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~displacementFluidMeshMotionCoupler()
        {}


    // Member Functions

        //- Return the total displacement of the points of the given patch
        //  from the reference points, points0, after the next mesh update
        static tmp<vectorField> totalPatchDisplacement
        (
            const polyPatch& patch,
            const pointField& points0,
            const vectorField& iterDisplacement
        );

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
