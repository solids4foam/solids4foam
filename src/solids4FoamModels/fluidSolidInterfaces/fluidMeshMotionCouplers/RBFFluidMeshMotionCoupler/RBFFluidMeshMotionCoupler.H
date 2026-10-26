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
    RBFFluidMeshMotionCoupler

Description
    Fluid mesh motion coupler for the radial basis function mesh motion
    solver, RBFMeshMotionSolver.

    The change in the displacement of the fluid interface patch points is
    interpolated to the patch faces and passed, with zero motion on all
    other patches, to the solver.

Author
    Philip Cardiff, UCD.

SourceFiles
    RBFFluidMeshMotionCoupler.C

\*---------------------------------------------------------------------------*/

#if !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)

#ifndef RBFFluidMeshMotionCoupler_H
#define RBFFluidMeshMotionCoupler_H

#include "fluidMeshMotionCoupler.H"
#include "RBFMeshMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class RBFFluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class RBFFluidMeshMotionCoupler
:
    public fluidMeshMotionCoupler
{
    // Private data

        //- RBF mesh motion solver
        RBFMeshMotionSolver& motionSolver_;

        //- Face motion of every patch of the fluid mesh, for the next update
        Field<vectorField> motion_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RBFFluidMeshMotionCoupler(const RBFFluidMeshMotionCoupler&);

        //- Disallow default bitwise assignment
        void operator=(const RBFFluidMeshMotionCoupler&);


public:

    //- Runtime type information
    TypeName("RBF");


    // Constructors

        //- Construct from components
        RBFFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~RBFFluidMeshMotionCoupler()
        {}


    // Member Functions

        //- Set the motion of the given fluid interface patch
        virtual void setPatchMotion
        (
            const label patchID,
            const vectorField& iterDisplacement
        );

        //- Pass the motion of all patches to the RBF solver
        virtual void commit();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)

// ************************************************************************* //
