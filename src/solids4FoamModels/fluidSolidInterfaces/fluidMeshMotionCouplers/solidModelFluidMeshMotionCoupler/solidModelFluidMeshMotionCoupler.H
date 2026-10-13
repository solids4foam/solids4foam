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
    solidModelFluidMeshMotionCoupler

Description
    Fluid mesh motion coupler for the solids4foamSolidModel mesh motion
    solver, which solves a solids4foam solid model on a static reference copy
    of the fluid mesh, e.g. with the diffusionElastic or diffusionHyperElastic
    mechanical laws.

    On each fluid interface patch, the total displacement from the reference
    points (see displacementFluidMeshMotionCoupler) is set on the solid model
    point displacement, pointD, and, interpolated to the faces, on the
    fixedDisplacement boundary condition of the solid model displacement, D,
    which must not set pointD itself:

    \verbatim
    interface
    {
        type                    fixedDisplacement;
        enforcePointPatchField  no;
        value                   uniform (0 0 0);
    }
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    solidModelFluidMeshMotionCoupler.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAM_COM

#ifndef solidModelFluidMeshMotionCoupler_H
#define solidModelFluidMeshMotionCoupler_H

#include "fluidMeshMotionCoupler.H"
#include "meshMotionSolidModelFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class solidModelFluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class solidModelFluidMeshMotionCoupler
:
    public fluidMeshMotionCoupler
{
    // Private data

        //- Solid model mesh motion solver
        meshMotionSolidModelFvMotionSolver& motionSolver_;

        //- Index of the reference mesh patch for each fluid mesh patch
        labelList referencePatchIDs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidModelFluidMeshMotionCoupler
        (
            const solidModelFluidMeshMotionCoupler&
        );

        //- Disallow default bitwise assignment
        void operator=(const solidModelFluidMeshMotionCoupler&);


public:

    //- Runtime type information
    TypeName("solidModel");


    // Constructors

        //- Construct from components
        solidModelFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~solidModelFluidMeshMotionCoupler()
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

#endif // OPENFOAM_COM

// ************************************************************************* //
