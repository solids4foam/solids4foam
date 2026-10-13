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
    subsetMeshFluidMeshMotionCoupler

Description
    Fluid mesh motion coupler for the foam-extend newSubsetMotionSolverFvMesh,
    whose velocity-based mesh motion solver moves a sub-mesh of the fluid
    mesh.

    The pointMotionU field of the sub-mesh is set as by the velocity coupler;
    the fluid interface patches have the same indices in the sub-mesh.

Author
    Philip Cardiff, UCD.

SourceFiles
    subsetMeshFluidMeshMotionCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef OPENFOAM_NOT_EXTEND

#ifndef subsetMeshFluidMeshMotionCoupler_H
#define subsetMeshFluidMeshMotionCoupler_H

#include "velocityFluidMeshMotionCoupler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class subsetMeshFluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class subsetMeshFluidMeshMotionCoupler
:
    public velocityFluidMeshMotionCoupler
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        subsetMeshFluidMeshMotionCoupler
        (
            const subsetMeshFluidMeshMotionCoupler&
        );

        //- Disallow default bitwise assignment
        void operator=(const subsetMeshFluidMeshMotionCoupler&);


public:

    //- Runtime type information
    TypeName("subsetMesh");


    // Constructors

        //- Construct from components
        subsetMeshFluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~subsetMeshFluidMeshMotionCoupler()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
