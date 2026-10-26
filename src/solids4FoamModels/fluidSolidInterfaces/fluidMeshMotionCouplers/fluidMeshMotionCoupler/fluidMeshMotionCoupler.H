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
    fluidMeshMotionCoupler

Description
    Virtual base class for passing the motion of the fluid-solid interfaces
    to the mesh motion solver of the fluid mesh.

    The fluid-solid interface prescribes, for each fluid interface patch, the
    change in the displacement of the patch points since the previous mesh
    update; a derived class sets this motion in the form its family of mesh
    motion solvers expects, e.g. as a velocity boundary condition on
    pointMotionU, or as a total displacement on pointDisplacement.

    The coupler is selected with the optional fluidMeshMotion sub-dictionary
    in fsiProperties:

    \verbatim
    fluidMeshMotion
    {
        type    auto;
    }
    \endverbatim

    where type auto, the default, selects the coupler from the mesh motion
    solver of the fluid mesh.

Author
    Philip Cardiff, UCD.

SourceFiles
    fluidMeshMotionCoupler.C
    fluidMeshMotionCouplerNew.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMeshMotionCoupler_H
#define fluidMeshMotionCoupler_H

#include "dynamicFvMesh.H"
#include "motionSolver.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidMeshMotionCoupler Declaration
\*---------------------------------------------------------------------------*/

class fluidMeshMotionCoupler
{
    // Private data

        //- Fluid mesh
        dynamicFvMesh& mesh_;

        //- Indices of the fluid interface patches
        const labelList patchIDs_;

        //- Coupler dictionary
        const dictionary dict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluidMeshMotionCoupler(const fluidMeshMotionCoupler&);

        //- Disallow default bitwise assignment
        void operator=(const fluidMeshMotionCoupler&);


public:

    //- Runtime type information
    TypeName("fluidMeshMotionCoupler");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidMeshMotionCoupler,
            dictionary,
            (
                dynamicFvMesh& mesh,
                const labelList& patchIDs,
                const dictionary& dict
            ),
            (mesh, patchIDs, dict)
        );


    // Constructors

        //- Construct from components
        fluidMeshMotionCoupler
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Selectors

        //- Select the coupler given by the type entry in dict, where type
        //  auto, the default, selects it from the mesh motion solver
        static autoPtr<fluidMeshMotionCoupler> New
        (
            dynamicFvMesh& mesh,
            const labelList& patchIDs,
            const dictionary& dict
        );


    // Destructor

        virtual ~fluidMeshMotionCoupler()
        {}


    // Member Functions

        // Access

            //- Return the fluid mesh
            dynamicFvMesh& mesh()
            {
                return mesh_;
            }

            //- Return the fluid mesh
            const dynamicFvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return the indices of the fluid interface patches
            const labelList& patchIDs() const
            {
                return patchIDs_;
            }

            //- Return the coupler dictionary
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return a pointer to the mesh motion solver of the mesh, or
            //  nullptr if it does not have one that can be found
            static motionSolver* lookupMotionSolver(const dynamicFvMesh& mesh);

            //- Return the type of coupler for the mesh motion solver of the
            //  mesh, or a FatalError if there is none
            static word autoSelectType(const dynamicFvMesh& mesh);

        // Edit

            //- Set the motion of the given fluid interface patch for the next
            //  mesh update, where iterDisplacement is the change in the
            //  displacement of the patch points since the previous mesh
            //  update, in the order of the local patch points
            virtual void setPatchMotion
            (
                const label patchID,
                const vectorField& iterDisplacement
            ) = 0;

            //- Complete the motion once it has been set on all the fluid
            //  interface patches, before the mesh is updated
            virtual void commit()
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
