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
    newtonCouplingInterface

Description
    Monolithic fluid-solid interaction solver with Newton-Raphson solution
    algorithm.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    newtonCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef newtonCouplingInterface_H
#define newtonCouplingInterface_H

#include "fluidSolidInterface.H"
#include "foamPetscSnesHelper.H"
#include "labelPairHashes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class newtonCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class newtonCouplingInterface
:
    public fluidSolidInterface,
    public foamPetscSnesHelper
{
    // Private data

        //- Scale factor for the mesh motion equations to preserve the condition
        //  number of the monlithic system
        //  We default this to be the reciprocal of the motion's shear modulus
        const scalar motionSystemScaleFactor_;

        //- Scale factor for the solid equations to preserve the condition
        //  number of the monlithic system
        //  We default this to be the reciprocal of the solid's shear modulus
        const scalar solidSystemScaleFactor_;

        //- Flag to enable fluid-to-solid coupling
        const Switch fluidToSolidCoupling_;

        //- Flag to enable solid-to-fluid coupling
        const Switch solidToFluidCoupling_;

        //- Flag to enable mesh-to-fluid coupling
        const Switch meshToFluidCoupling_;

        //- Flag to enable solid-to-mesh coupling
        const Switch solidToMeshCoupling_;

        //- Flag to extrapolate the displacement at the solid interface
        //  Alternatively, the adjacenet cell centre displacement is used
        //  This affects the residual calculation, whereas the Jacobian always
        //  uses the adjacenet cell centre displacement
        const Switch extrapolateSolidInterfaceDisplacement_;

        //- Pass the viscous stress from the fluid to the solid
        //  Defaults to yes
        const Switch passViscousStress_;

        //- Flag to indicate if the mesh motion should be included in the
        //  monolithic system; if not, mesh motion is solved sequentially with
        //  the monolithic fluid-solid system
        //  Defaults to "off". Monolithic mesh motion is still work in progress
        const Switch monolithicMeshMotion_;

        //- Number of sub-regions: 2 for fluid+solid, 3 for
        //  fluid + meshMotion + solid
        const label nRegions_;

        //- Pointer to sub-matrices
        Mat* subMatsPtr_;


    // Private Member Functions

        //- Return a reference to the mesh motion solver
        foamPetscSnesHelper& motion();

        //- Return a reference to the mesh motion solver cast as a solid
        //  model
        solidModel& motionSolid();

        //- Create the empty submatrices and nested matrix given the
        //  corresponding solution lengths and block sizes
        //  Note: we pass a reference to the subMats pointer to allow us to
        //  change (e.g. allocate) the pointer
        void createSubMatsAndMat
        (
            Mat& jac,
            Mat*& subMatsPtr,
            const labelPairList& nBlocksAndBlockSize,
            const labelPairHashSet& nullSubMats = labelPairHashSet()
        ) const;

        //- Initialise the Afm Jacobian matrix
        //  This represents the mesh motion coupling in the fluid equations
        label initialiseAfm
        (
            Mat Afm,
            const fvMesh& fluidMesh,
            const label fluidBlockSize,
            const label motionBlockSize,
            const bool twoD
        ) const;

        //- Initialise the Afs Jacobian matrix
        //  This represents the solid coupling in the fluid equations
        label initialiseAfs
        (
            Mat Afs,
            const fvMesh& fluidMesh,
            const label fluidBlockSize,
            const label solidBlockSize,
            const bool twoD
        ) const;

        //- Initialise the Ams Jacobian matrix
        //  This represents the solid coupling in the mesh motion equations
        label initialiseAms
        (
            Mat Ams,
            const fvMesh& fluidMesh, // same as motion mesh
            const label motioBlockSize,
            const label solidBlockSize,
            const bool twoD
        ) const;

        //- Initialise the Asf Jacobian matrix
        //  This represents the fluid coupling in the solid equations
        label initialiseAsf
        (
            Mat Afs,
            const fvMesh& solidMesh,
            const label fluidBlockSize,
            const label solidBlockSize,
            const bool twoD
        ) const;

        //- Insert Afm matrix coefficients
        //  This represents how the mesh motion affects the fluid equations
        //  This does NOT zero the entries first
        label formAfm
        (
            Mat Afm,
            const label fluidBlockSize,
            const label motionBlockSize,
            const bool twoD
        );

        //- Insert Afs matrix coefficients
        //  This represents how the solid affects the fluid equations
        //  This does NOT zero the entries first
        label formAfs
        (
            Mat Afs,
            const label fluidBlockSize,
            const label solidBlockSize,
            const bool twoD
        );

        //- Insert Ams matrix coefficients
        //  This represents how the solid affects the mesh motion equations
        //  This does NOT zero the entries first
        label formAms
        (
            Mat Ams,
            const label solidBlockSize,
            const label motionBlockSize,
            const bool twoD
        );

        //- Insert Asf matrix coefficients
        //  This represents how the fluid affects the solid equations
        //  This does NOT zero the entries first
        label formAsf
        (
            Mat Asf,
            const label fluidBlockSize,
            const label solidBlockSize,
            const bool twoD
        );

        //- Map the velocity from the solid to the fluid interface
        void mapInterfaceSolidUToFluidU();

        //- Map the displacement from the solid to the mesh motion interface
        //  This maps the current solid displacement interface field to the
        //  mesh motion displacement field. Both D and pointD are mapped
        void mapInterfaceSolidUToMeshMotionD();

        //- Predict the solution field at the start of a new time step
        void predict();

        //- Disallow default bitwise copy construct
        newtonCouplingInterface(const newtonCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const newtonCouplingInterface&);


public:

    //- Runtime type information
    TypeName("NewtonMonolithic");


    // Constructors

        //- Construct from components
        newtonCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~newtonCouplingInterface()
        {
            // Destroy each submatrix
            for (label i = 0; i < nRegions_*nRegions_; ++i)
            {
                MatDestroy(&subMatsPtr_[i]);
            }

            // Free the array of Mat handles
            delete[] subMatsPtr_;
            subMatsPtr_ = nullptr;
        }


    // Member Functions

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the interface
            virtual bool evolve();

            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            virtual label initialiseJacobian(Mat& jac);

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x);

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector
            virtual label formResidual
            (
                PetscScalar *f,        // Residual
                const PetscScalar *x   // Solution
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,              // Jacobian
                const PetscScalar *x  // Solution
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
