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
    newtonQuasiMonolithicCouplingInterface

Description
    Monolithic fluid-solid interaction solver with Newton-Raphson solution
    algorithm, based on the quasi-monolithic approach described by Jaiman and
    Joshi, 2022, Computational Mechanics of Fluid-Structure Interaction
    Computational Methods for Coupled Fluid-Structure Analysis.

    The approach implemented here is a finite volume version of the method
    described in Algorithm 2 of Chapter 6 in Jaiman and Joshi.

    This "quasi" monolithic approach solves for velocity and pressure in the
    fluid domain and velocity in the solid domain in the same linearised
    monolithic system.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    newtonQuasiMonolithicCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef newtonQuasiMonolithicCouplingInterface_H
#define newtonQuasiMonolithicCouplingInterface_H

#include "fluidSolidInterface.H"
#include "foamPetscSnesHelper.H"
#include "labelPairHashes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class newtonQuasiMonolithicCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class newtonQuasiMonolithicCouplingInterface
:
    public fluidSolidInterface,
    public foamPetscSnesHelper
{
    // Private data

        //- Scale factor for the fluid equations to preserve the condition
        //  number of the monlithic system
        const scalar fluidSystemScaleFactor_;

        //- Scale factor for the solid equations to preserve the condition
        //  number of the monlithic system
        const scalar solidSystemScaleFactor_;

        //- Flag to enable fluid-to-solid coupling
        // const Switch fluidToSolidCoupling_;

        //- Flag to enable mesh-to-fluid coupling
        // const Switch meshToFluidCoupling_;

        //- Flag to enable solid-to-mesh coupling
        // const Switch solidToMeshCoupling_;

        //- Flag to extrapolate the displacement at the solid interface
        //  Alternatively, the adjacenet cell centre displacement is used
        //  This affects the residual calculation, whereas the Jacobian always
        //  uses the adjacenet cell centre displacement
        const Switch extrapolateSolidInterfaceDisplacement_;

        //- Pass the viscous stress from the fluid to the solid
        //  Defaults to yes
        const Switch passViscousStress_;

        //- Use the Liu (2014) Eq.31 interface velocity condition
        //  U_fluid = (3/4)*U_solid + (1/2)*U_solid_old - (1/4)*U_solid_oldold
        //  This condition is from Liu et al., JCP 270, 2014, Section 4.
        //  Note: in the FVM-ALE context, this condition creates a mismatch
        //  between the mesh velocity (U_solid) and the fluid BC (Liu combo),
        //  violating the GCL for interface-adjacent cells. This causes
        //  pressure oscillations in practice. Defaults to no.
        const Switch liuInterfaceCondition_;

        //- Number of sub-regions: 2 == fluid + solid
        const label nRegions_;

        //- Pointer to sub-matrices
        Mat* subMatsPtr_;

        //- PETSc index set for the fluid unknowns
        mutable IS isFluid_;

        //- PETSc index set for the fluid velocity unknowns
        mutable IS isFluidVelocity_;

        //- PETSc index set for the fluid pressure unknowns
        mutable IS isFluidPressure_;

        //- PETSc index set for the motion unknowns
        //mutable IS isMotion_;

        //- PETSc index set for the solid unknowns
        mutable IS isSolid_;

        //- Flag to indicate the fluid index fields have been set
        mutable bool configuredNestedFluidSplit_;

        //- Pointer to the MatNest used for assembly (stored after first
        //  formJacobian call). Not owned — it is A_.m from the base class.
        Mat nestMat_;

        //- Converted monolithic MPIAIJ preconditioner matrix (owned)
        Mat Pmat_;

        //- Optional file for recording the time step sizes and number of
        //  nonlinear iterations
        autoPtr<OFstream> tsLogPtr_;

        // //- Store the old time value for resetting the time
        // scalar oldTimeValue_;

        // //- Number of consecutive failed SNES solves
        // label nConsecutiveFailedSolves_;

        // //- Maximum allowed number of consecutive failed SNES solves
        // const label maxAllowedConsecutiveFailedSolves_;


    // Private Member Functions

        //- Make the isFluid, isFluidVelocity, isFluidPressure, isMotion and
        //  isSolid index sets
        void makeIsFluidIsMotionIsSolid() const;

        //- Return isFluid pointer
        IS isFluid() const
        {
            if (isFluid_ == nullptr)
            {
                makeIsFluidIsMotionIsSolid();
            }

            return isFluid_;
        }

        //- Return isFluidVelocity pointer
        IS isFluidVelocity() const
        {
            if (isFluidVelocity_ == nullptr)
            {
                makeIsFluidIsMotionIsSolid();
            }

            return isFluidVelocity_;
        }

        //- Return isFluidPressure pointer
        IS isFluidPressure() const
        {
            if (isFluidPressure_ == nullptr)
            {
                makeIsFluidIsMotionIsSolid();
            }

            return isFluidPressure_;
        }

        //- Return isSolid pointer
        // IS isMotion() const
        // {
        //     if (isMotion_ == nullptr)
        //     {
        //         makeIsFluidIsMotionIsSolid();
        //     }

        //     return isMotion_;
        // }

        //- Return isSolid pointer
        IS isSolid() const
        {
            if (isSolid_ == nullptr)
            {
                makeIsFluidIsMotionIsSolid();
            }

            return isSolid_;
        }

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

        //- Initialise the Afs Jacobian matrix
        //  This represents the solid velocity coupling in the fluid equations
        label initialiseAfs
        (
            Mat Afs,
            const fvMesh& fluidMesh,
            const label fluidBlockSize,
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
        label formAfs
        (
            Mat Afs,
            const label fluidBlockSize,
            const label solidBlockSize,
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

        //- Map the velocity from the mesh motion to the fluid interface
        void mapInterfaceMotionUToFluidU();

        //- Set the fluid interface velocity using Liu (2014) Eq.31
        //  U_fluid = (3/4)*U_solid + (1/2)*U_solid_old - (1/4)*U_solid_oldold
        //  This replaces mapInterfaceMotionUToFluidU() when
        //  liuInterfaceCondition_ is enabled
        void setLiuInterfaceVelocity();

        //- Map the displacement and velocity from the solid to the mesh motion
        //  interface. Three fields are mapped: D, U and pointD
        void mapInterfaceSolidToMeshMotion();

        //- Reset the fields to the old time fields, e.g. after a failed time
        //  step
        // void resetFieldsToOldTime();

        //- Retrieve UFluid, p and USolid from the PETSc solution vector
        void retrieveSolution
        (
            Vec solution,
            volVectorField& UFluid,
            volScalarField& p,
            volVectorField& USolid,
            const label fluidBlockSize,
            const label solidBlockSize,
            const bool twoD
        ) const;

        //- Disallow default bitwise copy construct
        newtonQuasiMonolithicCouplingInterface
        (
            const newtonQuasiMonolithicCouplingInterface&
        );

        //- Disallow default bitwise assignment
        void operator=(const newtonQuasiMonolithicCouplingInterface&);


public:

    //- Runtime type information
    TypeName("NewtonQuasiMonolithic");


    // Constructors

        //- Construct from components
        newtonQuasiMonolithicCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~newtonQuasiMonolithicCouplingInterface()
        {
            if (subMatsPtr_ != nullptr)
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

            // Free the converted preconditioner matrix
            MatDestroy(&Pmat_);

            // Free index sets
            ISDestroy(&isFluid_);
            ISDestroy(&isFluidVelocity_);
            ISDestroy(&isFluidPressure_);
            // ISDestroy(&isMotion_);
            ISDestroy(&isSolid_);
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

            //- Called after SNESSetFromOptions(), before first SNESSolve()
            virtual void customiseSolver();

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector
            virtual label formResidual
            (
                Vec f,        // Residual
                const Vec x   // Solution
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,     // Jacobian
                const Vec x  // Solution
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
