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

        //- Pass the viscous stress from the fluid to the solid
        //  Defaults to yes
        const Switch passViscousStress_;

        //- Transfer the fluid interface forces to the solid with the
        //  transpose of the solid-to-fluid velocity transfer, instead of
        //  interpolating the traction, so that the total interface force and
        //  power are conserved on a non-conformal interface. The solid
        //  applies the traction on the undeformed face areas. Defaults to no
        const Switch conservativeTractionTransfer_;

        //- Report interface normal flux defect diagnostics
        const Switch interfaceFluxDiagnostics_;

        //- Offset subtracted from the mapped solid interface displacement on
        //  the local fluid patch faces, so that the fluid mesh follows the
        //  solid displacement relative to the start of the coupling
        autoPtr<vectorField> interfaceDOffsetPtr_;

        //- Offset subtracted from the mapped solid interface point
        //  displacement on the local fluid patch points
        autoPtr<vectorField> interfacePointDOffsetPtr_;

        //- Solid zone faces giving the value on each fluid zone face
        mutable autoPtr<labelListList> fluidFaceSolidFacesPtr_;

        //- Weights of the solid zone faces giving each fluid zone face value
        mutable autoPtr<scalarListList> fluidFaceSolidWeightsPtr_;

        //- Fluid zone faces giving the value on each solid zone face
        mutable autoPtr<labelListList> solidFaceFluidFacesPtr_;

        //- Weights of the fluid zone faces giving each solid zone face value
        mutable autoPtr<scalarListList> solidFaceFluidWeightsPtr_;

        //- Accumulated signed leakage from interface flux defect
        scalar accumulatedInterfaceFluxLeakage_;

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

        //- Calculate the interface face transfer addressing and weights in
        //  both directions from the interface-to-interface mapping. For a
        //  conformal interface (directMap), each face has one unit weight;
        //  for a non-conformal interface (AMI), the weights are the face
        //  overlap areas divided by the area of the receiving face
        void calcInterfaceFaceWeights() const;

        //- Return the solid zone faces giving each fluid zone face value
        const labelListList& fluidFaceSolidFaces() const;

        //- Return the weights of the solid zone faces for each fluid face
        const scalarListList& fluidFaceSolidWeights() const;

        //- Return the fluid zone faces giving each solid zone face value
        const labelListList& solidFaceFluidFaces() const;

        //- Return the weights of the fluid zone faces for each solid face
        const scalarListList& solidFaceFluidWeights() const;

        //- Transfer a solid zone face field to the fluid zone faces
        template<class Type>
        tmp<Field<Type>> solidZoneToFluidZone
        (
            const Field<Type>& solidZoneField
        ) const;

        //- Transfer a fluid zone face field to the solid zone faces
        template<class Type>
        tmp<Field<Type>> fluidZoneToSolidZone
        (
            const Field<Type>& fluidZoneField
        ) const;

        //- Apply the transpose of the solid-to-fluid zone face transfer to a
        //  fluid zone face field, e.g. to distribute face forces
        template<class Type>
        tmp<Field<Type>> fluidZoneToSolidZoneTranspose
        (
            const Field<Type>& fluidZoneField
        ) const;

        //- Return the local solid interface face areas on which the solid
        //  applies the interface traction
        tmp<scalarField> solidInterfaceMagSf() const;

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

        //- Set the normal gradient of the fluid pressure on the interface
        //  from the normal momentum equation, if the pressure condition
        //  there is fixedGradient (see
        //  newtonIcoFluid::updatePressureGradientBoundaries). A zeroGradient
        //  pressure instead gives a first-order error in the interface
        //  pressure and the added mass
        void setInterfacePressureGradient();

        //- Return current geometric interface mesh flux
        scalarField currentInterfaceMeshPhi() const;

        //- Return current fluid-boundary normal velocity flux
        scalarField currentInterfaceVelocityPhi() const;

        //- Report interface pressure flux defect diagnostics: the defect
        //  between the interface velocity flux and the mesh flux
        void reportInterfaceFluxDiagnostics();

        //- Map the displacement and velocity from the solid to the mesh motion
        //  interface. Three fields are mapped: D, U and pointD
        void mapInterfaceSolidToMeshMotion();

        //- Set the interface displacement offsets from the difference between
        //  the mapped solid interface displacement and the mesh motion
        //  interface displacement at the start of the coupling, which is
        //  zero for an undeformed start and recovers the offset on a restart
        void setInterfaceDisplacementOffsets();

        //- Restore the coupled fields and fluid mesh after a failed PETSc solve
        void restoreOldTimeState(const pointField& oldFluidPoints);

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

        //- Seed the PETSc solution vector from the current UFluid, p
        //  and USolid fields. Applies the inverse fluid pressure
        //  unknown scale (pHat = p/pressureUnknownScale) so the Vec
        //  stores the SNES-side unknowns consistently
        void seedSolution
        (
            Vec solution,
            const volVectorField& UFluid,
            const volScalarField& p,
            const volVectorField& USolid,
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
            ISDestroy(&isSolid_);
        }


    // Member Functions

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- The solid velocity is imposed on the fluid interface as part
            //  of the monolithic system
            virtual bool imposesFluidInterfaceVelocity() const
            {
                return true;
            }

            //- Evolve the interface
            virtual bool evolve();

            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            virtual label initialiseJacobian(Mat& jac);

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x);

            //- Called after the PETSc options are pushed, before SNESSolve()
            virtual void customiseSolver();

            //- Called after resetting the SNES/KSP/PC setup
            virtual void resetCustomSolverState();

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
