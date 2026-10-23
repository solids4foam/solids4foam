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
    newtonIcoFluid

Description
    Incompressible Newtonian laminar fluid model solved used a Jacobian-free
    Newton-Krylov approach via a PETSc SNES interface.

Author
    Philip Cardiff, UCD.

SourceFiles
    newtonIcoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef newtonIcoFluid_H
#define newtonIcoFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#ifdef OPENFOAM_COM
    #include "turbulentTransportModel.H"
    #include "CorrectPhi.H"
#elif defined(OPENFOAM_ORG)
    #include "singlePhaseTransportModel.H"
    #include "kinematicMomentumTransportModel.H"
#else
    #include "turbulenceModel.H"
#endif
#include "foamPetscSnesHelper.H"
#include "solidModel.H"
#include "stabilisationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------* \
                        Class newtonIcoFluid Declaration
\*---------------------------------------------------------------------------*/

class newtonIcoFluid
:
    public fluidModel,
    public foamPetscSnesHelper
{
    // Private data

        //- Pointer to the velocity surface field
        autoPtr<surfaceVectorField> Uf_;

        //- rAU field
        //  Smoothing coefficient in the pressure equation
        mutable autoPtr<surfaceScalarField> rAUfPtr_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
#ifdef OPENFOAM_ORG
        autoPtr<incompressible::momentumTransportModel> turbulence_;
#else
        autoPtr<incompressible::turbulenceModel> turbulence_;
#endif

        //- Fluid density [kg/m3]
        dimensionedScalar rho_;

        //- Momentum stabilisation model
        autoPtr<stabilisationModel> momentumStabilisationPtr_;

        //- Pressure stabilisation model
        autoPtr<stabilisationModel> pressureStabilisationPtr_;

        //- Scale factor applied to the pressure-equation residual and
        //  Jacobian rows relative to the momentum equations
        const scalar pressureScaleFactor_;

        //- Store pressure in the PETSc solution vector as p/pScale.
        //  When disabled (default) the PETSc Vec stores physical
        //  pressure unchanged
        const Switch scaleMixedPetScFields_;

        //- Scale type used for the PETSc pressure unknown.
        //  Valid options: "none", "user"/"scalar", "dynamicHead"
        const word pressureUnknownScaleType_;

        //- Pressure scale used for pHat = p/pressureUnknownScale_.
        //  Defaults to 1 (no-op); intended to give matrix-free SNES
        //  (-snes_mf_operator) comparable velocity and pressure
        //  magnitudes in the PETSc Vec
        scalar pressureUnknownScale_;

        //- Effective scale applied to the pressure-equation residual
        //  and the pressure-row Jacobian entries.
        //  Equals pressureScaleFactor_ today; introduced for symmetry
        //  with the solid-side scaling pattern
        scalar pressureEqnScale_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;

        //- Optional file for recording the time step sizes and number of
        //  nonlinear iterations
        autoPtr<OFstream> tsLogPtr_;


    // Viscous cache

        //- Mesh-dependent part of the least-squares surface fit of the wall
        //  velocity on one fixed-velocity patch (see viscousStressTranspose)
        class wallVelocityFit
        {
        public:

            //- Neighbour faces kept by the sharp-edge filter, per face
            labelListList nbrs;

            //- Tangential distance vectors d to the neighbours (in the order
            //  of nbrs), followed by the wedge images, per face
            List<vectorField> d;

            //- Inverse-square-distance weights w, in the same order
            List<scalarField> w;

            //- Unit normals
            vectorField n;

            //- Inverse of the regularised distance tensor; zero where the
            //  neighbours resolve no tangential direction
            tensorField invDdReg;
        };

        //- Surface-fit geometry of every fixed-velocity patch (unset for the
        //  other patches), for the mesh points fitPoints_
        mutable PtrList<wallVelocityFit> wallFits_;

        //- Wedge transformations, for the images of the wall faces across
        //  the wedge planes
        mutable List<tensor> wedgeT_;

        //- Mesh points the surface-fit geometry was evaluated for
        mutable pointField fitPoints_;

        //- The velocity gradient, the transpose part of the viscous stress,
        //  its divergence, the viscous Laplacian and the time derivative of
        //  the velocity, evaluated together once for each velocity state
        //  (updateViscousCache) and shared by the momentum residual, the
        //  wall pressure gradient and the interface traction
        mutable autoPtr<volTensorField> cachedGradUPtr_;
        mutable autoPtr<volTensorField> cachedNuDevGradUTPtr_;
        mutable autoPtr<volVectorField> cachedDivNuDevGradUTPtr_;
        mutable autoPtr<volVectorField> cachedLaplacianUPtr_;
        mutable autoPtr<volVectorField> cachedDdtUPtr_;

        //- The velocity state the cached fields were evaluated for: the
        //  internal and boundary values of U, the time index and the
        //  time-step size (and the mesh points, fitPoints_)
        mutable vectorField cachedUInternal_;
        mutable List<vectorField> cachedUBoundary_;
        mutable label cachedTimeIndex_;
        mutable scalar cachedDeltaT_;

        //- False after invalidateViscousCache, until the next evaluation
        mutable bool viscousCacheValid_;


    // Private Member Functions

        // Make the rAUf field
        void makeRAUf() const;

        // Return a const reference to the rAU field
        const surfaceScalarField& rAUf() const;

        // Return a const reference to the rAU field
        surfaceScalarField& rAUf();

        //- Correct phi function
        void correctPhi();

        //- True if the surface-fit geometry is for the current mesh points
        bool fitGeometryCurrent() const;

        //- Evaluate the surface-fit geometry (wallFits_, wedgeT_) for the
        //  current mesh points
        void updateWallVelocityFits() const;

        //- Set the wall value of the transpose stress on every boundary with
        //  fixed velocity from the boundary velocity, with the surface-fit
        //  geometry
        void setWallViscousStressTranspose
        (
            volTensorField& nuDevGradUT,
            const volScalarField& nuEff
        ) const;

        //- Transpose part of the viscous stress, nuEff*dev2(T(gradU)), as
        //  in the divDevReff of pimpleFluid. On boundaries with fixed velocity
        //  the boundary value is n*q, set from the wall velocity only:
        //      q = nuEff*((I - n n) & (G_s & n) - n*tr(G_s))
        //  with G_s the least-squares surface gradient of the wall velocity,
        //  i.e. n & grad(U)^T at fixed n with the normal derivative from
        //  continuity (zero velocity divergence)
        tmp<volTensorField> viscousStressTranspose
        (
            const volTensorField& gradU
        ) const;

        //- Evaluate the cached fields if the velocity, the mesh points, the
        //  time index or the time step changed since the last evaluation.
        //  Collective: the decision is agreed by all processors, as the
        //  evaluation communicates across processor boundaries
        void updateViscousCache() const;

        //- Cached fields for the current velocity state
        const volTensorField& cachedGradU() const;
        const volTensorField& cachedViscousStressTranspose() const;
        const volVectorField& cachedDivViscousStressTranspose() const;
        const volVectorField& cachedLaplacianU() const;
        const volVectorField& cachedDdtU() const;

        //- Restore the flow fields and mesh after a failed PETSc solve
        void restoreOldTimeState
        (
            const pointField& oldPoints,
            const bool meshMoved
        );

        //- Solve pressure equation
        void solvePEqn
        (
            const scalar& UUrf,
            const fvVectorMatrix& ddtUEqn,
            const fvVectorMatrix& HUEqn
        );

#ifdef USE_PETSC
        //- Copy the PETSc solution vector x into the U and p fields,
        //  applying the pressure unknown scale
        //  (p = pressureUnknownScale_ * pHat) and correcting boundary
        //  conditions. Vec is a pointer type, so it is taken by value
        void unpackSolution(const Vec x);

        //- Copy the current U and p fields into the PETSc solution
        //  vector x, applying the inverse pressure unknown scale
        //  (pHat = p/pressureUnknownScale_).
        //  Vec is a pointer type, so it is taken by value
        void packSolution(Vec x);

        //- Apply equation under-relaxation to the velocity equation
        //  using an implicit diagonal source term
        void relaxVelocityEqn(fvVectorMatrix& UEqn) const;
#endif

        //- Disallow default bitwise copy construct
        newtonIcoFluid(const newtonIcoFluid&);

        //- Disallow default bitwise assignment
        void operator=(const newtonIcoFluid&);


public:

    //- Runtime type information
    TypeName("newtonIcoFluid");

    // Constructors

        //- Construct from components
        newtonIcoFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~newtonIcoFluid()
        {}


    // Member Functions

        // Access

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Patch viscous force (N/m2) integrated over the reference mesh
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID, const solidModel& motion
            ) const;

            //- Return a reference to the turbulence model
#ifdef OPENFOAM_ORG
            incompressible::momentumTransportModel& turbulence()
#else
            incompressible::turbulenceModel& turbulence()
#endif
            {
                if (turbulence_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not set" << abort(FatalError);
                }

                return autoPtrRef(turbulence_);
            }

            //- Reference to the fluid density
            const dimensionedScalar& rho() const
            {
                return rho_;
            }

            //- Return const reference to the momentum stabilisation
            const stabilisationModel& momentumStabilisation() const
            {
                return momentumStabilisationPtr_();
            }

            //- Return const reference to the pressure stabilisation
            const stabilisationModel& pressureStabilisation() const
            {
                return pressureStabilisationPtr_();
            }

            //- Return the pressure-equation scaling factor
            scalar pressureScaleFactor() const
            {
                return pressureScaleFactor_;
            }

            //- Return the PETSc pressure unknown scale, i.e. the scale s
            //  in pHat = p / s used to store pressure in the SNES vector.
            //  Equals 1 (no-op) unless scaleMixedPetScFields is enabled
            scalar pressureUnknownScale() const
            {
                return pressureUnknownScale_;
            }



        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the fluid model
            virtual bool evolve();

            //- Return a const reference to the rAU field
            void clearRAUf();

            //- Mark the cached viscous fields as out of date. The cache also
            //  checks the velocity, the mesh points, the time index and the
            //  time step itself; this is for changes it cannot see, such as
            //  a new effective viscosity
            void invalidateViscousCache() const;

            //- Set the normal gradient of the (kinematic) pressure on a
            //  patch with fixedGradient pressure and fully prescribed
            //  velocity, at which the fluid moves with the boundary (a
            //  no-slip wall or a fluid-solid interface), from the normal
            //  momentum equation:
            //      dp/dn = n & (laplacian(nuEff, U)
            //                 + div(nuEff*dev2(T(grad(U)))) - dU/dt)
            //  where dU/dt is the acceleration of the boundary faces. The
            //  convective term, which vanishes at a no-slip wall that moves
            //  with the mesh, is neglected. A zeroGradient pressure at such a
            //  boundary gives first-order errors in the boundary pressure and
            //  in the pressure gradient of the cells next to it
            void updatePressureGradientBoundary(const label patchI);

            //- Call updatePressureGradientBoundary on every wall patch with
            //  fixedGradient pressure (not a condition derived from it) and
            //  fully prescribed velocity
            void updatePressureGradientBoundaries();

#ifdef USE_PETSC
            //- Initialise the Jacobian matrix
            //  This procedure should set the size and non-zero structure of
            //  the matrix
            virtual label initialiseJacobian(Mat& jac) override;

            //- Initialise the solution vector
            //  This procedure should set the size of the solution vector
            virtual label initialiseSolution(Vec& x) override;

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector
            virtual label formResidual
            (
                Vec f,         // Residual
                const Vec x    // Solution
            ) override
            {
                return formResidual(f, x, false);
            }

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector, with optional flux extrapolation
            virtual label formResidual
            (
                Vec f,         // Residual
                const Vec x,   // Solution
                const bool extrapolatedFlux
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,       // Jacobian
                const Vec x    // Solution
            ) override
            {
                return formJacobian(jac, x, false);
            }

            //- Form the Jacobian of the governing equation, with optional
            //  flux extrapolation
            virtual label formJacobian
            (
                Mat jac,                    // Jacobian
                const Vec x,                // Solution
                const bool extrapolatedFlux
            );

            //- Apply a SIMPLE-type segregated preconditioner
            //  Approximately solves y = M^{-1} x using one SIMPLE iteration
            //  when PETSc uses the `physics` preconditioner
            virtual label precondition
            (
                Vec y,         // Output: y = M^{-1} x
                const Vec x    // Input: vector to be preconditioned
            ) override;


#endif // USE_PETSC
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
