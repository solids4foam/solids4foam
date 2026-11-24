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

    Optional entries in newtonIcoFluidCoeffs:

        // Include the convection term (default yes). Without it, the model
        // solves the unsteady Stokes equations
        convection      no;

        // High-order (k-exact least-squares) residual, in the format of the
        // solid models' highOrderCoeffs with velocity in place of
        // displacement. The stress and the volume flux are integrated over
        // the face quadrature points of the reconstructions of U and p, and
        // the second-order Jacobian is kept as the preconditioner. It needs
        // convection off, a static mesh, a serial run and a pressure
        // stabilisation that supports it, e.g. alpha
        highOrderCoeffs
        {
            highOrderResidual yes;

            velocity
            {
                type                    kExactLeastSquares;
                polynomialOrder         3;
                faceStencilExtraCells   20;
                weightFunctionCoeffs
                {
                    type    radiallySymmetricExponential;
                    k       6;
                }
            }

            pressure
            {
                // As for velocity
            }
        }

    An optional momentumSource field (m/s^2) in the start time directory is
    added to the momentum equation, e.g. for a manufactured solution. Its cell
    values are the cell averages of the body force per unit mass, and they do
    not change in time.

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
#include "leastSquaresReconstruction.H"


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

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;

        //- Include the convection term; without it, the model solves the
        //  unsteady Stokes equations
        const Switch convection_;

        //- Use the high-order (k-exact least-squares) residual
        Switch highOrderResidual_;

        //- Optional momentum source per unit mass, read from the start time
        //  directory; its cell values are cell averages and do not change in
        //  time
        autoPtr<volVectorField> momentumSourcePtr_;

        //- Optional file for recording the time step sizes and number of
        //  nonlinear iterations
        autoPtr<OFstream> tsLogPtr_;


    // Private Member Functions

        // Make the rAUf field
        void makeRAUf() const;

        // Return a const reference to the rAU field
        const surfaceScalarField& rAUf() const;

        // Return a const reference to the rAU field
        surfaceScalarField& rAUf();

        //- Return the highOrderCoeffs sub-dictionary for the given field
        //  role (velocity or pressure)
        const dictionary& highOrderCoeffs(const word& fieldRole) const;

        //- Return true for the patches on which the field fixes the value,
        //  which are included in the least-squares reconstruction stencils
        template<class Type>
        static boolList fixedValuePatchMask
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        );

        //- Calculate the face traction (the stress integrated over each face
        //  with the face quadrature, divided by the face area) and the face
        //  volume flux from the high-order reconstructions of U and p
        void highOrderFaceFluxes
        (
            const volVectorField& U,
            const volScalarField& p,
            surfaceVectorField& traction,
            surfaceScalarField& flux
        ) const;

        //- Correct phi function
        void correctPhi();

        //- Solve pressure equation
        void solvePEqn
        (
            const scalar& UUrf,
            const fvVectorMatrix& ddtUEqn,
            const fvVectorMatrix& HUEqn
        );

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

            //- Return true if the high-order residual is active
            bool highOrderResidual() const
            {
                return highOrderResidual_;
            }

            //- Return the least-squares reconstruction scheme of U
            const leastSquaresScheme& velocityLeastSquares() const;

            //- Return the least-squares reconstruction scheme of p
            const leastSquaresScheme& pressureLeastSquares() const;


        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the fluid model
            virtual bool evolve();

            //- Return a const reference to the rAU field
            void clearRAUf();

#ifdef USE_PETSC
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
                Vec f,         // Residual
                const Vec x    // Solution
            );

            //- Calculate the residual (PETSc array) given the solution (PETSc
            //  array) vector, given a reference configuration (solid model
            //  representing the fluid domain motion)
            virtual label formResidual
            (
                Vec f,         // Residual
                const Vec x,    // Solution
                const solidModel& motion
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,       // Jacobian
                const Vec x    // Solution
            );
#endif // USE_PETSC
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
