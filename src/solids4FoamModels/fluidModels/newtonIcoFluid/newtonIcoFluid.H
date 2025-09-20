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

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;


    // Private Member Functions

        // Make the rAUf field
        void makeRAUf() const;

        // Return a const reference to the rAU field
        const surfaceScalarField& rAUf() const;

        // Return a const reference to the rAU field
        surfaceScalarField& rAUf();

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


        // Edit

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
