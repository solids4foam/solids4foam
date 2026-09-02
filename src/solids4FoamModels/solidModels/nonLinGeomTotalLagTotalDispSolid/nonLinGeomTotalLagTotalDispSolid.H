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
    nonLinGeomTotalLagTotalDispSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and total displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    More details in:
    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

Author
    Philip Cardiff, UCD.

SourceFiles
    nonLinGeomTotalLagTotalDispSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomTotalLagTotalDispSolid_H
#define nonLinGeomTotalLagTotalDispSolid_H

#include "solidModel.H"
#include "mechanicalConstitutiveLawManager.H"
#include "foamPetscSnesHelper.H"
#include "hofvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class nonLinGeomTotalLagTotalDispSolid Declaration
\*---------------------------------------------------------------------------*/

class nonLinGeomTotalLagTotalDispSolid
:
    public solidModel,
    public foamPetscSnesHelper
{
    // Private data

        //- Total deformation gradient
        volTensorField F_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Acceleration
        volVectorField A_;

        //- Take the implicit stiffness from the mechanicalConstitutiveLaw
        //  framework rather than from the legacy mechanicalModel.
        //  Declared before impK_ because makeImpK() consults it
        const Switch useMechanicalConstitutiveLawManager_;

        //- The mechanicalConstitutiveLaw framework, built on demand
        mutable autoPtr<mechanicalConstitutiveLawManager> mechanicalManagerPtr_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        const surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Reciprocal of the bulk modulus
        //  Used in the pressure equation when solvePressure is enabled
        mutable autoPtr<volScalarField> rKappaPtr_;

        //- Time rate of pressure
        //  This is used by the optional linear predictor
        mutable autoPtr<volScalarField> dpdtPtr_;

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;


    // Private Member Functions

        //- The mechanicalConstitutiveLaw framework, built on first use
        mechanicalConstitutiveLawManager& mechanicalManager() const;

        //- Update the cell-centred stress, from the framework when
        //  useMechanicalConstitutiveLawManager_ is set and from the legacy
        //  mechanicalModel otherwise
        void correctStress();

        //- Construct impK from the framework when
        //  useMechanicalConstitutiveLawManager_ is set, and from the legacy
        //  mechanicalModel otherwise
        tmp<volScalarField> makeImpK() const;

        //- Calculate the current face-area vectors
        tmp<surfaceVectorField> currentSf() const;

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Allocate the reciprocal of the bulk modulus
        void makeRKappa() const;

        //- Return const access to the reciprocal of the bulk modulus
        const volScalarField& rKappa() const;

        //- Enforce traction boundary conditions on the provided traction field
        void enforceTractionBoundaries
        (
            surfaceVectorField& force,             // force at the faces
            const volVectorField& D,               // displacement field
            const surfaceVectorField& nCurrent,    // deformed unit normals
            const surfaceScalarField& magSfCurrent // deformed area magnitudes
        ) const;

        //- Evolve the solid solver using the segregated implicit approach
        bool evolveImplicitSegregated();

        //- Evolve the solid solver using the PETSc nonlinear solver (SNES)
        //  A compact approximate Jacobian (same one as used by the segregated
        //  approach) is provided to PETSc. A Jacobian-free Krylov-Newton method
        //  works well with this approach
        bool evolveSnes();

        //- Disallow default bitwise copy construct
        nonLinGeomTotalLagTotalDispSolid
        (
            const nonLinGeomTotalLagTotalDispSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const nonLinGeomTotalLagTotalDispSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("nonLinearGeometryTotalLagrangianTotalDisplacement");

    // Constructors

        //- Construct from components
        nonLinGeomTotalLagTotalDispSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~nonLinGeomTotalLagTotalDispSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                return D();
            }

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            //- Let the framework know the time step has ended, so a law
            //  can report or commit anything it keeps per step
            virtual void updateTotalFields() override;

            virtual bool evolve();

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

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,       // Jacobian
                const Vec x    // Solution
            );
#endif // USE_PETSC

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
