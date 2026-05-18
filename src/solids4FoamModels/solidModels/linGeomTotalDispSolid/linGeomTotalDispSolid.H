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
    linGeomTotalDispSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement is the primary
    unknown.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.

SourceFiles
    linGeomTotalDispSolid.C

\*---------------------------------------------------------------------------*/

#ifndef linGeomTotalDispSolid_H
#define linGeomTotalDispSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "foamPetscSnesHelper.H"
#include "hofvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linGeomTotalDispSolid Declaration
\*---------------------------------------------------------------------------*/

class linGeomTotalDispSolid
:
    public solidModel,
    public foamPetscSnesHelper
{
    // Private data

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Reciprocal of the bulk modulus
        //  Used in the pressure Jacobian when solvePressure is enabled
        const volScalarField rKappa_;

        //- Acceleration field
        //  Only used by the explicit solution scheme
        volVectorField A_;

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Cached face unit-normal vectors (= mesh.Sf()/mesh.magSf()).
        //  Purely geometric: independent of D and p. Rebuilt only when
        //  the mesh has moved between timesteps
        mutable autoPtr<surfaceVectorField> nfPtr_;

        //- Time index when nfPtr_ was last refreshed (-1 = never built).
        //  Used to detect mesh motion between timesteps
        mutable label nfTimeIndex_;

        //- Time index when solidModel::rAUf() was last computed
        //  (-1 = never built)
        mutable label rAUfTimeIndex_;

        //- deltaT for which solidModel::rAUf() was last computed.
        //  rAUf depends only on impKf_, rho, mesh and deltaT (not on
        //  the values of D or p), so it is safe to reuse across PETSc
        //  matrix-free residual evaluations within a Newton solve
        mutable scalar rAUfDeltaT_;

        //- Store pressure in the PETSc solution vector as p/pScale
        const Switch scaleMixedPetScFields_;

        //- Scale type used for the PETSc pressure unknown
        const word pressureUnknownScaleType_;

        //- Pressure scale used for pHat = p/pressureUnknownScale_
        scalar pressureUnknownScale_;

        //- User-facing dimensionless scale factor applied to the
        //  pressure-equation residual and Jacobian rows relative to
        //  the momentum equations
        const scalar pressureScaleFactor_;

        //- If true, multiply pressureScaleFactor_ by a representative
        //  2*mu so that the pressure-equation residual has units
        //  comparable to the momentum equation
        const Switch pressureScaleByTwoMu_;

        //- Representative value of 2*mu used as the physical scale of
        //  the pressure equation. Computed once at construction
        scalar twoMuRef_;

        //- Effective scale applied to the pressure residual and the
        //  pressure-row Jacobian entries.
        //  pressureEqnScale_ = pressureScaleFactor_
        //                       * (pressureScaleByTwoMu_ ? twoMuRef_ : 1)
        scalar pressureEqnScale_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;

        //- Field to provide fvOptions with the correct dimensions
        volScalarField ds_;

    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Return the cached face unit-normal field (Sf()/magSf()).
        //  Built lazily and refreshed only when the mesh has moved
        //  between timesteps. Independent of D and p, so safe under
        //  PETSc matrix-free perturbation
        const surfaceVectorField& nf() const;

        //- Update solidModel::rAUf() if the mesh or deltaT changed
        //  since it was last computed. The result depends only on
        //  impKf_, rho, mesh and deltaT (not on the values of D or p),
        //  so it can safely be reused across MFFD residual evaluations
        //  within a single Newton solve. Does nothing when
        //  solvePressure() is disabled
        void updateRAUfIfStale();

#ifdef USE_PETSC
        //- Copy the PETSc solution vector x into the D field (and p
        //  when solvePressure() is active), applying the pressure
        //  unknown scale (p = pressureUnknownScale_ * pHat) and
        //  correcting boundary conditions. Centralises the mapping
        //  used by formResidual, formJacobian and evolveSnes.
        //  Vec is a pointer type, so it is taken by value
        void unpackSolution(const Vec x);

        //- Copy the current D field (and, when solvePressure() is
        //  active, p) into the PETSc solution vector x, applying the
        //  inverse pressure unknown scale (pHat = p/pressureUnknownScale_).
        //  Vec is a pointer type, so it is taken by value
        void packSolution(Vec x);
#endif

        //- Enforce traction boundary conditions on the provided traction field
        void enforceTractionBoundaries
        (
            surfaceVectorField& traction,
            const volVectorField& D,      // displacement field
            const surfaceVectorField& n   // unit normals field
        ) const;

        //- Evolve the solid solver using the segregated implicit approach
        bool evolveImplicitSegregated();

        //- Evolve the solid solver using the PETSc nonlinear solver (SNES)
        //  A compact approximate Jacobian (same one as used by the segregated
        //  approach) is provided to PETSc. A Jacobian-free Krylov-Newton method
        //  works well with this approach
        bool evolveSnes();

        //- Evolve the solid solver using the explicit approach
        bool evolveExplicit();

        //- Disallow default bitwise copy construct
        linGeomTotalDispSolid(const linGeomTotalDispSolid&);

        //- Disallow default bitwise assignment
        void operator=(const linGeomTotalDispSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("linearGeometryTotalDisplacement");

    // Constructors

        //- Construct from components
        linGeomTotalDispSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~linGeomTotalDispSolid()
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

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the solid solver and solve the mathematical model
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

            //- Apply a physics-based preconditioner
            virtual label precondition
            (
                Vec y,         // Output: y = M^{-1} x
                const Vec x    // Input: vector to be preconditioned
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
