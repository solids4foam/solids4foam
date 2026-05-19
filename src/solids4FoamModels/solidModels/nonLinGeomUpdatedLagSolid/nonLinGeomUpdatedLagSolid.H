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
    nonLinGeomUpdatedLagSolid

Description
    Mathematical model for finite strains and rotations, where an updated
    Lagrangian and incremental displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    More details in:
    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

Author
    Philip Cardiff, UCD.

SourceFiles
    nonLinGeomUpdatedLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomUpdatedLagSolid_H
#define nonLinGeomUpdatedLagSolid_H

#include "solidModel.H"
#include "foamPetscSnesHelper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class nonLinGeomUpdatedLagSolid Declaration
\*---------------------------------------------------------------------------*/

class nonLinGeomUpdatedLagSolid
:
    public solidModel,
    public foamPetscSnesHelper
{
    // Private data

        //- Total deformation gradient
        volTensorField F_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volTensorField relFinv_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volScalarField relJ_;

        //- Density (updated configuration)
        volScalarField rho_;

        //- Acceleration
        volVectorField A_;

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

        //- Time index when solidModel::rAUf() was last computed
        //  (-1 = never built)
        mutable label rAUfTimeIndex_;

        //- deltaT for which solidModel::rAUf() was last computed.
        //  rAUf depends only on impKf_, rho, mesh and deltaT (not on
        //  the values of DD or p), so it is safe to reuse across PETSc
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

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Allocate the reciprocal of the bulk modulus
        void makeRKappa() const;

        //- Return const access to the reciprocal of the bulk modulus
        const volScalarField& rKappa() const;

        //- Update solidModel::rAUf() if the mesh or deltaT changed
        //  since it was last computed. The result depends only on
        //  impKf_, rho, mesh and deltaT (not on the values of DD or p),
        //  so it can safely be reused across MFFD residual evaluations
        //  within a single Newton solve. Does nothing when
        //  solvePressure() is disabled
        void updateRAUfIfStale();

#ifdef USE_PETSC
        //- Copy the PETSc solution vector x into the DD field (and p
        //  when solvePressure() is active), refresh the dependent
        //  kinematic fields (D, gradDD, relF, relFinv, F, relJ, J,
        //  sigma) and correct boundary conditions. Applies the
        //  pressure unknown scale (p = pressureUnknownScale_ * pHat).
        //  Centralises the mapping used by formResidual, formJacobian
        //  and evolveSnes. Vec is a pointer type, so it is taken by
        //  value
        void unpackSolution(const Vec x);

        //- Copy the current DD field (and, when solvePressure() is
        //  active, p) into the PETSc solution vector x, applying the
        //  inverse pressure unknown scale (pHat = p/pressureUnknownScale_).
        //  Vec is a pointer type, so it is taken by value
        void packSolution(Vec x);
#endif

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
        nonLinGeomUpdatedLagSolid(const nonLinGeomUpdatedLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const nonLinGeomUpdatedLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::UPDATED_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("nonLinearGeometryUpdatedLagrangian");

    // Constructors

        //- Construct from components
        nonLinGeomUpdatedLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~nonLinGeomUpdatedLagSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for DD
                return DD();
            }

            //- This solid model solves for DD instead of D
            virtual bool incremental() const
            {
                return true;
            }

            //- Does this model move the mesh?
            virtual bool movingMesh() const
            {
                return true;
            }

        // Edit

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
                Vec f,        // Residual
                const Vec x   // Solution
            );

            //- Form the Jacobian of the governing equation
            virtual label formJacobian
            (
                Mat jac,              // Jacobian
                const Vec x  // Solution
            );
#endif // USE_PETSC

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
