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
    Philip Cardiff, UCD.  All rights reserved.

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

        //- Acceleration field
        //  Only used by the explicit solution scheme
        volVectorField A_;

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

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

            //- Calculate the momentum equation residual
            //  This should be const but is not because of the mechanical law
            //  The can be fixed when the mechanical models are refactored
            virtual tmp<vectorField> residualMomentum
            (
                const volVectorField& D
            );

            //- Calculate the Jacobian for the momentum equation
            //  This Jacobian is the derivative of residualMomentum with respect
            //  to the solution vector (D)
            //  This should be const but is not because of the mechanical law
            //  The can be fixed when the mechanical models are refactored
            //  This can also be an approximation of the true Jacobian, e.g. a
            //  compact molecule Laplacian
            //  We initially adopt this approach for testin
            virtual tmp<sparseMatrix> JacobianMomentum
            (
                const volVectorField& D
            );

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
