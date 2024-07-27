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
#include "sparseMatrix.H"
#include "globalIndex.H"
#ifdef USE_PETSC
    #include <petscksp.h>
    #include <petscsnes.h>
#endif


// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

#ifdef USE_PETSC

// Forward declaration of the Foam::solidModels::linGeomTotalDispSolid class
namespace Foam
{
    namespace solidModels
    {
        class linGeomTotalDispSolid;
    }
}

// User data "context" for PETSc functions
typedef struct appCtx
{
    // Reference to the solid model object
    Foam::solidModels::linGeomTotalDispSolid& solMod_;

    // Constructor
    appCtx
    (
        Foam::solidModels::linGeomTotalDispSolid& solMod
    )
    :
        solMod_(solMod)
    {}
} appCtx;

#endif // ifdef USE_PETSC


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
    public solidModel
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

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Helper object for converting between local (on processor) and global
        //  cell indices
        globalIndex globalCells_;

        //- Flag to set whether the solver stops when PETSs SNES returns an
        //  error
        const Switch stopOnPetscError_;

        //- Pointer to the PETSc SNES object
        autoPtr<SNES> snesPtr_;

        //- Pointer to the PETSc solution vector
        autoPtr<Vec> xPtr_;

        //- Pointer to the PETSc matrix
        autoPtr<Mat> APtr_;

        //- PETSc user context
        autoPtr<appCtx> snesUserPtr_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Evolve the solid solver using the segregated implicit approach
        bool evolveImplicitSegregated();

        //- Evolve the solid solver using the PETSc nonlinear solver (SNES)
        //  A compact approximate Jacobian (same one as used by the segregated
        //  approach) is provided to PETSc. A Jacobian-free Krylov-Newton method
        //  works well with this approach
        bool evolveSnes();

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

            //- Const access to the global index object
            const globalIndex& globalCells() const
            {
                return globalCells_;
            }


        // Edit

            //- Calculate the momentum equation residual
            //  This should be const but is not because of the mechanical law
            //  The can be fixed when the mechanical models are refactored
            tmp<vectorField> residualMomentum
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
            tmp<sparseMatrix> JacobianMomentum
            (
                const volVectorField& D
            );

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- This function is called at the end of the simulation
            virtual void end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
