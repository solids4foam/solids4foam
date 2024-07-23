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
    linGeomTotalDispSnesSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement is the primary
    unknown.

    The stress is calculated by the run-time selectable mechanical law.

    The PETSc SNES (Scalable Nonlinear Equations Solvers) nonlinear solver is
    used to solve the governing equation. SNES allows run-time selection of the
    solution algorithm.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    linGeomTotalDispSnesSolid.C

\*---------------------------------------------------------------------------*/

#ifndef linGeomTotalDispSnesSolid_H
#define linGeomTotalDispSnesSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "twoDPointCorrector.H"
#include "sparseMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linGeomTotalDispSnesSolid Declaration
\*---------------------------------------------------------------------------*/

class linGeomTotalDispSnesSolid
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

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Disallow default bitwise copy construct
        linGeomTotalDispSnesSolid(const linGeomTotalDispSnesSolid&);

        //- Disallow default bitwise assignment
        void operator=(const linGeomTotalDispSnesSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("linearGeometryTotalDisplacementSnes");

    // Constructors

        //- Construct from components
        linGeomTotalDispSnesSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~linGeomTotalDispSnesSolid()
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

            //- Is the model two-dimensional
            const bool& twoD() const
            {
                return twoD_;
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
