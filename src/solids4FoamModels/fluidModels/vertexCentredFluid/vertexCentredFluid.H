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
    vertexCentredFluid

Description
    Mathematical model for an incompressible, laminar Newtonian fluid where
    discretisation is performed using a vertex-centred finite volume approach.
    The governing equations are integrated over a dual mesh. The dual mesh
    is automatically constructed from the primary mesh using the meshDualiser.

    An implicit Newton-Raphson algorithm is employed where the full Jacobian is
    constructed.

    As an initial stepping stone, pressure (pointP) is prescribed as we only
    solve for the velocity.

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

SourceFiles
    vertexCentredFluid.C

\*---------------------------------------------------------------------------*/

#ifndef vertexCentredFluid_H
#define vertexCentredFluid_H

#include "fluidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "sparseMatrix.H"
#include "GeometricField.H"
#include "meshDual.H"
#include "globalPointIndices.H"
#include "Enum.H"
#include "twoDPointCorrector.H"
#include "pointVolInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class vertexCentredFluid Declaration
\*---------------------------------------------------------------------------*/

class vertexCentredFluid
:
    public fluidModel
{
private:

    // Private Data

        //- Point mesh
        const pointMesh& pMesh_;

        //- Dual mesh
        const meshDual dualMesh_;

        //- Compact gradient coefficient
        //  1: use compact gradient stencil
        //  0: use larger stencil
        //  0 < zeta < 1: blended scheme
        const scalar zeta_;

        //- Compact stencil flag for the pressure Laplacian
        const Switch compactStencil_;

        //- Pressure smoothing scaling factor
        const scalar pressureSmoothingCoeff_;

        //- Dynamic viscosity
        const dimensionedScalar mu_;

        //- Bulk (penalty) modulus
        const dimensionedScalar k_;

        //- Density
        const dimensionedScalar rho_;

        //- Implicit stiffness used by segregated algorithm
        //  Coefficient of the Laplacian term. Defined at the dual mesh faces.
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        // mutable autoPtr<surfaceScalarField> dualImpKfPtr_;

        //- Flag to enable full Newton-Raphson approach, where the Jacobian
        //  matrix is updated once per outer iteration. Alternatively, the
        //  stiffness matrix is only updated once per time-step
        const Switch fullNewton_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

        //- Object for correcting the displacement field in 2-D
        //  This ensures the front and back planes have the same displacements
        twoDPointCorrector twoDCorrector_;

        //- Field of flags indicating the degrees of freedom where the value is
        //  set
        //  Currently the size of this list is equal to the number of points
        boolList fixedUDofs_;
        boolList fixedPDofs_;

        //- Values for fixed degrees of freedom
        //  Currently the size of this list is equal to the number of points
        pointField fixedUDofValues_;
        scalarField fixedPDofValues_;

        //- Fixed directions for fixed degrees of freedom
        //  Calculated as sqr(n) where n is a fixed direction, e.g. I for fully
        //  fixed
        symmTensorField fixedUDofDirections_;

        //- Scale factor for the fixed degree of freedom equations to improve the
        //  matrix conditioning
        const scalar fixedUDofScale_;
        const scalar fixedPDofScale_;

        //- Point velocity
        pointVectorField pointU_;

        //- Point pressure field
        pointScalarField pointP_;

        //- Point density field
        //pointScalarField pointRho_;

        //- Point volume field
        //  Dual mesh cell volumes associated with the primary mesh points
        pointScalarField pointVol_;

        // //- Displacement gradient field at the dual mesh faces
        // surfaceTensorField dualGradDf_;

        //- Solution standard tolerance
        const scalar solutionTol_;

        //- Maximum number of outer correctors
        const int nCorr_;

        //- Write frequency for residuals information
        const int infoFrequency_;

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;

        //- Interpolator from points to cells
        pointVolInterpolation pointVolInterp_;


    // Private Member Functions

        //- Calculate the residual of the momentum equation
        tmp<vectorField> residualU
        (
            const pointVectorField& pointU, // velocity
            const pointScalarField& pointP  // pressure
        );

        //- Calculate the residual of the pressure equation
        tmp<scalarField> residualP
        (
            const pointVectorField& pointU, // velocity
            const pointScalarField& pointP  // pressure
        );

        //- Set velocity fixed degree of freedom lists
        void setFixedUDofs
        (
            const pointVectorField& pointU,
            boolList& fixedUDofs,
            pointField& fixedUDofValues,
            symmTensorField& fixedUDofDirections
        ) const;

        //- Set pressure fixed degree of freedom lists
        void setFixedPDofs
        (
            const pointScalarField& pointP,
            boolList& fixedPDofs,
            scalarField& fixedPDofValues
        ) const;

        //- Enforce exact flux on boundaries, e.g. walls and inlets
        void enforceFluxBoundaries
        (
            const pointVectorField& pointU,
            surfaceScalarField& dualFlux,
            const fvMesh& mesh,
            const labelListList& pointToDualFaces
        ) const;

        //- Enforce extact tractions on traction boundaries
        void enforceTractionBoundaries
        (
            const pointScalarField& pointP,
            const pointVectorField& pointU,
            surfaceVectorField& dualTraction,
            const fvMesh& mesh,
            const labelListList& pointToDualFaces
        ) const;

        //- Check convergence of Newton-Raphson loop
        bool converged
        (
            const label iCorr,
            scalar& initResidualU,
            scalar& initResidualP,
            const scalar res,
            const label nInterations,
            const pointVectorField& pointU,
            const pointScalarField& pointP,
            const Field<RectangularMatrix<scalar>>& pointUPcorr
        ) const;

        //- Disallow default bitwise copy construct
        vertexCentredFluid(const vertexCentredFluid&);

        //- Disallow default bitwise assignment
        void operator=(const vertexCentredFluid&);


public:

    //- Runtime type information
    TypeName("vertexCentredFluid");

    // Constructors

        //- Construct from components
        vertexCentredFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~vertexCentredFluid();


    // Member Functions

        // Access

            //- Const-access to the dual mesh
            const meshDual& dualMesh() const
            {
                return dualMesh_;
            }


        // Edit

            //- Update the size of the time-step
            //virtual void setDeltaT(Time& runTime);

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const
            {
                notImplemented("patchViscousForce is work in progress!");

                // Keep the compiler happy
                return tmp<vectorField>(new vectorField());
            }

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const
            {
                notImplemented("patchPressureForce is work in progress!");

                // Keep the compiler happy
                return tmp<scalarField>(new scalarField());
            }

        // Write

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
