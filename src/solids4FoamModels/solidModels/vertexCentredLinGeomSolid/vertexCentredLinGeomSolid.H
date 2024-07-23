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
    vertexCentredLinGeomSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement is the primary
    unknown.

    Discretisation is performed using a vertex-centred finite volume approach
    where the governing equations are integrated over a dual mesh. The dual mesh
    is automatically constructed from the primary mesh using the meshDualiser.

    Three solution are available, selected via solutionAlgorithm option:
      - implicit Newton-Raphson
      - implicit segregated
      - explicit
      - PETSc SNES (Scalable Nonlinear Equations Solvers)

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

SourceFiles
    vertexCentredLinGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef vertexCentredLinGeomSolid_H
#define vertexCentredLinGeomSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "sparseMatrix.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "globalPointIndices.H"
#include "Enum.H"
#include "twoDPointCorrector.H"
#include "pointVolInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                  Class vertexCentredLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class vertexCentredLinGeomSolid
:
    public solidModel
{
public:

        //- Solution algorithms
        enum solutionAlgorithm
        {
            PETSC_SNES,
            IMPLICIT_COUPLED,
            IMPLICIT_SEGREGATED,
            EXPLICIT
        };

        //- Solution algorithm names
        static const Enum<solutionAlgorithm> solutionAlgorithmNames_;


private:

    // Private Data

        //- Mechanical model for dual mesh
        autoPtr<dualMechanicalModel> dualMechanicalPtr_;

        //- Solution algorithm, e.g. coupled, segregated, explicit, SNES
        const solutionAlgorithm solutionAlgorithm_;

        //- Implicit stiffness used by segregated algorithm
        //  Coefficient of the Laplacian term. Defined at the dual mesh faces.
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        mutable autoPtr<surfaceScalarField> dualImpKfPtr_;

        // Compact edge gradient factor
        const scalar zeta_;

        //- Flag to enable full Newton-Raphson approach, where the stiffness
        //  matrix is updated once per outer iteration. Alternatively, the
        //  stiffness matrix is only updated once per time-step
        const Switch fullNewton_;

        //- Flag to enable steadyState
        Switch steadyState_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

        //- Object for correcting the displacement field in 2-D
        //  This ensures the front and back planes have the same displacements
        twoDPointCorrector twoDCorrector_;

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;

        //- Field of flags indicating the degrees of freedom where the value is
        //  set
        //  Currently the size of this list is equal to the number of points
        boolList fixedDofs_;

        //- Values for fixed degrees of freedom
        //  Currently the size of this list is equal to the number of points
        pointField fixedDofValues_;

        //- Fixed directions for fixed degrees of freedom
        //  Calculated as sqr(n) where n is a fixed direction, e.g. I for fully
        //  fixed
        symmTensorField fixedDofDirections_;

        //- Scale factor for the fixed degree of freedom equations to improve the
        //  matrix conditioning
        const scalar fixedDofScale_;

        //- Point velocity
        pointVectorField pointU_;

        //- Point acceleration field
        pointVectorField pointA_;

        //- Point density field
        pointScalarField pointRho_;

        //- Point volume field
        //  Dual mesh cell volumes associated with the primary mesh points
        pointScalarField pointVol_;

        //- Global point volume field
        //  Dual mesh cell volumes associated with the primary mesh points,
        //  where point volumes on processor boundaries are summed
        pointScalarField pointGlobalVol_;

        //- Divergence of sigma at each primary point
        pointVectorField pointDivSigma_;

        //- Displacement gradient field at the dual mesh faces
        surfaceTensorField dualGradDf_;

        //- Stress field at the dual mesh faces
        surfaceSymmTensorField dualSigmaf_;

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;

        //- Interpolator from points to cells
        pointVolInterpolation pointVolInterp_;


    // Private Member Functions

        //- Update the divergence of stress based on the current pointD field
        void updatePointDivSigma
        (
            const pointVectorField& pointD,
            surfaceTensorField& dualGradDf,
            surfaceSymmTensorField& dualSigmaf,
            pointVectorField& pointDivSigma
        );

        //- Set fixed degree of freedom lists
        void setFixedDofs
        (
            const pointVectorField& pointD,
            boolList& fixedDofs,
            pointField& fixedDofValues,
            symmTensorField& fixedDofDirections
        ) const;

        //- Enforce extact tractions on traction boundaries
        void enforceTractionBoundaries
        (
            const pointVectorField& pointD,
            surfaceVectorField& dualTraction,
            const fvMesh& mesh, // corresponding to pointD mesh
            const labelListList& pointToDualFaces
        ) const;

        //- Check convergence of Newton-Raphson loop
        bool converged
        (
            const label iCorr,
            scalar& initResidual,
            const scalar res,
            const label nInterations,
            const pointVectorField& pointD,
            const vectorField& pointDcorr
        ) const;

        //- Calculate line search slope for the given eta parameter
        scalar calculateLineSearchSlope
        (
            const scalar eta,
            const vectorField& pointDcorr,
            pointVectorField& pointD,
            surfaceTensorField& dualGradDf,
            surfaceSymmTensorField& dualSigmaf,
            const scalar zeta
        );

        scalar calculateLineSearchFactor
        (
            const scalar rTol, // Slope reduction tolerance
            const int maxIter, // Maximum number of line search iterations
            const vectorField& pointDcorr, // Point displacement correction
            const vectorField& source, // Linear system source
            const scalar zeta // Discretisation parameter
        );

        //- Make dualImpKf field
        void makeDualImpKf() const;

        //- Return dualImpK field
        const surfaceScalarField& dualImpKf() const;

        //- Predict displacement field using the velocity field
        void predict();

        //- Evolve the solid solver using the PETSc "SNES" nonlinear solver
        //  SNES (Scalable Nonlinear Equations Solvers) is a general interface
        //  to a large variety on nonlinear solution algorithms
        virtual bool evolveSnes();

        //- Evolve the solid solver using the implicit coupled approach
        virtual bool evolveImplicitCoupled();

        //- Evolve the solid solver using the implicit segregated approach
        virtual bool evolveImplicitSegregated();

        //- Evolve the solid solver using the explicit approach
        virtual bool evolveExplicit();

        //- Disallow default bitwise copy construct
        vertexCentredLinGeomSolid(const vertexCentredLinGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const vertexCentredLinGeomSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("vertexCentredLinearGeometry");

    // Constructors

        //- Construct from components
        vertexCentredLinGeomSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~vertexCentredLinGeomSolid();


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                notImplemented("For vertex-centred approach");

                // This model solves for D
                return D();
            }

            //- Const-access the fixed DOFs list
            const globalPointIndices& gPointIndices() const
            {
                return globalPointIndices_;
            }

            //- Is the model two-dimensional
            const boolList& fixedDofs() const
            {
                return fixedDofs_;
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
                const pointVectorField& pointD
            );

            //- Calculate the Jacobian for the momentum equation
            //  This Jacobian is the derivative of residualMomentum with respect
            //  to the solution vector (pointD)
            //  This should be const but is not because of the mechanical law
            //  The can be fixed when the mechanical models are refactored
            tmp<sparseMatrix> JacobianMomentum
            (
                const pointVectorField& pointD
            );

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const
            {
                notImplemented("For vertex-centred approach");

                // Keep compiler happy
                return tmp<vectorField>(new vectorField());
            }

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label interfaceI,
                const label patchID,
                const vectorField& faceZoneTraction
            );

        // Write

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
