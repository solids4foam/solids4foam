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
    vertexCentredNonLinGeomTotalLagSolid

Description
    Mathematical model where where non linear geometry is modelled i.e. large strains
    and the total displacement is the primary unknown.

    Discretisation is performed using a vertex-centred finite volume approach
    where the governing equations are integrated over a dual mesh. The dual mesh
    is automatically constructed from the primary mesh using the meshDualiser.

    Three solution are available, selected via solutionAlgorithm option:
      - implicit Newton-Raphson
      - implicit segregated
      - explicit

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

SourceFiles
    vertexCentredNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef vertexCentredNonLinGeomTotalLagSolid_H
#define vertexCentredNonLinGeomTotalLagSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "matNN.H"
#ifdef OPENFOAM_COM
    #include "pointVolInterpolation.H"
#endif
#ifdef USE_PETSC
    #include "foamPetscSnesHelper.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
       Class vertexCentredNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class vertexCentredNonLinGeomTotalLagSolid
:
    public solidModel
#ifdef USE_PETSC
    , public foamPetscSnesHelper
#endif
{
private:

    // Private Data

        //- Mechanical model for dual mesh
        autoPtr<dualMechanicalModel> dualMechanicalPtr_;

        //- Implicit stiffness used by segregated algorithm
        //  Coefficient of the Laplacian term. Defined at the dual mesh faces.
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        mutable autoPtr<surfaceScalarField> dualImpKfPtr_;

        //- Flag to include the geometric stiffness
        const Switch useGeometricStiffness_;

        //- Number of unknowns per cell, used by PETSc SNES solver
        const label blockSize_;

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
        vectorField fixedDofDirectionsVec_;

        //- Scale factor for the fixed degree of freedom equations to improve the
        //  matrix conditioning
        const scalar fixedDofScale_;

#ifdef USE_PETSC
        //- Pointer to PETSc indexing object used for identifying fixed degrees
        //  of freedom
        mutable IS fixedDofRowsISPtr_;
#endif

        //- Point pressure
        pointScalarField pointP_;

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

        //- Deformation gradient field at the dual mesh faces
        surfaceTensorField dualFf_;

        //- Inverse of the deformation gradient field at the dual mesh faces
        surfaceTensorField dualFinvf_;

        //- Jacobian of the deformation gradient field at the dual mesh faces
        surfaceScalarField dualJf_;

        //- Pressure at the dual mesh faces
        surfaceScalarField dualPf_;

        //- Cell pressure
        volScalarField volP_;

#ifdef OPENFOAM_COM
        //- Interpolator from points to cells
        pointVolInterpolation pointVolInterp_;
#endif


    // Private Member Functions

        //- Update the divergence of stress based on the current pointD field
        void updatePointDivSigma
        (
            const pointVectorField& pointD,
            surfaceTensorField& dualGradDf,
            surfaceTensorField& dualFf,
            surfaceTensorField& dualFinvf,
            surfaceScalarField& dualJf,
            surfaceSymmTensorField& dualSigmaf,
            pointVectorField& pointDivSigma
        );

        //- Set fixed degree of freedom lists
        void setFixedDofs
        (
            const pointVectorField& pointD,
            boolList& fixedDofs,
            pointField& fixedDofValues,
            symmTensorField& fixedDofDirections,
            vectorField& fixedDofDirectionsVec
        ) const;

        //- Enforce extact tractions on traction boundaries
        void enforceTractionBoundaries
        (
            const pointVectorField& pointD,
            surfaceVectorField& dualTraction,
            const surfaceVectorField& dualDeformedNormals,
            const fvMesh& mesh, // corresponding to pointD mesh
            const labelListList& pointToDualFaces
        ) const;

#ifdef USE_PETSC
        //- Make fixedDofRowsIS object
        void makeFixedDofRowsIS() const;

        //- Return fixedDofRowsIS object
        IS fixedDofRowsIS() const
        {
            if (fixedDofRowsISPtr_ == nullptr)
            {
                makeFixedDofRowsIS();
            }

            return fixedDofRowsISPtr_;
        }
#endif

        //- Geometric stiffness (tangent)
        //  This is the contribution to the stiffness matrix from geometric
        //  nonlinearity
        // The undeformed surface areas and current displacement gradient are
        // required
        void geometricStiffnessField
        (
            List<mat39>& geometricStiffness,
            const surfaceVectorField& SfUndef,
            const surfaceTensorField& gradDRef
        ) const;

        //- Make dualImpKf field
        void makeDualImpKf() const;

        //- Return const dualImpK field
        const surfaceScalarField& dualImpKf() const;

        //- Return dualImpK field
        surfaceScalarField& dualImpKf();

        //- Predict displacement field using the velocity field
        void predict();

        //- Evolve the solid solver using the PETSc SNES nonlinear solver
        virtual bool evolveSnes();

        //- Evolve the solid solver using the explicit approach
        virtual bool evolveExplicit();

        //- Disallow default bitwise copy construct
        vertexCentredNonLinGeomTotalLagSolid
        (
            const vertexCentredNonLinGeomTotalLagSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const vertexCentredNonLinGeomTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("vertexCentredNonLinTotalLagGeometry");

    // Constructors

        //- Construct from components
        vertexCentredNonLinGeomTotalLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~vertexCentredNonLinGeomTotalLagSolid();


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
#endif // USE_PETSC

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
