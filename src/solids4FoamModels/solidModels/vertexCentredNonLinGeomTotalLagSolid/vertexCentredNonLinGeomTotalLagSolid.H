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

Author
    Philip Cardiff, UCD.  All rights reserved.

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
#include "sparseMatrix.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "globalPointIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class vertexCentredLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class vertexCentredNonLinGeomTotalLagSolid
:
    public solidModel
{
    // Private Data

        //- Mechanical model for dual mesh
        autoPtr<dualMechanicalModel> dualMechanicalPtr_;

        //- Flag to enable full Newton-Raphson approach, where the stiffness
        //  matrix is updated once per outer iteration. Alternatively, the
        //  stiffness matrix is only updated once per time-step
        const Switch fullNewton_;

        //- Flag to enable steadyState
        Switch steadyState_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

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

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;


    // Private Member Functions

        //- Update the source vector for the linear system
        void updateSource
        (
            vectorField& source,
            const labelList& dualCellToPoint
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
        
        Foam::tmp<Foam::Field<Foam::RectangularMatrix<Foam::scalar>>>
		geometricStiffnessField
		(
			const surfaceVectorField SfUndef, //Undeformed surface area vector field
			const surfaceTensorField gradDRef //Reference gradD 
		) const; 

        //- Disallow default bitwise copy construct
        vertexCentredNonLinGeomTotalLagSolid(const vertexCentredNonLinGeomTotalLagSolid&);

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
