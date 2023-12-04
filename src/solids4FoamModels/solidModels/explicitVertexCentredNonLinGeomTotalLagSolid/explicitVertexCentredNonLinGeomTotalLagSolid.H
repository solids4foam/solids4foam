/*---------------------------------------------------------------------------* \
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
    explicitVertexCentredNonLinGeomTotalLagSolid

Description
    Mathematical model where nonnlinear geometry is assumed using a total
    Lagrangian formulation. The point total displacement is the primary unknown.

    Discretisation is performed using a vertex-centred finite volume approach
    where the governing equations are integrated over a dual mesh. The dual mesh
    is automatically constructed from the primary mesh using the meshDualiser.

    An explicit solution algorithm is employed in contrast to the implicit
    solution algorithm used in vertexCentredLinGeomSolid.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    explicitVertexCentredNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitVertexCentredNonLinGeomTotalLagSolid_H
#define explicitVertexCentredNonLinGeomTotalLagSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "globalPointIndices.H"
#ifdef OPENFOAM_COM
    #include "pointVolInterpolation.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
         Class explicitVertexCentredNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitVertexCentredNonLinGeomTotalLagSolid
:
    public solidModel
{
    // Private Data

        //- Store stiffness field: only used by pointSolidContact boundaries for
        //  calculating the penalty factor
        const volScalarField impK_;

        //- Use artificial compressibility approach to calculate the pressure
        const Switch useArtificialCompressibility_;

        //- Use deformation gradient involution equation
        const Switch useDeformationGradientInvolution_;

        //- Compact edge gradient factor
        const scalar zeta_;

        //- Damping coefficient
        const dimensionedScalar dampingCoeff_;

        //- Mechanical model for dual mesh
        autoPtr<dualMechanicalModel> dualMechanicalPtr_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

        //- Bulk modulus at the dual mesh faces
        surfaceScalarField dualKf_;

        //- Point velocity
        pointVectorField pointU_;

        //- Point acceleration field
        pointVectorField pointA_;

        //- Point density field
        pointScalarField pointRho_;

        //- Point volume field
        //  Dual mesh cell volumes associated with the primary mesh points
        pointScalarField pointVol_;

        //- Divergence of sigma at each primary point
        pointVectorField pointDivSigma_;

        //- Displacement gradient field at the dual mesh faces
        surfaceTensorField dualGradDf_;

        //- Deformation gradient field at the dual mesh faces
        surfaceTensorField dualFf_;

        //- Inverse deformation gradient field at the dual mesh faces
        surfaceTensorField dualFinvf_;

        //- Jacobian of the deformation gradient field at the dual mesh faces
        surfaceScalarField dualJf_;

        //- Hydrostatic pressure at the dual mesh faces
        surfaceScalarField dualPf_;

        //- Stress field at the dual mesh faces
        surfaceSymmTensorField dualSigmaf_;

        //- Number of processors which share this point
        // Interior points will be 1, but points on processor boundaries will be
        //  at least 2
        pointScalarField nProcPerPoint_;

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;

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
            surfaceScalarField& dualPf,
            surfaceSymmTensorField& dualSigmaf,
            pointVectorField& pointDivSigma
        );

        //- Enforce extact tractions on traction boundaries
        void enforceTractionBoundaries
        (
            const pointVectorField& pointD,
            surfaceVectorField& dualTraction,
            const surfaceVectorField& dualDeformedNormals,
            const fvMesh& mesh, // corresponding to pointD mesh
            const labelListList& pointToDualFaces
        ) const;

        //- Disallow default bitwise copy construct
        explicitVertexCentredNonLinGeomTotalLagSolid(const explicitVertexCentredNonLinGeomTotalLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const explicitVertexCentredNonLinGeomTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("explicitVertexCentredNonLinearGeometryTotalLagrangian");

    // Constructors

        //- Construct from components
        explicitVertexCentredNonLinGeomTotalLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~explicitVertexCentredNonLinGeomTotalLagSolid();


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
