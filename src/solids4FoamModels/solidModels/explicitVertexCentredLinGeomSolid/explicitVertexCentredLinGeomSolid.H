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
    explicitVertexCentredLinGeomSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement is the primary
    unknown.

    Discretisation is performed using a vertex-centred finite volume approach
    where the governing equations are integrated over a dual mesh. The dual mesh
    is automatically constructed from the primary mesh using the meshDualiser.

    An explicit solution algorithm is employed in contrast to the implicit
    solution algorithm used in vertexCentredLinGeomSolid.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    explicitVertexCentredLinGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitVertexCentredLinGeomSolid_H
#define explicitVertexCentredLinGeomSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "globalPointIndices.H"
#ifdef OPENFOAMESI
    #include "pointVolInterpolation.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class explicitVertexCentredLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitVertexCentredLinGeomSolid
:
    public solidModel
{
    // Private Data

        //- Compact edge gradient factor
        const scalar zeta_;

        //- Damping coefficient
        const dimensionedScalar dampingCoeff_;

        //- Mechanical model for dual mesh
        autoPtr<dualMechanicalModel> dualMechanicalPtr_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

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

        //- Stress field at the dual mesh faces
        surfaceSymmTensorField dualSigmaf_;

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;

#ifdef OPENFOAMESI
        //- Interpolator from points to cells
        pointVolInterpolation pointVolInterp_;
#endif


    // Private Member Functions

        //- Update the divergence of stress based on the current pointD field
        void updatePointDivSigma
        (
            const pointVectorField& pointD,
            surfaceTensorField& dualGradDf,
            surfaceSymmTensorField& dualSigmaf,
            pointVectorField& pointDivSigma
        );

        //- Enforce extact tractions on traction boundaries
        void enforceTractionBoundaries
        (
            const pointVectorField& pointD,
            surfaceVectorField& dualTraction,
            const fvMesh& mesh, // corresponding to pointD mesh
            const labelListList& pointToDualFaces
        ) const;

        //- Disallow default bitwise copy construct
        explicitVertexCentredLinGeomSolid(const explicitVertexCentredLinGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const explicitVertexCentredLinGeomSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("explicitVertexCentredLinearGeometry");

    // Constructors

        //- Construct from components
        explicitVertexCentredLinGeomSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~explicitVertexCentredLinGeomSolid();


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
