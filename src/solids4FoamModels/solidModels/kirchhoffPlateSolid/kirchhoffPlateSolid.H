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
    kirchhoffPlateSolid

Description
    Kirchhoff plate solver.

    Based on "rotation-free formulation" of Torlak (2006).

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    kirchhoffPlateSolid.C

\*---------------------------------------------------------------------------*/

#ifndef kirchhoffPlateSolid_H
#define kirchhoffPlateSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "areaFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class kirchhoffPlateSolid Declaration
\*---------------------------------------------------------------------------*/

class kirchhoffPlateSolid
:
    public solidModel
{
    // Private data

        //- Finite Area mesh
        faMesh aMesh_;

        //- Transverse displacement
        //  displacement prependicular to the plane of the plate
        areaScalarField w_;

        //- Transverse displacement volume field for visualisation
        volScalarField wVf_;

        //- Moment sum
        areaScalarField M_;

        //- Moment sum volume field for visualisation
        volScalarField MVf_;

        //- Net pressure field on the plate (perpendicular to the plate normal)
        areaScalarField p_;

        //- Net pressure field on the plate volume field
        volScalarField pVf_;

        //- Angle of rotation
        areaVectorField theta_;

        //- Angle of rotation volume field for visualisation
        volVectorField thetaVf_;

        //- Gradient of angle of rotation
        areaTensorField gradTheta_;

        //- Plate density
        dimensionedScalar rho_;

        //- Plate Young's modulus
        dimensionedScalar E_;

        //- Plate Poisson's ratio
        dimensionedScalar nu_;

        //- Plate thickness
        const dimensionedScalar h_;

        //- Plate bending stiffness
        dimensionedScalar bendingStiffness_;

        //- areaPatch index
        //  This is the patch on the fvMesh mesh that corresponds to the
        //  areaMesh
        //  Note: it is assumed that this is one patch
        mutable label areaPatchID_;

        //- areaShadowPatch index
        //  This is the patch on the fvMesh mesh that is directly opposite from
        //  the areaPatch
        //  Note: it is assumed that the fvMesh is one cell thick from the
        //  areaPatch to the areaShadowPatch
        mutable label areaShadowPatchID_;


    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
            const SolverPerformance<scalar>& solverPerfM,
            const SolverPerformance<scalar>& solverPerfw,
#else
            const lduSolverPerformance& solverPerfM,
            const lduSolverPerformance& solverPerfw,
#endif
            const areaScalarField& M,
            const areaScalarField& w
        );

        //- Return a reference to the areaPatch
        const fvPatch& areaPatch() const;

        //- Return a reference to the areaShadowPatch
        const fvPatch& areaShadowPatch() const;

        //- Calculate and set areaPatchID and shadowPatchID
        void calcAreaPatches() const;

        //- Map an area geometric field to a volume geometric field, where we
        //  assume the volume mesh to be one layer thick i.e. there is a
        //  one-to-one relationship between the number of area faces to the
        //  number of volume cells
        template<class Type>
        void mapAreaFieldToSingleLayerVolumeField
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Disallow default bitwise copy construct
        kirchhoffPlateSolid(const kirchhoffPlateSolid&);

        //- Disallow default bitwise assignment
        void operator=(const kirchhoffPlateSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("kirchhoffPlate");

    // Constructors

        //- Construct from components
        kirchhoffPlateSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~kirchhoffPlateSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                notImplemented(type());

                // Keep compiler happy
                return D();
            }

            //- Plate bending stiffness
            const dimensionedScalar& bendingStiffness() const
            {
                return bendingStiffness_;
            }

            //- Plate Poisson's ration
            const dimensionedScalar& nu() const
            {
                return nu_;
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
            ) const;

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label interfaceI,
                const label patchID,
                const vectorField& faceZoneTraction
            );

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "kirchhoffPlateSolidTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
