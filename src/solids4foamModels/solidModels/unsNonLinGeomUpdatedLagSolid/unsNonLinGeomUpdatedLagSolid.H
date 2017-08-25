/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    nonLinGeomUpdatedLagSolid

Description
    Mathematical model for finite strains and rotations, where an updated
    Lagrangian and incremental displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    "uns" in the class name refers to the "unstructured" and indicates that
    the face tangential gradient are calculated using a face-Gauss approach.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    unsNonLinGeomUpdatedLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsNonLinGeomUpdatedLagSolid_H
#define unsNonLinGeomUpdatedLagSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class unsNonLinGeomUpdatedLagSolid Declaration
\*---------------------------------------------------------------------------*/

class unsNonLinGeomUpdatedLagSolid
:
    public solidModel
{
    // Private data

        //- Displacement increment field
        volVectorField DD_;

        //- Total displacement field
        volVectorField D_;

        //- Velocity field
        volVectorField U_;

        //- Point mesh
        pointMesh pMesh_;

        //- Point displacement increment field
        pointVectorField pointDD_;

        //- Point total displacement field
        pointVectorField pointD_;

        //- Cauchy (true) stress field
        volSymmTensorField sigma_;

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Cell displacement increment gradient
        volTensorField gradDD_;

        //- Surface displacement increment gradient
        surfaceTensorField gradDDf_;

        //- Cell displacement gradient
        volTensorField gradD_;

        //- Total deformation gradient
        volTensorField F_;

        //- Total deformation gradient surface field
        surfaceTensorField Ff_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Jacobian of total deformation gradient surface field
        surfaceScalarField Jf_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Relative deformation gradient surface field
        surfaceTensorField relFf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volTensorField relFinv_;

        //- Relative Jacobian: surface field
        surfaceTensorField relFinvf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volScalarField relJ_;

        //- Relative Jacobian: surface field
        surfaceScalarField relJf_;

        //- Density (updated configuration)
        volScalarField rho_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        const surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Solution standard tolerance
        const scalar solutionTol_;

        //- Solution tighter tolerance
        const scalar alternativeTol_;

        //- Material law tolerance
        const scalar materialTol_;

        //- Write frequency for residuals information
        const int infoFrequency_;

        //- Maximum number of momentum correctors
        const int nCorr_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField g_;

        //- Number of times the maximum number of correctors was reached in the
        //  momentum equation
        int maxIterReached_;

        //- Stabilise hydrostatic pressure
        //  Apply smoothing to the relative Jacobian calculation
        const Switch stabilisePressure_;


    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfD
        );

        //- Move the mesh to the deformed configuration
        void moveMesh(const pointField& oldPoints);

        //- Create global face zone point maps
        void findGlobalFaceZones();

        //- Disallow default bitwise copy construct
        unsNonLinGeomUpdatedLagSolid(const unsNonLinGeomUpdatedLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const unsNonLinGeomUpdatedLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::UPDATED_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("unsNonLinearGeometryUpdatedLagrangian");

    // Constructors

        //- Construct from components
        unsNonLinGeomUpdatedLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~unsNonLinGeomUpdatedLagSolid()
        {}


    // Member Functions

        // Access

            //- Return const access to the cell displacement field
            const volVectorField& D() const
            {
                return D_;
            }

            //- Return non-const access to thecell displacement field
            volVectorField& D()
            {
                return D_;
            }

            //- Return point mesh
            const pointMesh& pMesh() const
            {
                return pMesh_;
            }

            //- Stress tensor
            const volSymmTensorField& sigma() const
            {
                return sigma_;
            }

            //- Return point displacement field
            virtual const pointVectorField& pointD() const
            {
                return pointD_;
            }

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                const label zoneID
            ) const;

            //- Face zone velocity surface (tangential) gradient
            virtual tmp<tensorField> faceZoneSurfaceGradientOfVelocity
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone current normal
            virtual tmp<vectorField> faceZoneNormal
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone primitive patch in current configuration
            virtual tmp<vectorField> currentFaceZonePoints
            (
                const label zoneID
            ) const;

        // Edit

            //- Set traction on the specified patch
            //  This function enables the FSI procedure to pass a traction from
            //  fluid to the solid
            virtual void setTraction
            (
                const label patchID,
                const vectorField& traction
            );

            //- Set pressure at specified patch
            virtual void setPressure
            (
                const label patchID,
                const scalarField& pressure
            );

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            //  Given the user-specified boundary patch traction, this function
            //  should return the surface normal gradient
            //  This function is used by traction-type boundary conditions to
            //  set the boundary gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- End of the program
            virtual void end();

            //- Return const acess to the solution (displacement field)
            virtual const volVectorField& solution() const
            {
                return D_;
            }

            //- Return non-const access to the solution (displacement field)
            virtual volVectorField& solution()
            {
                return D_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
