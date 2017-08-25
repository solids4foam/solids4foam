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

    More details in:
    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    nonLinGeomUpdatedLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomUpdatedLagSolid_H
#define nonLinGeomUpdatedLagSolid_H

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
                        Class nonLinGeomUpdatedLagSolid Declaration
\*---------------------------------------------------------------------------*/

class nonLinGeomUpdatedLagSolid
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

        //- Cell displacement increment gradient
        volTensorField gradDD_;

        //- Cell displacement gradient
        volTensorField gradD_;

        //- Total deformation gradient
        volTensorField F_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volTensorField relFinv_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volScalarField relJ_;

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


    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfD
        );

        //- Create global face zone point maps
        void findGlobalFaceZones();

        //- Disallow default bitwise copy construct
        nonLinGeomUpdatedLagSolid(const nonLinGeomUpdatedLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const nonLinGeomUpdatedLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::UPDATED_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("nonLinearGeometryUpdatedLagrangian");

    // Constructors

        //- Construct from components
        nonLinGeomUpdatedLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~nonLinGeomUpdatedLagSolid()
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
