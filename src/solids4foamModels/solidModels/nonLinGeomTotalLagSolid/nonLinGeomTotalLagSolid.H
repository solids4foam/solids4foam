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
    nonLinGeomTotalLagSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and total displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    nonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinGeomTotalLagSolid_H
#define nonLinGeomTotalLagSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class nonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class nonLinGeomTotalLagSolid
:
    public solidModel
{
    // Private data

        //- Total displacement field
        volVectorField D_;

        //- Velocity field
        volVectorField U_;

        //- Point mesh
        pointMesh pMesh_;

        //- Point displacement field
        pointVectorField pointD_;

        //- Cauchy (true) stress field
        volSymmTensorField sigma_;

        //- Cell-to-point interpolation
        volPointInterpolation volToPoint_;

        //- Cell displacement gradient
        volTensorField gradD_;

        //- Total deformation gradient
        volTensorField F_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Density (reference configuration)
        const volScalarField rho_;

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

        //- D equation under-relaxation factor
        const scalar DEqnRelaxFactor_;

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

        //- Number of times the maximum number of correctors was reached in the
        //  momentum equation
        int maxIterReached_;


    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
//            const lduMatrix::solverPerformance& solverPerfD
            const solverPerformance& solverPerfD
        );

        //- Disallow default bitwise copy construct
        nonLinGeomTotalLagSolid(const nonLinGeomTotalLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const nonLinGeomTotalLagSolid&);


public:

    //- Runtime type information
    TypeName("nonLinearGeometryTotalLagrangian");

    // Constructors

        //- Construct from components
        nonLinGeomTotalLagSolid(dynamicFvMesh& mesh);


    // Destructor

        virtual ~nonLinGeomTotalLagSolid()
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
            // const pointMesh& pMesh() const
            // {
            //     return pMesh_;
            // }

            //- Stress tensor
            const volSymmTensorField& sigma() const
            {
                return sigma_;
            }

            //- Strain tensor
            // const volSymmTensorField& epsilon() const
            // {
            //     return epsilon_;
            // }

            //- Return point displacement field
             virtual const pointVectorField& pointD() const
             {
                 return pointD_;
             }

            //- Return multi-material interface pointer
//             const leastSquaresVolPointInterpolation& volToPoint() const
            const volPointInterpolation& volToPoint() const
             {
                 return volToPoint_;
             }

            //- Return point displacement field
             virtual vector pointU(label pointID) const;

            //- Patch point displacement increment
             virtual tmp<vectorField> patchPointDisplacementIncrement
             (
                 const label patchID
             ) const;

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

            // //- Set pressure at specified patch
             virtual void setPressure
             (
                 const label patchID,
                 const scalarField& pressure
             );

            //- Set traction at specified patch
             virtual void setTraction
             (
                 const label patchID,
                 const label zoneID,
                 const vectorField& zoneTraction
             );

            //- Set pressure at specified patch
             virtual void setPressure
             (
                 const label patchID,
                 const label zoneID,
                 const scalarField& zonePressure
             );

            //- Predict traction at specified zone
             virtual tmp<vectorField> predictTraction
             (
                 const label patchID,
                 const label zoneID
             );

            //- Predict pressure at specified zone
             virtual tmp<scalarField> predictPressure
             (
                 const label patchID,
                 const label zoneID
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
