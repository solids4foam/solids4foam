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
    poroLinGeomSolid

Description
    Mathematical model for poro-elasticity, where linear geometry is assumed
    i.e. small strains and small rotations are assumed.

    The poro-pressure equation is solved.

    The stress is calculated by the run-time selectable mechanical law.

    The approach is based on the elastoPlasticBiotFoam in the miniGeotechFoam
    toolbox shared by Tian Tang, found at:
    https://bitbucket.org/tiantang/minigeotechfoam

    More reading at:

    T.Tang, O.Hededal & P.Cardiff (2014). On finite volume method implementation
    of poro-elasto-plasticity soil model. International Journal for
    Numerical and Analytical Methods in Geomechanics, 10.1002/nag.2361.

    T.Tang & O.Hededal (2014). Simulation of pore pressure accumulation under
    cyclic loading using finite volume method. Proceedings of the 8th European
    Conference on Numerical Methods in Geotechnical Engineering (NUMGE14),
    Volume 2, Pages 1301-1306.


Author
    Tian Tang, DTU.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    poroLinGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef poroLinGeomSolid_H
#define poroLinGeomSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class poroLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class poroLinGeomSolid
:
    public solidModel
{
    // Private data

        //- Displacement field
        volVectorField D_;

        //- Velocity field
        volVectorField U_;

        //- Point mesh
        pointMesh pMesh_;

        //- Point displacement field
        pointVectorField pointD_;

        //- Strain field
        volSymmTensorField epsilon_;

        //- Stress field
        volSymmTensorField sigma_;

        //- Cell-to-point interpolation
        newLeastSquaresVolPointInterpolation volToPoint_;

        //- Cell displacement gradient
        volTensorField gradD_;

        //- Density
        const volScalarField rho_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Poro pressure field
        volScalarField p_;

        //- Cell pressure gradient
        volVectorField gradp_;

        //- Hydraulic conductivity
        //  Note: Tang et al (2014) refers to this term as permeability
        dimensionedScalar hydraulicConductivity_;

        //- Specific weight of water
        //  Note: Tang et al (2014) refers to this term as water density
        dimensionedScalar gammaWater_;

        //- Porosity
        dimensionedScalar porosity_;

        //- Degree of saturation
        dimensionedScalar saturation_;

        //- Bulk modulus of water
        dimensionedScalar KWater_;

        //- Reciprocal of the effective bulk modulus of water and air mixture
        dimensionedScalar rKprime_;

        //- D equation under-relaxation factor
        const scalar DEqnRelaxFactor_;

        //- p equation under-relaxation factor
        const scalar pEqnRelaxFactor_;

        //- Solution standard tolerance
        const scalar solutionTol_;

        //- Solution tighter tolerance
        const scalar alternativeTol_;

        //- Material law tolerance
        const scalar materialTol_;

        //- Write frequency for residuals information
        const int infoFrequency_;

        //- Maximum number of momentum-pressure correctors
        const int nCorr_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField g_;

        //- Number of times the maximum number of correctors was reached in the
        //  momentum-pressure corrector loop
        int maxIterReached_;

    // Private Member Functions

        //- Check if the momentum and pressure equations are converged
        bool converged
        (
            const int iCorr,
            const lduMatrix::solverPerformance& solverPerfD,
            const lduMatrix::solverPerformance& solverPerfp
        );

        //- Disallow default bitwise copy construct
        poroLinGeomSolid(const poroLinGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const poroLinGeomSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("poroLinearGeometry");

    // Constructors

        //- Construct from components
        poroLinGeomSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~poroLinGeomSolid()
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

            //- Return const access to the cell pressure field
            const volScalarField& p() const
            {
                return p_;
            }

            //- Return non-const access to the cell pressure field
            volScalarField& p()
            {
                return p_;
            }

            //- Return point mesh
            const pointMesh& pMesh() const
            {
                return pMesh_;
            }

            //- Gradient of displacement
            const volTensorField& gradD() const
            {
                return gradD_;
            }

            //- Stress tensor
            const volSymmTensorField& sigma() const
            {
                return sigma_;
            }

            //- Strain tensor
            const volSymmTensorField& epsilon() const
            {
                return epsilon_;
            }

            //- Implicit stiffness
            const volScalarField& impK() const
            {
                return impK_;
            }

            //- Reciprocal of the implicit stiffness
            const volScalarField& rImpK() const
            {
                return rImpK_;
            }

            //- Return point displacement field
            virtual const pointVectorField& pointD() const
            {
                return pointD_;
            }

            //- Return multi-material interface pointer
            const newLeastSquaresVolPointInterpolation& volToPoint() const
            {
                return volToPoint_;
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
