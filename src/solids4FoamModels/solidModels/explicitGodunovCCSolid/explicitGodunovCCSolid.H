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
    explicitGodunovCCSolid

Description
    A solid mechanics solver based on a Total Lagrangian mixed formulation
    comprising of conservation laws for linear momentum and deformation
    gradient of the system.

Original Author


Integrated with solids4Foam
    Khoder Alhamwi Alshaar 

SourceFiles
    explicitGodunovCCSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitGodunovCC_H
#define explicitGodunovCC_H

#include "solidModel.H"
// #include "volFields.H"
// #include "surfaceFields.H"
// #include "pointFields.H"
// #include "uniformDimensionedFields.H"
// #include "mechanicalEnergies.H"

#include "operations.H"
#include "solidMaterialModel.H"
#include "mechanics.H"
#include "gradientSchemes.H"
#include "interpolationSchemes.H"
#include "angularMomentum.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitGodunovCCSolid
:
    public solidModel
{
    // Private data
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

    //-----------------------------------

    Time& runTime_;

    // Reading dictionaries 
    // Mechanical properties
    IOdictionary mechanicalProperties_;

    // Control dictionary
    IOdictionary controlDict_;

    // Finite volume solution dictionary
    IOdictionary fvSolution_;

    // Stabilisation parameter for near incompressibility
    const scalar& beta_;

    // Angular momentum conservation
    const word angularMomentumConservation_;

    // Optional linear momentum damping coefficient
    const dimensionedScalar dampingCoeff_;

    //Creating mesh parameters 

    // Mesh-related fields
    operations op_;

    // Point mesh
    // pointMesh pMesh(mesh); // available pMesh() 

    // Material face area
    const surfaceScalarField& magSf_;

    // Material face area normal vector
    const surfaceVectorField& Sf_;
    
    // Minimum edge length
    const dimensionedScalar h_;

    // Boundary patches
    const polyBoundaryMesh& bm_;
    const label& symmetricPatchID_;
    const label& symmetricXpatchID_;
    const label& symmetricYpatchID_;
    const label& symmetricZpatchID_;
   
    //Creating mesh coordinate fields

    // Material cell center coordinates
    const volVectorField& C_;

    // Spatial cell center coordinates
    volVectorField x_;

    // Spatial nodal coordinates
    pointVectorField xN_;

    // Material nodal coordinates
    pointVectorField XN_;   

    // Spatial face center coordinates
    surfaceVectorField xF_;

    //Creating mesh normal fields

    // Material normals
    const surfaceVectorField N_;

    // Spatial normals
    surfaceVectorField n_;

    //Creating linear momentum fields
    // Cell linear momentum
    volVectorField lm_;

    // Nodal linear momentum
    pointVectorField lmN_;

    //Creating strain measure fields 

    // Deformation gradient tensor
    volTensorField F_;

    // Cofactor of deformation
    volTensorField H_;

    // Jacobian of deformation
    volScalarField J_;

    // Creating constitutive model

    // Solid model class
    solidMaterialModel model_;

    // Density
    const dimensionedScalar& rho_;

    // Pressure
    volScalarField p_;

    // First Piola Kirchhoff stress tensor
    volTensorField P_;
    volVectorField Px_;
    volVectorField Py_;
    volVectorField Pz_;

    // Creating fields for wave speeds

    // Continuum mechanics class
    mechanics mech_;

    // Longitudinal wave speed
    volScalarField Up_;

    // Wave speed for time increment
    volScalarField Up_time_;

    // Shear wave speed
    volScalarField Us_;

    // Creating fields for gradient  
    // Gradient class
    gradientSchemes grad_;

    // Gradient of cell linear momentum
    volTensorField lmGrad_;

    // Gradients of first Piola Kirchhoff stress tensor
    volTensorField PxGrad_;
    volTensorField PyGrad_;
    volTensorField PzGrad_;

    //Creating fields for reconstruction 
    // Reconstruction of linear momentum
    surfaceVectorField lm_M_;
    surfaceVectorField lm_P_;

    // Reconstruction of PK1 stresses
    surfaceTensorField P_M_;
    surfaceTensorField P_P_;

    // Reconstruction of traction
    surfaceVectorField t_M_;
    surfaceVectorField t_P_;

    // Creating fields for riemann solver
    // Surface tensors
    surfaceTensorField S_lm_;
    surfaceTensorField S_t_;

    // Contact traction
    surfaceVectorField tC_;

    // Contact linear momentum
    surfaceVectorField lmC_;

    // Volumetric fields
    volVectorField t_b_;
    volVectorField lm_b_;

    // Creating fields for the constrained procedure
    // Constrained class
    interpolationSchemes interpolate_;

    // Cell-averaged linear momentum
    volVectorField lmR_;

    // Local gradient of cell-averaged linear momentum
    volTensorField lmRgrad_;


    // Creating fields for angular momentum
    // Angular momentum class
    angularMomentum am_;

    // RHS of linear momentum equation
    volVectorField rhsLm_;
    volVectorField rhsLm1_;

    // RHS of angular momentum equation
    volVectorField rhsAm_;


    // Creating fields for post-processing
    // Nodal displacement field
    pointVectorField uN_;

    //Creating variables for time 
    // Time increment
    dimensionedScalar deltaT_;

    // Runge-Kutta stage
    scalarList RKstages_;



    // Private Member Functions
    
        //- Update the stress field
        void updateStress();

        //- Disallow default bitwise copy construct
        explicitGodunovCCSolid(const explicitGodunovCCSolid&);

        //- Disallow default bitwise assignment
        void operator=(const explicitGodunovCCSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("explicitGodunovCC");

    // Constructors

        //- Construct from components
        explicitGodunovCCSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~explicitGodunovCCSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
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
            ) const;

            //- End of the program
            virtual void end()
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
