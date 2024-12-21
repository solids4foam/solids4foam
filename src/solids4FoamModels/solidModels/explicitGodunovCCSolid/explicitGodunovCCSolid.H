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

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit
    Philip Cardiff, UCD - solids4foam toolbox

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project. 
    Additional modifications include the integration of a dual-time stepping 
    algorithm. 

    Note: The angular momentum algorithm has not been adapted to the new 
    time discretization scheme. As a result, it may not provide accurate 
    solutions if used in its current form. 

Class
    explicitGodunovCCSolid

Description
    A solid model based on the explicitSolidDynamics toolkit.
    Solves a first-order system of total Lagrangian mixed formulation
    comprising conservation laws for linear momentum and deformation
    gradient.

SourceFiles
    explicitGodunovCCSolid.C

\*---------------------------------------------------------------------------*/


#ifndef explicitGodunovCCSolid_H
#define explicitGodunovCCSolid_H

#include "solidModel.H"
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
                        Class explicitGodunovCCSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitGodunovCCSolid
:
    public solidModel
{
    // Private data
        Time& runTime_;

    // Reading dictionaries 
    // Mechanical properties
    IOdictionary mechanicalProperties_;

    // Control dictionary
    IOdictionary controlDict_;

    scalar incompressibilityCoefficient_; 
    // Stabilisation parameter for near incompressibility
    const scalar& beta_;

    // Angular momentum conservation
    const word angularMomentumConservation_;

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
  
    //Creating mesh coordinate fields

    // Material cell center coordinates
    const volVectorField& C_;

    // Spatial cell center coordinates
    volVectorField x_;
    // material cell center coordinates
    const volVectorField X_;

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
    // Time increment (pesudo)
    dimensionedScalar pDeltaT_;

    // Runge-Kutta stage
    scalarList RKstages_;

    //limiter fields
    volVectorField phi_lm_;
    volTensorField phi_P_;

    // Private Member Functions
    //- Check for solution convergence
    bool converged
    (
        const int iCorr,
        const dimensionedScalar pDeltaT,
        const GeometricField<vector, fvPatchField, volMesh>& vf
    );

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

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
