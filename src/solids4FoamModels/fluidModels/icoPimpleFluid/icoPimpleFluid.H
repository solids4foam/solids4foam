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
    icoPimpleFluid

Description
    icoPimpleFoam fluid model: Same as icoFluid but with support for the "elastic"
    Robin boundary conditions for partitioned FSI as described in:

        Tuković et al. Added Mass Partitioned Fluid–Structure Interaction
        Solver Based on a Robin Boundary Condition for Pressure. OpenFOAM.
        2019;1–22.

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    icoPimpleFluid.C

\*---------------------------------------------------------------------------*/

#ifndef icoPimpleFluid_H
#define icoPimpleFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class icoPimpleFluid Declaration
\*---------------------------------------------------------------------------*/

class icoPimpleFluid
:
    public fluidModel
{
    // Private data

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence_;

        //- Density
        dimensionedScalar rho_;

        //- Acceleration field
        volVectorField ddtU_;

        //- Viscous laplacian field
        volVectorField laplNuU_;
    
        //- Velocity field
        surfaceVectorField Uf_;

        //- rAUf field
        volScalarField rAU_;
    
        //- rAUf field
        surfaceScalarField rAUf_;

        //- Pressure time series
        interpolationTable<scalar> refPressureSeries_;

        //- Apply temporal consistency correction
        Switch ddtCorr_;
    
        //- Apply new procedure for phi adjustment
        Switch adjustPhi_;
    
        //- Solve energy equation
        Switch solveEnergyEq_;
    
        //- Temperature field
        autoPtr<volScalarField> TPtr_;

        //- Thermal conductivity field
        autoPtr<volScalarField> lambdaEffPtr_;

    // Private Member Functions

        //- Reset flux if mesh changed
        void updateRobinFsiInterfaceFlux();

        //- Update phi for Robin FSI boundary conditions
        void updateRobinFsiInterface();

        //- Correct phi after solving pressure equation
        void correctRobinFsiInterfaceFlux();

        //- Solve energy equation if requred
        void solveEnergyEq();

        //- Return Boussinesq momentu source term
        tmp<volVectorField> boussinesqMomentumSource() const;
    
        //- Disallow default bitwise copy construct
        icoPimpleFluid(const icoPimpleFluid&);

        //- Disallow default bitwise assignment
        void operator=(const icoPimpleFluid&);


public:

    //- Runtime type information
    TypeName("icoPimpleFluid");


    // Constructors

        //- Construct from components
        icoPimpleFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~icoPimpleFluid()
        {}


    // Member Functions

        // Access

            //- Density
            const dimensionedScalar& rho()
            {
                return rho_;
            }

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

        // Edit
			
            //- Evolve the fluid model
            void evolveConsistent();

            //- Evolve the fluid model
            void evolveInconsistent();

            //- Evolve the fluid model
            virtual bool evolve();

            // //- Set temperature at specified patch
            // virtual void setTemperature
            // (
            //     const label interfaceI,
            //     const label patchID,
            //     const scalarField& faceZoneTemperature
            // );

            //- Set temperature and heat flux at specified patch
            virtual void setTemperatureAndHeatFlux
            (
                fvPatchScalarField& patch,
                const scalarField& temperature,
                const scalarField& heatFlux
            );

            //- Set temperature and heat flux at specified patch
            virtual void setTemperatureAndHeatFlux
            (
                const label interfaceI,
                const label patchID,
                const scalarField& faceZoneTemperature,
                const scalarField& faceZoneHeatFlux
            );

            //- Set temperature and heat flux at specified patch
            virtual void setEqInterHeatTransferCoeff
            (
                fvPatchScalarField& patch,
                const scalarField& heatTransferCoeff
            );
    
            //- Set eq interface heat transfer coeff at specified patch
            virtual void setEqInterHeatTransferCoeff
            (
                const label interfaceI,
                const label patchID,
                const scalarField& faceZoneHeatTransferCoeff
            );

            //- Patch temperature
            virtual tmp<scalarField> patchTemperature
            (
                const label patchID
            ) const;
    
            //- Patch heat flux (W/m2)
            virtual tmp<scalarField> patchHeatFlux(const label patchID) const;
    
            //- Patch heat transfer coefficient
            virtual tmp<scalarField> patchHeatTransferCoeff
            (
                const label patchID
            ) const;

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
