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
    pimpleFluid

Description
    Same as pimpleDyMFoam in foam-extend-4.1 but with an improve ddtPhi
    correction and support for the "elastic" Robin boundary conditions for
    partitioned FSI as described in:

        Tuković et al. Added Mass Partitioned Fluid–Structure Interaction
        Solver Based on a Robin Boundary Condition for Pressure. OpenFOAM.
        2019;1–22.

    In addition, thermal functions are implemented to allow conjugate heat
    transfer.

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    pimpleFluid.foamextend.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleFluid_H
#define pimpleFluid_H

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
                     Class pimpleFluid Declaration
\*---------------------------------------------------------------------------*/

class pimpleFluid
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

        //- rAU field
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
        pimpleFluid(const pimpleFluid&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleFluid&);


public:

    //- Runtime type information
    TypeName("pimpleFluid");


    // Constructors

        //- Construct from components
        pimpleFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~pimpleFluid()
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
