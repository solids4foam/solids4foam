/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    pimpleDyMFoam fluid model from foam-extend-4.1 with support for the
    "elastic" Robin boundary conditions for partitioned FSI as described in:

        Tuković et al. Added Mass Partitioned Fluid–Structure Interaction
        Solver Based on a Robin Boundary Condition for Pressure. OpenFOAM.
        2019;1–22.

Author
    Hrvoje Jasak, Wikki Ltd. All rights reserved.
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

        //- Fluid density [kg/m3]
        dimensionedScalar rho_;

        //- Acceleration field
        volVectorField ddtU_;

        //- Velocity field
        surfaceVectorField Uf_;

        //- rAU field
        volScalarField rAU_;

        //- rAUf field
        surfaceScalarField rAUf_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- Correct phi flag
        const bool correctPhi_;

        //- Check Courant number flag
        const bool checkMeshCourantNo_;

        //- Cumulative continuity errors
        scalar sumLocalContErr_;
        scalar globalContErr_;
        scalar cumulativeContErr_;


    // Private Member Functions

        //- Reset flux if mesh changed
        void updateRobinFsiInterfaceFlux();

        //- Update phi for Robin FSI boundary conditions
        void updateRobinFsiInterface();

        //- Correct phi after solving pressure equation
        void correctRobinFsiInterfaceFlux();

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
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
