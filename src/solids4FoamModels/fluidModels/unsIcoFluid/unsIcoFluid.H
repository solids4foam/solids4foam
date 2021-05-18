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
    unsIcoFluid

Description
    unsIcoFoam fluid model: Same as icoFluid but with support for the "elastic"
    Robin boundary conditions for partitioned FSI as described in:

        Tuković et al. Added Mass Partitioned Fluid–Structure Interaction
        Solver Based on a Robin Boundary Condition for Pressure. OpenFOAM.
        2019;1–22.

Author
    Hrvoje Jasak, Wikki Ltd. All rights reserved.
    Zeljko Tukovic, FSB Zagreb. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    unsIcoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef unsIcoFluid_H
#define unsIcoFluid_H

#include "fluidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class unsIcoFluid Declaration
\*---------------------------------------------------------------------------*/

class unsIcoFluid
:
    public fluidModel
{
    // Private data

        // Transport properties dictionary
        IOdictionary transportProperties_;

        //- Kinematic viscosity
        dimensionedScalar nu_;

        //- Density
        dimensionedScalar rho_;

		//- Acceleration field
		volVectorField ddtU_;

		//- Velocity field
		surfaceVectorField Uf_;

        //- rAUf field
        surfaceScalarField rAUf_;

    // Private Member Functions

        //- Reset flux if mesh changed
        void updateRobinFsiInterfaceFlux();

        //- Update phi for Robin FSI boundary conditions
        void updateRobinFsiInterface();

        //- Correct phi after solving pressure equation
        void correctRobinFsiInterfaceFlux();

        //- Disallow default bitwise copy construct
        unsIcoFluid(const unsIcoFluid&);

        //- Disallow default bitwise assignment
        void operator=(const unsIcoFluid&);


public:

    //- Runtime type information
    TypeName("unsIcoFluid");


    // Constructors

        //- Construct from components
        unsIcoFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~unsIcoFluid()
        {}


    // Member Functions

        // Access

            //- Return kinematic viscosity
            const dimensionedScalar& nu() const
            {
                return nu_;
            }

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
