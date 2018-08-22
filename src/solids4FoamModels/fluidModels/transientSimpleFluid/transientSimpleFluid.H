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
    transientSimpleFluid

Description
    transientSimpleFoam (pimple) fluid model

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    transientSimpleFluid.C

\*---------------------------------------------------------------------------*/

#ifndef transientSimpleFluid_H
#define transientSimpleFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class transientSimpleFluid Declaration
\*---------------------------------------------------------------------------*/

class transientSimpleFluid
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


    // Private Member Functions

        //- Disallow default bitwise copy construct
        transientSimpleFluid(const transientSimpleFluid&);

        //- Disallow default bitwise assignment
        void operator=(const transientSimpleFluid&);

public:

    //- Runtime type information
    TypeName("transientSimpleFluid");

    // Constructors

        //- Construct from components
        transientSimpleFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~transientSimpleFluid()
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

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
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
