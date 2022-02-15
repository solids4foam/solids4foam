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
    pisoFluid

Description
    pisoFoam fluid model

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    pisoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef pisoFluid_H
#define pisoFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#ifdef OPENFOAMESI
    #include "turbulentTransportModel.H"
    #include "fvOptions.H"
#elif OPENFOAMFOUNDATION
    #include "kinematicMomentumTransportModel.H"
#else
    #include "turbulenceModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class pisoFluid Declaration
\*---------------------------------------------------------------------------*/

class pisoFluid
:
    public fluidModel
{
    // Private data

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
#ifdef OPENFOAMFOUNDATION
        autoPtr<incompressible::momentumTransportModel> turbulence_;
#else
        autoPtr<incompressible::turbulenceModel> turbulence_;
#endif

        //- Density
        dimensionedScalar rho_;

#ifdef OPENFOAMESI
        //- Create finite volume options
        fv::options& fvOptions_;
#endif

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pisoFluid(const pisoFluid&);

        //- Disallow default bitwise assignment
        void operator=(const pisoFluid&);


public:

    //- Runtime type information
    TypeName("pisoFluid");

    // Constructors

        //- Construct from components
        pisoFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~pisoFluid()
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
