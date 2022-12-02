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
    pimpleOversetFluid

Description
    Transient solver for incompressible, flow of Newtonian fluids
    with dynamic mesh and overset mesh using the PIMPLE
    (merged PISO-SIMPLE) algorithm.
    Currently this is inconsistent version of the solver.


Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    pimpleOversetFluid.C

\*---------------------------------------------------------------------------*/

#if FOAMEXTEND

#ifndef pimpleOversetFluid_H
#define pimpleOversetFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class pimpleOversetFluid Declaration
\*---------------------------------------------------------------------------*/

class pimpleOversetFluid:    public fluidModel{
    // Private data
        singlePhaseTransportModel laminarTransport_;  //- Transport model
        autoPtr<incompressible::turbulenceModel> turbulence_;      //- Turbulence model
       
        dimensionedScalar rho_; //- Density
        label pRefCell_;  //- Reference pressure cell
        scalar pRefValue_;  //- Reference pressure value

    // Private Member Functions
        pimpleOversetFluid(const pimpleOversetFluid&);  //- Disallow default bitwise copy construct
        void operator=(const pimpleOversetFluid&);//- Disallow default bitwise assignment

public:

    TypeName("pimpleOversetFluid");  //- Runtime type information
    // Constructors
        //- Construct from components
        pimpleOversetFluid (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );
    // Destructor
        virtual ~pimpleOversetFluid(){}

    // Member Functions
        // Access
            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce  (const label patchID ) const;
            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce  (const label patchID) const;
            //-Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff(const label zoneID,const label patchID) const;

        // Edit
            //- Evolve the fluid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif // FOAMEXTEND

// ************************************************************************* //
