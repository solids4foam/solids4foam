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
    buoyantPimpleFluid

Description
    Transient solver for buoyant, turbulent flow of compressible fluids for
    ventilation and heat-transfer.

    Turbulence is modelled using a run-time selectable compressible RAS or
    LES model.

Author
    Danial Khazaei, NIT. All rights reserved.

SourceFiles
    buoyantPimpleFluid.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAMESIORFOUNDATION

#ifndef buoyantPimpleFluid_H
#define buoyantPimpleFluid_H

#include "fluidModel.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#ifdef OPENFOAMFOUNDATION
    #include "radiationModel.H"
#endif
#include "fvOptions.H"
#include "pimpleControl.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class buoyantPimpleFluid Declaration
\*---------------------------------------------------------------------------*/

class buoyantPimpleFluid
:
    public fluidModel
{
    // Private data

        //- Pressure field
        volScalarField p_rgh_;

        //- Thermophysical model
        autoPtr<rhoThermo> pThermo_;

        //- Ref. of thermophysical model
        rhoThermo& thermo_;

        //- Fluid density
        volScalarField rho_;

        //- const Ref. of psi field
        const volScalarField& psi_;

        //- dpdt field
        volScalarField dpdt_;

        //- Kinetic energy field
        volScalarField K_;

        //- Gravitional acceleration
        const uniformDimensionedVectorField g_;

        //- Reference height field
        uniformDimensionedScalarField hRef_;

        //- gravity*reference height
        dimensionedScalar ghRef_;

        //- gravity*height
        volScalarField gh_;

        //- gravity*height: face-interpolated
        surfaceScalarField ghf_;

        //- Finite volume options: source terms, constrains and corrections
        fv::options& fvOptions_;

        //- Turbulent Transport model
        autoPtr<compressible::turbulenceModel> turbulence_;

#ifdef OPENFOAMFOUNDATION
        //- Create radiation model
        autoPtr<radiationModel> radiation_;
#endif

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;


    // Private Member Functions

        //- Solve density equation
        void solveRhoEqn();

        //- Solve velocity equation
        tmp<fvVectorMatrix> solveUEqn();

        //- Solve energy equation
        void solveEEqn();

        //- Solve pressure equation
        void solvePEqn
        (
            const fvVectorMatrix& UEqn
        );

        //- Disallow default bitwise copy construct
        buoyantPimpleFluid(const buoyantPimpleFluid&);

        //- Disallow default bitwise assignment
        void operator=(const buoyantPimpleFluid&);


public:

    //- Runtime type information
    TypeName("buoyantPimpleFluid");

    // Constructors

        //- Construct from components
        buoyantPimpleFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~buoyantPimpleFluid()
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

#endif // OPENFOAMESIORFOUNDATION

// ************************************************************************* //
