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
    sonicLiquidFluid

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible liquid. Assumes a linear equation of state between 
	density and pressure related by the bulk modulus.

Author
    Iago Lessa de Oliveira, FEIS/UNESP/UCD.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    sonicLiquidFluid.C

\*---------------------------------------------------------------------------*/

#ifndef sonicLiquidFluid_H
#define sonicLiquidFluid_H

#include "fluidModel.H"
#include "pimpleControl.H"
#include "uniformDimensionedFields.H"
#ifdef OPENFOAMESIORFOUNDATION
	#include "pressureControl.H"
	#include "fvOptions.H"
	#include "fvcSmooth.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class sonicLiquidFluid Declaration
\*---------------------------------------------------------------------------*/

class sonicLiquidFluid
:
    public fluidModel
{
    // Private data

        // Transport properties dictionary
        IOdictionary transportProperties_;

        //- Dynamic viscosity
        dimensionedScalar mu_;

		//- Thermodynamic properties
		IOdictionary thermodynamicProperties_;

		//- Ref. fluid density
		const dimensionedScalar rho0_;

		//- Ref. pressure
		const dimensionedScalar p0_;

        //- Fluid bulk modulus field
        const dimensionedScalar K_;

        //- const Ref. of psi field
        const dimensionedScalar psi_;

		//- Density offset
        const volScalarField rhoO_;

        //- Fluid density
        volScalarField rho_;

#ifdef OPENFOAMESIORFOUNDATION
		//- Creates and initialises the momentum field rhoUf
		autoPtr<surfaceVectorField> rhoUf_;

		//- Finite volume options: source terms, constrains and corrections
		fv::options& fvOptions_;

#endif
        //- rAU field
        volScalarField rAU_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

    // Private Member Functions

        //- Calculates and prints the continuity errors
        void compressibleContinuityErrs();

        //- Solve density equation
        void solveRhoEqn();

        //- Phi correction for mesh motion
        void CorrectFlux();

        //- Print compressible Courant number
		void compressibleCourantNo();

        //- Solve pressure equation
		void solvePEqn
		(
			const fvVectorMatrix& UEqn
		);

        //- Disallow default bitwise copy construct
        sonicLiquidFluid(const sonicLiquidFluid&);

        //- Disallow default bitwise assignment
        void operator=(const sonicLiquidFluid&);


public:

    //- Runtime type information
    TypeName("sonicLiquidFluid");

    // Constructors

        //- Construct from components
        sonicLiquidFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~sonicLiquidFluid()
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
