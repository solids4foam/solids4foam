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
    interFluid

Description
    interFoam fluid model based on interDyMFoam from foam-extend-4.1.

Author
    Based on interFoam.
    Ported by Philip Cardiff, UCD.

SourceFiles
    interFluid.C

\*---------------------------------------------------------------------------*/

#ifndef interFluid_H
#define interFluid_H

#include "fluidModel.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class interFluid Declaration
\*---------------------------------------------------------------------------*/

class interFluid
:
    public fluidModel
{
    // Private data

        //- pd
        volScalarField pd_;

        //- Volume fraction of phase 1
        volScalarField alpha1_;

        //- Transport model
        twoPhaseMixture twoPhaseProperties_;

        //- Density of phase 1
        const dimensionedScalar& rho1_;

        //- Density of phase 2
        const dimensionedScalar& rho2_;

        //- Average density field
        volScalarField rho_;

        //- Mass flux
        surfaceScalarField rhoPhi_;

        //- Interface
        interfaceProperties interface_;

        //- pd reference cell index
        label pdRefCell_;

        //- pd reference value
        scalar pdRefValue_;

        //- p reference value
        scalar pRefValue_;

        //- Correct phi flag
        const bool correctPhi_;

        //- Check Courant number flag
        const bool checkMeshCourantNo_;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence_;

        //- Pressure correction boundary condition types
        //  Needed by correctPhi.H
        wordList pcorrTypes_;

        //- Cumulative continuity errors
        scalar sumLocalContErr_;
        scalar globalContErr_;
        scalar cumulativeContErr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interFluid(const interFluid&);

        //- Disallow default bitwise assignment
        void operator=(const interFluid&);


public:

    //- Runtime type information
    TypeName("interFluid");

    // Constructors

        //- Construct from components
        interFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~interFluid()
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
