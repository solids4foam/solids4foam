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
    pimpleFoam fluid model from OpenFOAM-9.

Author
    Ported by Philip Cardiff, UCD.

SourceFiles
    pimpleFluid.foamextend.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleFluid_H
#define pimpleFluid_H

#include "fluidModel.H"
#include "localEulerDdtScheme.H"
#include "pressureReference.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"

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

        //- Local time step flag
        bool LTS_;

        //- Reciprocal of the local time-step field
        tmp<volScalarField> trDeltaT_;

        //- Pointer to the velocity surface field
        autoPtr<surfaceVectorField> Uf_;

        //- Pressure reference object
        pressureReference pressureReference_;

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
        autoPtr<incompressible::momentumTransportModel> turbulence_;

        //- Fluid density [kg/m3]
        dimensionedScalar rho_;

        //- Correct phi flag
        const bool correctPhi_;

        //- Check Courant number flag
        const bool checkMeshCourantNo_;

        // Move mesh in the outer corrector loop flag
        const bool moveMeshOuterCorrectors_;

        //- Cumulative continuity errors
        scalar cumulativeContErr_;


    // Private Member Functions

        //- Correct phi function
        void correctPhi();

        //- Solve pressure equation
        void solvePEqn
        (
            const scalar& UUrf,
            const fvVectorMatrix& ddtUEqn,
            const fvVectorMatrix& HUEqn
        );

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
