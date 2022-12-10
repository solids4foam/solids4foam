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
    interFluid

Description
    interFoam fluid model based on interFoam from OpenFOAM-v2012.

Author
    Based on interFoam.
    Ported by Philip Cardiff, UCD.

SourceFiles
    interFluid.C

\*---------------------------------------------------------------------------*/

#ifndef interFluid_H
#define interFluid_H

#include "fluidModel.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
// #include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

/*---------------------------------------------------------------------------*\
                        Class interFluid Declaration
\*---------------------------------------------------------------------------*/

class interFluid
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

        //- Transport model
        immiscibleIncompressibleTwoPhaseMixture mixture_;

        //- Average density field
        volScalarField rho_;

        //- Mass flux
        surfaceScalarField rhoPhi_;

        //- Incompressible turbulence model
        autoPtr<transportModelType> turbulence_;

        //- Pressure without rho*g*h component
        volScalarField p_rgh_;

        //- Gravity related fields
        uniformDimensionedScalarField hRef_;
        dimensionedScalar ghRef_;
        volScalarField gh_;
        surfaceScalarField ghf_;

        //- p reference cell index
        label pRefCell_;

        //- p reference value
        scalar pRefValue_;

        //- Reciprocal of the momentum central coefficient
        tmp<volScalarField> rAU_;

        //- Correct phi flag
        const bool correctPhi_;

        //- Check Courant number flag
        const bool checkMeshCourantNo_;

        // Move mesh in the outer corrector loop flag
        const bool moveMeshOuterCorrectors_;

        //- Cumulative continuity errors
        scalar cumulativeContErr_;

        //- MULES compressed flux is registered in case scalarTransport FO needs
        //  it
        surfaceScalarField alphaPhiUn_;

        //- alpha flux header
        IOobject alphaPhi10Header_;

        //- alpha restart flag
        const bool alphaRestart_;

        // MULES flux from previous time-step
        surfaceScalarField alphaPhi10_;

        // MULES Correction
        tmp<surfaceScalarField> talphaPhi1Corr0_;


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
