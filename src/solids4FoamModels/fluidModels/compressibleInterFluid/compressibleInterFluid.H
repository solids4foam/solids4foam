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
    compressibleInterFluid

Description
    compressibleInterDyMFoam fluid model from OpenFOAM-v2012.

Author
    Ported by Philip Cardiff, UCD.

SourceFiles
    compressibleInterFluid.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAMESI

#ifndef compressibleInterFluid_H
#define compressibleInterFluid_H

#include "fluidModel.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"
//#include "CorrectPhi.H"
#include "CorrectPhiUpper.H"
#include "fvcSmooth.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class compressibleInterFluid Declaration
\*---------------------------------------------------------------------------*/

class compressibleInterFluid
:
    public fluidModel
{
    // Private data

        //- Cumulative continuity errors
        scalar cumulativeContErr_;

        //- Adjustable time step switch
        bool adjustTimeStep_;

        //- Max Courant number
        scalar maxCo_;

        //- Max time step
        scalar maxDeltaT_;

        //- Correct phi flag
        const bool correctPhi_;

        //- Check Courant number flag
        const bool checkMeshCourantNo_;

        // Move mesh in the outer corrector loop flag
        const bool moveMeshOuterCorrectors_;

        //- Local time step flag
        bool LTS_;

        //- Reciprocal of the local time-step field
        tmp<volScalarField> trDeltaT_;

        //- Pressure without rho*g*h component
        volScalarField p_rgh_;

        //- Transport model
        twoPhaseMixtureThermo mixture_;

        //- Alpha field references
        volScalarField& alpha1_;
        volScalarField& alpha2_;

        //- Density field references
        volScalarField& rho1_;
        volScalarField& rho2_;
        volScalarField rho_;

        //- Minium pressure
        dimensionedScalar pMin_;

        //- Gravity related fields
        uniformDimensionedScalarField hRef_;
        dimensionedScalar ghRef_;
        volScalarField gh_;
        surfaceScalarField ghf_;

        //- Mass flux
        surfaceScalarField rhoPhi_;

        //- dgdt field
        volScalarField dgdt_;

        //- alpha flux header
        IOobject alphaPhi10Header_;

        //- alpha restart flag
        const bool alphaRestart_;

        //- MULES flux from previous time-step
        surfaceScalarField alphaPhi10_;

        //- MULES Correction
        tmp<surfaceScalarField> talphaPhi1Corr0_;

        //- Turbulence model
        compressibleInterPhaseTransportModel turbulence_;

        //- Kinetic energy field
        volScalarField K_;

        //- MRF
        IOMRFZoneList MRF_;

        //- Pointer to the velocity surface field
        autoPtr<surfaceVectorField> Uf_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        compressibleInterFluid(const compressibleInterFluid&);

        //- Disallow default bitwise assignment
        void operator=(const compressibleInterFluid&);


public:

    //- Runtime type information
    TypeName("compressibleInterFluid");

    // Constructors

        //- Construct from components
        compressibleInterFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~compressibleInterFluid()
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

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the fluid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif

// ************************************************************************* //
