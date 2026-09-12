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
    thermalSolid

Description
    Heat-conduction-only solid model: the transient heat equation is solved for
    the temperature field T, using the run-time selectable thermal model for the
    specific heat capacity and thermal conductivity.

    No mechanical problem is solved: there is no momentum equation, no
    displacement, strain or stress solution, and tractionBoundarySnGrad()
    returns a zero field. The mechanical law is only consulted for the density
    used in the heat capacity term.

    This model is intended for the solid side of conjugate heat transfer, where
    it provides the functions required for explicit CHT coupling.

Author
    Philip Cardiff, UCD.
    Tukovic Zeljko, FSB.

SourceFiles
    thermalSolid.C

\*---------------------------------------------------------------------------*/

#ifndef thermalSolid_H
#define thermalSolid_H

#include "solidModel.H"
#include "thermalModel.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "SolverPerformance.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalSolid Declaration
\*---------------------------------------------------------------------------*/

class thermalSolid
:
    public solidModel
{
    // Private data

        //- Thermal model
        thermalModel thermal_;

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField k_;

        //- Temperature
        volScalarField T_;

        //- Temperature gradient
        volVectorField gradT_;

        //- Absolute temperature tolerance (in degrees)
        const scalar absTTol_;

    // Private Member Functions

        //- Check if the heat equation is converged
        bool converged
        (
            const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
            const SolverPerformance<scalar>& solverPerfT,
#else
            const lduSolverPerformance& solverPerfT,
#endif
            const volScalarField& T
        );

        //- Disallow default bitwise copy construct
        thermalSolid(const thermalSolid&);

        //- Disallow default bitwise assignment
        void operator=(const thermalSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }

public:

    //- Runtime type information
    TypeName("thermalSolid");

    // Constructors

        //- Construct from components
        thermalSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~thermalSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                return D();
            }

            //- Return const access to the cell temperature field
            const volScalarField& T() const
            {
                return T_;
            }

            //- Return non-const access to the cell temperature field
            volScalarField& T()
            {
                return T_;
            }

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Set temperature and heat flux at specified patch
            void setTemperatureAndHeatFlux
            (
                fvPatchScalarField& patch,
                const scalarField& temperature,
                const scalarField& heatFlux
            );

            //- Set temperature and heat flux at specified patch
            virtual void setTemperatureAndHeatFlux
            (
                const label interfaceI,
                const label patchID,
                const scalarField& faceZoneTemperature,
                const scalarField& faceZoneHeatFlux
            );

            //- Set equivalent interface heat transfer coeff
            void setEqInterHeatTransferCoeff
            (
                fvPatchScalarField& temperaturePatch,
                const scalarField& HTC
            );

            //- Set equivalent interface heat transfer coeff
            virtual void setEqInterHeatTransferCoeff
            (
                const label interfaceI,
                const label patchID,
                const scalarField& faceZoneHTC
            );

            //- Face zone temperature
            virtual tmp<scalarField> faceZoneTemperature
            (
                const label interfaceI
            ) const;

            //- Face zone heat flux
            virtual tmp<scalarField> faceZoneHeatFlux
            (
                const label interfaceI
            ) const;

            //- Face zone heat transfer coefficient
            virtual tmp<scalarField> faceZoneHeatTransferCoeff
            (
                const label interfaceI
            ) const;

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
