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
    thermalLinGeomSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the heat equation is solved.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    thermalLinGeomSolid.C

\*---------------------------------------------------------------------------*/

#ifndef thermalLinGeomSolid_H
#define thermalLinGeomSolid_H

#include "linGeomSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalLinGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class thermalLinGeomSolid
:
    public linGeomSolid
{
    // Private data

        //- Temperature field
        volScalarField T_;

        //- Referemce stress-free temperature field
        volScalarField T0_;

        //- Cell temperature gradient
        volVectorField gradT_;

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField k_;

        //- T equation under-relaxation factor
        const scalar TEqnRelaxFactor_;

        //- Solution tolerance for outer non-orthogonality correctors
        const scalar solutionTol_;

        //- Alternative tolerance for outer non-orthogonality correctors
        const scalar alternativeTol_;

        //- Write frequency for residuals information
        const int infoFrequency_;

        //- Maximum number of thermal correctors
        const int nCorr_;

        //- Number of times the maximum number of correctors was reached
        int maxIterReached_;

    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduMatrix::solverPerformance& solverPerfD
        );

        //- Disallow default bitwise copy construct
        thermalLinGeomSolid(const thermalLinGeomSolid&);

        //- Disallow default bitwise assignment
        void operator=(const thermalLinGeomSolid&);


public:

    //- Runtime type information
    TypeName("thermalLinearGeometry");

    // Constructors

        //- Construct from components
        thermalLinGeomSolid(dynamicFvMesh& mesh);


    // Destructor

        virtual ~thermalLinGeomSolid()
        {}


    // Member Functions

        // Access

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

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            //  Given the user-specified boundary patch traction, this function
            //  should return the surface normal gradient
            //  This function is used by traction-type boundary conditions to
            //  set the boundary gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- End of the program
            virtual void end();

            //- Return const acess to the solution (temperature field)
            // virtual const volVectorField& solution() const
            // {
            //     return T_;
            // }

            //- Return non-const access to the solution (temperature field)
            // virtual volVectorField& solution()
            // {
            //     return T_;
            // }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
