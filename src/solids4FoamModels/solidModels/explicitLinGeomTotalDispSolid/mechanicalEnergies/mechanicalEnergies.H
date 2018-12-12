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
    mechanicalEnergies

Description
    Class used to calculate the energies involved in a solid mechanics
    simulation, including the external work, internal strain energy, and
    kinetic energy.

    The class also calculates a linear bulk viscous pressure term, which can
    be used for stablisation in explicit solution algorithms.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    mechanicalEnergies.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalEnergies_H
#define mechanicalEnergies_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class mechanicalEnergies Declaration
\*---------------------------------------------------------------------------*/

class mechanicalEnergies
{
    // Private data

        // Reference to mesh
        const fvMesh& mesh_;

        //- Linear bulk viscosity coefficient
        //  Defaults to 0.06, as in Abaqus
        const scalar linearBulkViscosityCoeff_;

        //- Bulk viscosity stablisation pressure
        autoPtr<surfaceScalarField> viscousPressurePtr_;

        //- Volumetric strain
        //  This is used to calculate the linear bulk damping pressure
        autoPtr<volScalarField> epsilonVolPtr_;

        //- External work field used to store the integral of the boundary work
        volScalarField externalWorkField_;

        //- Internal energy field used to store the integral of the internal
        //  energy
        volScalarField internalEnergyField_;

        //- Energy dissipated due to the bulk viscosity stablisation term
        volScalarField bulkViscosityEnergyField_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mechanicalEnergies(const mechanicalEnergies&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalEnergies&);

public:

    //- Runtime type information
    TypeName("mechanicalEnergies");

    // Constructors

        //- Construct from components
        mechanicalEnergies
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        virtual ~mechanicalEnergies()
        {}


    // Member Functions

        //- Viscous pressure
        const surfaceScalarField& viscousPressure
        (
            const volScalarField& rho,
            const surfaceScalarField& waveSpeed,
            const volTensorField& gradD
        );

        //- Volumetric strain
        const volScalarField& epsilonVol
        (
            const volTensorField& gradD
        );

        //- Check global energy conservation
        void checkEnergies
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& DD,
            const volSymmTensorField& sigma,
            const volTensorField& gradD,
            const volTensorField& gradDD,
            const surfaceScalarField& waveSpeed
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
