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
    AitkenCouplingInterface

Description
    Strong Dirichlet-Neumann coupling with Aitken accelerated under-relaxation.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    AitkenCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef AitkenCouplingInterface_H
#define AitkenCouplingInterface_H

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class AitkenCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class AitkenCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- Fixed under-relaxation factor
        const scalar relaxationFactor_;

        //- Predict solid
        const Switch predictSolid_;

        //- List of Aitken under-relaxation factors
        List<scalar> aitkenRelaxationFactors_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        AitkenCouplingInterface(const AitkenCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const AitkenCouplingInterface&);


public:

    //- Runtime type information
    TypeName("Aitken");


    // Constructors

        //- Construct from components
        AitkenCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~AitkenCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Calculate interface displacement
            virtual void updateDisplacement();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
