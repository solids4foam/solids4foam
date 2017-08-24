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
    weakCouplingInterface

Description
    Weak Dirichlet-Neumann coupling, where no outer FSI correctors are
    performed each time-step.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    weakCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef weakCouplingInterface_H
#define weakCouplingInterface_H

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class weakCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class weakCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- Solid zone traction (pressure + vicous)
        vectorField solidZoneTraction_;

        //- Solid zone traction (pressure + vicous)
        vectorField solidZoneTractionPrev_;

        //- Solid zone traction (pressure + vicous)
        vectorField predictedSolidZoneTraction_;

        //- Fixed under-relaxation factor
        const scalar relaxationFactor_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        weakCouplingInterface(const weakCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const weakCouplingInterface&);


public:

    //- Runtime type information
    TypeName("weakCoupling");


    // Constructors

        //- Construct from components
        weakCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~weakCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Initialize fields
            virtual void initializeFields();

            //- Calculate interface displacement
            virtual void updateWeakDisplacement();

            //- Update interface traction
            virtual void updateWeakTraction();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
