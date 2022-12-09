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

        //- Ptr list of solid zones traction (pressure + vicous)
        //PtrList<vectorIOField> solidZonesTractionPtrList_;

        //- List of solid zones traction (pressure + vicous)
        //List<vectorField> solidZonesTractionPrev_;

        //- Ptr list of solid zones traction (pressure + vicous)
        //PtrList<vectorIOField> predictedSolidZonesTractionPtrList_;


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
