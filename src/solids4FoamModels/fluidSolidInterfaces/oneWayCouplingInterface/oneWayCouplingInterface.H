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
    oneWayCouplingInterface

Description
    Pseduo coupling scheme class that reads the fluid velocity and pressure
    fields from a pre-run fluid case; this then allows a one-way fluid-solid
    interaction simulation.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    oneWayCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oneWayCouplingInterface_H
#define oneWayCouplingInterface_H

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class oneWayCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class oneWayCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- List of solid zones traction (pressure + vicous)
        List<vectorField> solidZonesTraction_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oneWayCouplingInterface(const oneWayCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const oneWayCouplingInterface&);


public:

    //- Runtime type information
    TypeName("oneWayCoupling");


    // Constructors

        //- Construct from components
        oneWayCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~oneWayCouplingInterface()
        {}


    // Member Functions

        // Edit

            //- Evolve the interface
            virtual bool evolve();

            //- Update interface traction
            virtual void updateTraction();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
