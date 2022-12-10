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
    fixedRelaxationCouplingInterface

Description
    Strong Dirichlet-Neumann coupling with fixed under-relaxation.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    fixedRelaxationCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRelaxationCouplingInterface_H
#define fixedRelaxationCouplingInterface_H

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class fixedRelaxationCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class fixedRelaxationCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- Fixed under-relaxation factor
        const scalar relaxationFactor_;

        //- Predict solid
        const bool predictSolid_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fixedRelaxationCouplingInterface
        (
            const fixedRelaxationCouplingInterface&
        );

        //- Disallow default bitwise assignment
        void operator=(const fixedRelaxationCouplingInterface&);


public:

    //- Runtime type information
    TypeName("fixedRelaxation");


    // Constructors

        //- Construct from components
        fixedRelaxationCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~fixedRelaxationCouplingInterface()
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
