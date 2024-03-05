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
    Foam::simpleCrackerFvMesh

Description
    Dynamic mesh class that calls the updateCrack function of the simple
    cohesiveZoneIncremental mode-I boundary conditions.

    Note: uses const_cast to call updateCrack.

SourceFiles
    simpleCrackerFvMesh.C

Author
    Philip Cardiff, UCD. All rights reserved, UT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef simpleCrackerFvMesh_H
#define simpleCrackerFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleCrackerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class simpleCrackerFvMesh
:
    public dynamicFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        simpleCrackerFvMesh(const simpleCrackerFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const simpleCrackerFvMesh&);


public:

    //- Runtime type information
    TypeName("simpleCrackerFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit simpleCrackerFvMesh(const IOobject& io);


    // Destructor

        virtual ~simpleCrackerFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
