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
    solidPolyMeshLduAddressing

Description
    Addressing for solidPolyMesh allowing implicit treatment of cell point cell
    neighbours.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef solidPolyMeshLduAddressing_H
#define solidPolyMeshLduAddressing_H

#include "lduAddressing.H"
#include "solidPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class solidPolyMeshLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class solidPolyMeshLduAddressing
:
    public lduAddressing
{
    // Private data

        //- Lower as a subList of allOwner
        labelList lowerAddr_;

        //- Upper as a reference to neighbour
        labelList upperAddr_;

        //- Patch addressing as a list of sublists
        labelListList patchAddr_;

        //- Patch field evaluation schedule
        const lduSchedule& patchSchedule_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidPolyMeshLduAddressing
        (
            const solidPolyMeshLduAddressing&
        );

        //- Disallow default bitwise assignment
        void operator=(const solidPolyMeshLduAddressing&);


public:

    // Constructors

        //- Construct from components
        solidPolyMeshLduAddressing(const solidPolyMesh& mesh);


    // Destructor

        virtual ~solidPolyMeshLduAddressing()
        {}


    // Member Functions

        //- Return number of interfaces
        virtual label nPatches() const
        {
            return patchAddr_.size();
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        virtual const unallocLabelList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        virtual const unallocLabelList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return patch addressing
        virtual const unallocLabelList& patchAddr(const label i) const
        {
            return patchAddr_[i];
        }

        // Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
