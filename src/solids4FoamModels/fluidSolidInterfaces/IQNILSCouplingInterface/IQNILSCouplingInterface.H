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
    IQNILSCouplingInterface

Description
    Strong IQN-ILS Dirichlet-Neumann coupling, as developed by Degroote et al.

    A fluid solid interaction solver with IQN-ILS coupling algorithm.
    J. Degroote, K.-J. Bathe and J. Vierendeels.
    Performance of a new partitioned procedure versus a monolithic
    procedure in fluid-solid interaction. Computers & Solids

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    IQNILSCouplingInterface.C

\*---------------------------------------------------------------------------*/

#ifndef IQNILSCouplingInterface_H
#define IQNILSCouplingInterface_H

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

/*---------------------------------------------------------------------------*\
                        Class IQNILSCouplingInterface Declaration
\*---------------------------------------------------------------------------*/

class IQNILSCouplingInterface
:
    public fluidSolidInterface
{
    // Private data

        //- Fixed under-relaxation factor
        const scalar relaxationFactor_;

        //- Coupling reuse
        label couplingReuse_;

        //- Predict solid
        const bool predictSolid_;

        //- List of coupling field V
        List<DynamicList<vectorField> > fluidPatchesPointsV_;

        //- List of coupling field W
        List<DynamicList<vectorField> > fluidPatchesPointsW_;

        //- List of coupling field T
        List<DynamicList<scalar> > fluidPatchesPointsT_;


    // Private Member Functions

        //- Reuse coupling
        label couplingReuse() const;

        //- Disallow default bitwise copy construct
        IQNILSCouplingInterface(const IQNILSCouplingInterface&);

        //- Disallow default bitwise assignment
        void operator=(const IQNILSCouplingInterface&);


public:

    //- Runtime type information
    TypeName("IQNILS");


    // Constructors

        //- Construct from components
        IQNILSCouplingInterface
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~IQNILSCouplingInterface()
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
