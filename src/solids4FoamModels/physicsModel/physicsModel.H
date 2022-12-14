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
    physicsModel

Description
    Virtual base class for physics models, such as fluid, solid and fluid-solid
    solution procedures.

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    physicsModel.C
    newPhysicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef physicsModel_H
#define physicsModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class physicsModel Declaration
\*---------------------------------------------------------------------------*/

class physicsModel
{
    // Private data

        //- Physics model dictionary
        IOdictionary dict_;

        //- Reference to the master time object
        Time& runTime_;

        //- Fluid mesh pointer
        autoPtr<dynamicFvMesh> fluidMeshPtr_;

        //- Solid mesh pointer
        autoPtr<dynamicFvMesh> solidMeshPtr_;

        //- Flag to indicate whether info statements should be printed to
        //  std out or not. Defaults to true but can be modified by the
        //  printInfo() access function
        bool printInfo_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        physicsModel(const physicsModel&);

        //- Disallow default bitwise assignment
        void operator=(const physicsModel&);


public:

    //- Runtime type information
    TypeName("physicsModel");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            physicsModel,
            physicsModel,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );


    // Constructors

        //- Construct from components
        physicsModel
        (
            const word& type,
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
       );


    // Selectors

        //- Select constructor
        static autoPtr<physicsModel> New
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~physicsModel();


    // Member Functions

        // Access

            //- Return const reference to the time
            const Time& runTime() const
            {
                return runTime_;
            }

            //- Return non-const reference to the time
            Time& runTime()
            {
                return runTime_;
            }

            //- Return physicsProperties dictionary
            const dictionary& physicsProperties() const
            {
                return dict_;
            }

            //- Return reference to the printInfo flag
            bool& printInfo()
            {
                return printInfo_;
            }

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime)
            {}

            //- Evolve the mathematical model
            virtual bool evolve() = 0;

            //- This function lets the physics models know the end of the
            //  time-step has been reached
            //  For example, solid mechanical laws may use this
            virtual void updateTotalFields()
            {}

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- This function is called at the end of the simulation
            virtual void end()
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
