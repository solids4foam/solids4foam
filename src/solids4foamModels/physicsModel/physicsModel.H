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

#include "foamTime.H"
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

    // Private Member Functions

        //- Disallow default bitwise copy construct
        physicsModel(const physicsModel&);

        //- Disallow default bitwise assignment
        void operator=(const physicsModel&);


public:

    //- Runtime type information
    TypeName("physicsModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            physicsModel,
            fluid,
            (
                Time& runTime
            ),
            (runTime)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            physicsModel,
            solid,
            (
                Time& runTime
            ),
            (runTime)
        );

            declareRunTimeSelectionTable
        (
            autoPtr,
            physicsModel,
            fluidSolidInteraction,
            (
                Time& runTime
            ),
            (runTime)
        );

    // Constructors

        //- Construct from components
        physicsModel
        (
            const word& type,
            Time& runTime
        );


    // Selectors

        //- Select constructor
        static autoPtr<physicsModel> New
        (
            Time& runTime
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

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT()
            {}

            //- Evolve the mathematical model
            virtual bool evolve() = 0;

            //- Write fields
            virtual void writeFields(const Time& runTime)
            {
                runTime.write();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
