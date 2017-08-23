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

\*---------------------------------------------------------------------------*/

#include "physicsModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(physicsModel, 0);
    defineRunTimeSelectionTable(physicsModel, fluid);
    defineRunTimeSelectionTable(physicsModel, solid);
    defineRunTimeSelectionTable(physicsModel, fluidSolidInteraction);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicsModel::physicsModel
(
    const word& type,
    Time& runTime
)
:
    dict_
    (
        IOobject
        (
            "physicsProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    fluidMeshPtr_(),
    solidMeshPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::physicsModel::~physicsModel()
{}

// * * * * * * * * * * * * * * * * Member Members * * * * * * * * * * * * * //

Foam::autoPtr<Foam::physicsModel> Foam::physicsModel::New(Time& runTime)
{
    // Method
    // Read the physicsProperties dictionary and lookup the type
    // If the type is fluid then we create a run-time selectable fluidModel
    // If the type is solid then we create a run-time selectable solidModel
    // If the type is fluidSolidInterface then we create a run-time selectable
    // fluidSolidInterface

    // Read the type
    word physicsTypeName;
    {
        // Read dictionary and ensure it is deleted before the model is
        // created otherwise the dictionary is entered in the database twice
        IOdictionary physicsProperties
        (
            IOobject
            (
                "physicsProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        physicsProperties.lookup("type")
            >> physicsTypeName;
    }

    Info<< "Selecting physicsModel " << physicsTypeName << endl;

    // Currently, the physicsModel can be: fluid, solid or fluidSolidInteraction
    // though more can be added in the future as needed

    if (physicsTypeName == "fluid")
    {
        // Read the fluidModel type
        word fluidTypeName;
        {
            // Read dictionary and ensure it is deleted before the model is
            // created otherwise the dictionary is entered in the database twice
            IOdictionary fluidProperties
            (
                IOobject
                (
                    "fluidProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            fluidProperties.lookup("fluidModel")
                >> fluidTypeName;
        }

        fluidConstructorTable::iterator cstrIter =
            fluidConstructorTablePtr_->find(fluidTypeName);

        if (cstrIter == fluidConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "physicsModel::New(Time&)"
            )   << "Unknown fluidModel type " << fluidTypeName
                << endl << endl
                << "Valid fluidModel types are :" << endl
                << fluidConstructorTablePtr_->toc()
                << exit(FatalError);
        }

        return autoPtr<physicsModel>(cstrIter()(runTime));
    }
    else if (physicsTypeName == "solid")
    {
        // Read the solidModel type
        word solidTypeName;
        {
            // Read dictionary and ensure it is deleted before the model is
            // created otherwise the dictionary is entered in the database twice
            IOdictionary solidProperties
            (
                IOobject
                (
                    "solidProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            solidProperties.lookup("solidModel")
                >> solidTypeName;
        }

        solidConstructorTable::iterator cstrIter =
            solidConstructorTablePtr_->find(solidTypeName);

        if (cstrIter == solidConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "physicsModel::New(Time&)"
            )   << "Unknown solidModel type " << solidTypeName
                << endl << endl
                << "Valid solidModel types are :" << endl
                << solidConstructorTablePtr_->toc()
                << exit(FatalError);
        }

        return autoPtr<physicsModel>(cstrIter()(runTime));
    }
    else if (physicsTypeName == "fluidSolidInteraction")
    {
        // Read the fluidSolidInteractionModel type
        word fluidSolidInteractionTypeName;
        {
            // Read dictionary and ensure it is deleted before the model is
            // created otherwise the dictionary is entered in the database twice
            IOdictionary fluidSolidInteractionProperties
            (
                IOobject
                (
                    "fluidSolidInteractionProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            fluidSolidInteractionProperties.lookup("fluidSolidInterface")
                >> fluidSolidInteractionTypeName;
        }

        fluidSolidInteractionConstructorTable::iterator cstrIter =
            fluidSolidInteractionConstructorTablePtr_->find
            (
                fluidSolidInteractionTypeName
            );

        if (cstrIter == fluidSolidInteractionConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "physicsModel::New(Time&)"
            )   << "Unknown fluidSolidInteractionModel type "
                << fluidSolidInteractionTypeName
                << endl << endl
                << "Valid fluidSolidInteractionModel types are :" << endl
                << fluidSolidInteractionConstructorTablePtr_->toc()
                << exit(FatalError);
        }

        return autoPtr<physicsModel>(cstrIter()(runTime));
    }
    else
    {
        FatalErrorIn
        (
            "Foam::autoPtr<Foam::physicsModel> Foam::physicsModel::New\n"
            "(\n"
            "    Time& runTime\n"
            ")"
        )   << "Unknown physicsModel type " << physicsTypeName << endl << endl
            << "Valid physicalModel types are : "
            << "fluid, solid, fluidSolidInteraction" << exit(FatalError);
    }

    // Keep the compiler happy
    fluidConstructorTable::iterator cstrIter;
    return autoPtr<physicsModel>(cstrIter()(runTime));
}


// ************************************************************************* //
