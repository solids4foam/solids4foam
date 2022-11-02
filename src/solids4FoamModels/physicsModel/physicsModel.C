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
#ifdef OPENFOAMFOUNDATION
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(physicsModel, 0);
    defineRunTimeSelectionTable(physicsModel, physicsModel);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicsModel::physicsModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    dict_
    (
        IOobject
        (
            "physicsProperties",
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    fluidMeshPtr_(),
    solidMeshPtr_(),
    printInfo_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::physicsModel::~physicsModel()
{}

// * * * * * * * * * * * * * * * * Member Members * * * * * * * * * * * * * //

Foam::autoPtr<Foam::physicsModel> Foam::physicsModel::New
(
    Time& runTime,
    const word& region
)
{
    // Read the model type
    word physicsModelTypeName;
    {
        // Read dictionary and ensure it is deleted before the model is
        // created otherwise the dictionary is entered in the database twice
        IOdictionary physicsProperties
        (
            IOobject
            (
                "physicsProperties",
                bool(region == dynamicFvMesh::defaultRegion)
              ? fileName(runTime.caseConstant())
              : fileName(runTime.caseConstant()/region),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        physicsProperties.lookup("type")
            >> physicsModelTypeName;
    }

    // For backwards compatibility, update names
    // This means the user can equivalently select "solid" or "solidModel", etc.
    if (physicsModelTypeName == "solid")   {  physicsModelTypeName = "solidModel";   }
    else if (physicsModelTypeName == "fluid") { physicsModelTypeName = "fluidModel";  }
    else if (physicsModelTypeName == "fluidSolidInteraction")  { physicsModelTypeName = "fluidSolidInterface";  }

    Info<< "Selecting physicsModel " << physicsModelTypeName << endl;
#if OPENFOAM > 2005    
    auto* cstrIter = physicsModelConstructorTable(physicsModelTypeName);
    if (!cstrIter) {
            FatalErrorIn ( "physicsModel::New(Time&)")   << "Unknown physicsModel type " << physicsModelTypeName
                << endl << endl  << "Valid physicsModel types are :" << endl   << physicsModelConstructorTablePtr_->toc()
                << exit(FatalError);
        }
     return autoPtr<physicsModel>(cstrIter(runTime, region));   
#else    
    physicsModelConstructorTable::iterator cstrIter =
        physicsModelConstructorTablePtr_->find(physicsModelTypeName);

    if (cstrIter == physicsModelConstructorTablePtr_->end()) {
        FatalErrorIn( "physicsModel::New(Time&)" )   << "Unknown physicsModel type " << physicsModelTypeName
            << endl << endl   << "Valid physicsModel types are :" << endl   << physicsModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    return autoPtr<physicsModel>(cstrIter()(runTime, region));
#endif    
}


void Foam::physicsModel::writeFields(const Time& runTime)
{
    runTime.write();
}


// ************************************************************************* //
