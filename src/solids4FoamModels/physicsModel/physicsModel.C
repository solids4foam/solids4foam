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

\*---------------------------------------------------------------------------*/

#include "physicsModel.H"
#ifdef OPENFOAM_ORG
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
    // NB: dictionary must be unregistered to avoid adding to the database

    IOdictionary props
    (
        IOobject
        (
            "physicsProperties",
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // Do not register
        )
    );

    word modelType(props.lookup("type"));

    // For backwards compatibility, update names
    // This means the user can equivalently select "solid" or "solidModel", etc.
    if (modelType == "solid")
    {
        modelType = "solidModel";
    }
    else if (modelType == "fluid")
    {
        modelType = "fluidModel";
    }
    else if (modelType == "fluidSolidInteraction")
    {
        modelType = "fluidSolidInterface";
    }

    Info<< "Selecting physicsModel " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = physicsModelConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            props,
            "physicsModel",
            modelType,
            *physicsModelConstructorTablePtr_
        ) << exit(FatalIOError);
    }

#else
    physicsModelConstructorTable::iterator cstrIter =
        physicsModelConstructorTablePtr_->find(modelType);

    if (cstrIter == physicsModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "physicsModel::New(Time&)"
        )   << "Unknown physicsModel type " << modelType
            << endl << endl
            << "Valid physicsModel types are :" << endl
            << physicsModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<physicsModel>(ctorPtr(runTime, region));
}


void Foam::physicsModel::writeFields(const Time& runTime)
{
    runTime.write();
}


// ************************************************************************* //
