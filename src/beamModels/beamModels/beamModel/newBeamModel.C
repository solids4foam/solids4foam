/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamModel.H"
//#include "foamTime.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Stop with a clear message if the beam fields W or Theta still use the
// boundary conditions renamed in #61. The old names cannot be kept as aliases,
// as they clash with the solids4foam classes of the same name (issue #56).
// The fields are read as plain dictionaries before the beam model reads them:
// an old name is otherwise either an unknown type or, with solids4foam loaded,
// silently built by the solids4foam class.
// Temporary migration help, to be removed in a later release.
void checkRenamedBeamBoundaryConditions
(
    const Foam::Time& runTime,
    const Foam::word& region
)
{
    using namespace Foam;

    const HashTable<word> renamed
    ({
        {"fixedDisplacement", "fixedBeamDisplacement"},
        {"fixedRotation", "fixedBeamRotation"}
    });

    for (const word fieldName : {"W", "Theta"})
    {
        IOobject fieldHeader
        (
            fieldName,
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.timeName())
          : fileName(runTime.timeName()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // A missing field is reported when the beam model reads it
        if (!fieldHeader.typeHeaderOk<volVectorField>(false))
        {
            continue;
        }

        autoPtr<ISstream> isPtr
        (
            fileHandler().NewIFstream(fieldHeader.objectPath())
        );

        const dictionary fieldDict(*isPtr);

        const dictionary* bfDictPtr = fieldDict.findDict("boundaryField");

        if (!bfDictPtr)
        {
            continue;
        }

        for (const entry& patchEntry : *bfDictPtr)
        {
            if (!patchEntry.isDict())
            {
                continue;
            }

            const word patchType
            (
                patchEntry.dict().getOrDefault<word>("type", word::null)
            );

            const auto iter = renamed.cfind(patchType);

            if (iter.good())
            {
                FatalIOErrorInFunction(patchEntry.dict())
                    << "Patch " << patchEntry.keyword() << " of beam field "
                    << fieldName << " uses the boundary condition type "
                    << patchType << ", which has been renamed to "
                    << iter.val() << "." << nl
                    << "Replace '" << patchType << "' with '" << iter.val()
                    << "' in " << fieldHeader.objectPath() << nl
                    << "and in the W and Theta files of any other time or"
                    << " processor directories." << nl
                    << "See https://github.com/solids4foam/beamFoam"
                    << "#breaking-change-beam-boundary-condition-names"
                    << exit(FatalIOError);
            }
        }
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::beamModel> Foam::beamModel::New
(
    Time& runTime,
    const word& region
)
{
    checkRenamedBeamBoundaryConditions(runTime, region);

    word beamModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the flow is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary beamProperties
        (
            IOobject
            (
                "beamProperties",
                bool(region == dynamicFvMesh::defaultRegion)
              ? fileName(runTime.caseConstant())
              : fileName(runTime.caseConstant()/region),
                // runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        beamProperties.lookup("beamModel")
            >> beamModelTypeName;
    }

    Info<< "Selecting beamModel " << beamModelTypeName << endl;

    // changes to be compiled with OFv23
    auto* cstrIter = dictionaryConstructorTable(beamModelTypeName);

    if (!cstrIter)
    {
        FatalErrorIn
	(
            "beamModel::New(Time&, const word&)"
        )   << "Unknown beamModel type " << beamModelTypeName
            << endl << endl
            << "Valid beamModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }



    //dictionaryConstructorTable::iterator cstrIter =
    //    dictionaryConstructorTablePtr_->find(beamModelTypeName);

    //if (cstrIter == dictionaryConstructorTablePtr_->end())
    //{
    //    FatalErrorIn
    //    (
    //        "beamModel::New(Time&, const word&)"
    //    )   << "Unknown beamModel type " << beamModelTypeName
    //        << endl << endl
    //        << "Valid beamModel types are :" << endl
    //        << dictionaryConstructorTablePtr_->toc()
    //        << exit(FatalError);
    //}

    return autoPtr<beamModel>(cstrIter(runTime, region));
}


// ************************************************************************* //
