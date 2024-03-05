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

#include "interfaceToInterfaceMapping.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceToInterfaceMapping, 0);
    defineRunTimeSelectionTable(interfaceToInterfaceMapping, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceToInterfaceMapping::interfaceToInterfaceMapping
(
    const word& type,
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalPatchA,
    const globalPolyPatch& globalPatchB
)
:
    dict_(dict),
    patchA_(patchA),
    patchB_(patchB),
    globalPatchA_(globalPatchA),
    globalPatchB_(globalPatchB)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceToInterfaceMapping::~interfaceToInterfaceMapping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceToInterfaceMapping>
Foam::interfaceToInterfaceMapping::New
(
    const word& modelType,
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& zoneA,
    const globalPolyPatch& zoneB
)
{
    Info<< "Selecting interfaceToInterfaceMapping " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "interfaceToInterfaceMapping",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "interfaceToInterfaceMapping::New(Time&, const word&)"
        )   << "Unknown interfaceToInterfaceMapping type " << modelType
            << endl << endl
            << "Valid interfaceToInterfaceMapping types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<interfaceToInterfaceMapping>
    (
        ctorPtr(modelType, dict, patchA, patchB, zoneA, zoneB)
    );
}


void Foam::interfaceToInterfaceMapping::checkFieldSizes
(
    const label& fromZoneSize,
    const label& toZoneSize,
    const label fromFieldSize,
    const label toFieldSize
) const
{
    if (fromFieldSize != fromZoneSize)
    {
        FatalErrorIn
        (
            "void Foam::interfaceToInterfaceMapping::checkFieldSizes\n"
            "(\n"
            "    const label fromFieldSize, const label toFieldSize\n"
            ") const\n"
        )   << "fromField is the wrong size!" << nl
            << "fromField size: " << fromFieldSize
            << ", fromZone size: " << fromZoneSize
            << abort(FatalError);
    }

    if (toFieldSize != toZoneSize)
    {
        FatalErrorIn
        (
            "void Foam::interfaceToInterfaceMapping::checkFieldSizes\n"
            "(\n"
            "    const label fromFieldSize, const label toFieldSize\n"
            ") const\n"
        )   << "toField is wrong size!" << nl
            << "toField size: " << toFieldSize
            << ", toZone size: " << toZoneSize
            << abort(FatalError);
    }
}


// ************************************************************************* //
