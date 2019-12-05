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
    const word& type,
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& zoneA,
    const globalPolyPatch& zoneB
)
{
    Info<< "Selecting interfaceToInterfaceMapping " << type << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "interfaceToInterfaceMapping::New(Time&, const word&)"
        )   << "Unknown interfaceToInterfaceMapping type " << type
            << endl << endl
            << "Valid interfaceToInterfaceMapping types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<interfaceToInterfaceMapping>
    (
        cstrIter()(type, dict, patchA, patchB, zoneA, zoneB)
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
