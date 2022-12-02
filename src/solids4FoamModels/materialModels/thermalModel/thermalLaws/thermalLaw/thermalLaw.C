/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "thermalLaw.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalLaw, 0);
defineRunTimeSelectionTable(thermalLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalLaw::thermalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh)
{}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

autoPtr<thermalLaw> thermalLaw::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word thermalTypeName(dict.lookup("type"));

    Info<< "Selecting thermal model " << thermalTypeName << endl;

#if OPENFOAM > 2205
    auto* cstrIter = dictionaryConstructorTable(thermalTypeName);
    
    if (!cstrIter){
        FatalIOErrorIn(
            "thermalLaw::New(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown thermalLaw type "
            << thermalTypeName << endl << endl
            << "Valid  thermalLaws are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }
    return autoPtr<thermalLaw>(cstrIter(name, mesh, dict));
#else    
    
     dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermalTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "thermalLaw::New(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown thermalLaw type "
            << thermalTypeName << endl << endl
            << "Valid  thermalLaws are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<thermalLaw>(cstrIter()(name, mesh, dict));
#endif
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
