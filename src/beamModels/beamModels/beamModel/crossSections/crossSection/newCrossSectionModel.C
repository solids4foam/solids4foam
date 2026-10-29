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

InClass
    crossSectionModel

\*---------------------------------------------------------------------------*/

#include "crossSectionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<crossSectionModel> crossSectionModel::New
(
    const word& name,
    const dictionary& dict
)
{
    Info<< "    Cross-section model: " << name << endl;

    auto* cstrIter = dictionaryConstructorTable(name);

    if (!cstrIter)
    {
        FatalIOErrorIn
        (
            "crossSectionModel::New(\n"
            "    const word& name,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown crossSectionModel type "
            << name << endl << endl
            << "Valid  crossSectionModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<crossSectionModel>
    (
        cstrIter
        (
            name,
            dict
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
