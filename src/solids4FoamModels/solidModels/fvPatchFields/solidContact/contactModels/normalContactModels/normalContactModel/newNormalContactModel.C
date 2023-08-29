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

InClass
    normalContactModel

\*---------------------------------------------------------------------------*/

#include "normalContactModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<normalContactModel> normalContactModel::New
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
{
    Info<< "    Normal contact model: " << name << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "normalContactModel::New(\n"
            "    const word& name,\n"
            "    const fvPatch& patch,\n"
            "    const dictionary& dict,\n"
            "    const label masterPatchID,\n"
            "    const label slavePatchID,\n"
            "    const standAlonePatch& masterFaceZonePatch,\n"
            "    const standAlonePatch& slaveFaceZonePatch\n"
            ")",
            dict
        )   << "Unknown normalContactModel type "
            << name << endl << endl
            << "Valid  normalContactModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<normalContactModel>
    (
        cstrIter()
        (
            name,
            patch,
            dict,
            masterPatchID,
            slavePatchID,
            masterFaceZonePatch,
            slaveFaceZonePatch
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
