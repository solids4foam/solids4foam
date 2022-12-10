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

#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const solidModel& lookupSolidModel(const objectRegistry& obReg)
{
    return lookupSolidModel(obReg, obReg.name());
}


const solidModel& lookupSolidModel
(
    const objectRegistry& obReg,
    const word& obName
)
{
    if (obReg.foundObject<solidModel>("solidModel_" + obName))
    {
        return obReg.lookupObject<solidModel>
        (
            "solidModel_" + obName
        );
    }
    else if (obReg.parent().foundObject<solidModel>("solidModel_" + obName))
    {
        return obReg.parent().lookupObject<solidModel>
        (
            "solidModel_" + obName
        );
    }

    FatalErrorIn
    (
        "const solidModel& lookupSolidModel(const objectRegistry& obReg)"
    )   << "Could not find " << word("solidModel_" + obName) << nl << nl
        << "solidModels in the objectRegistry: "
        << obReg.names<solidModel>() << nl << nl
        << "solidModels in the parent objectRegistry:"
        << obReg.parent().names<solidModel>() << abort(FatalError);

    // Keep the compiler happy
    return obReg.lookupObject<solidModel>("none");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
