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

Description
    Lookup the current solidModel from the given objectRegistry and return it as
    a reference.

    If the solidModel is not found in the current registry then it will check
    parent registry.

    It assumes that the solidModel is called solidModel_<REGION_NAME> where
    <REGION_NAME> is taken from obReg.name(). e.g. solidModel_region0 or
    solidModel_solid.

    This facilitates multiple solidModels to co-exist.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef lookupSolidModel_H
#define lookupSolidModel_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const solidModel& lookupSolidModel(const objectRegistry& obReg);

const solidModel& lookupSolidModel
(
    const objectRegistry& obReg,
    const word& obName
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
