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

Description
    Lookup the current beamModel from the given objectRegistry and return it as
    a reference.

    If the beamModel is not found in the current registry then it will check
    parent registry.

    It assumes that the beamModel is called beamModel_<REGION_NAME> where
    <REGION_NAME> is taken from obReg.name(). e.g. beamModel_region0 or
    beamModel_solid.

    This facilitates multiple beamModels to co-exist.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef lookupBeamModel_H
#define lookupBeamModel_H

#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const beamModel& lookupBeamModel(const objectRegistry& obReg);

const beamModel& lookupBeamModel
(
    const objectRegistry& obReg,
    const word& obName
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
