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

Class
    Foam::cellZoneInterface

Description
    Functions related to faces straddling cell zones.

Author
    Philip Cardiff, UCD.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef cellZoneInterface_H
#define cellZoneInterface_H

#include "boolList.H"
#include "compatibilityFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //- Return a boolList indicating, for each internal face, if it
    //  straddles two cell zones (true) or not (false)
    tmp<Field<bool>> cellZoneInterface
    (
        const fvMesh& mesh, const bool debug = true
    );

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
