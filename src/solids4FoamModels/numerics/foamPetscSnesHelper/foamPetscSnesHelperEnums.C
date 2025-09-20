/*---------------------------------------------------------------------------* \
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

#include "foamPetscSnesHelper.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#ifdef OPENFOAM_COM
    const Enum<foamPetscSnesHelper::solutionLocation>
    foamPetscSnesHelper::solutionLocationNames_
    ({
        {
            foamPetscSnesHelper::solutionLocation::CELLS,
            "cells"
        },
        {
            foamPetscSnesHelper::solutionLocation::POINTS,
            "points"
        },
        {
            foamPetscSnesHelper::solutionLocation::NONE,
            "none"
        },
    });
#else
    template<>
    const char* NamedEnum<foamPetscSnesHelper::solutionLocation, 3>::names[] =
    {
     	"cells",
        "points",
	"none"
    };
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
