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

#include "uniformTranslation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{
    defineTypeNameAndDebug(uniformTranslation, 0);
    addToRunTimeSelectionTable
    (
        immersedBodyMotion,
        uniformTranslation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBodyMotions::uniformTranslation::uniformTranslation
(
    const dictionary& dict
)
:
    immersedBodyMotion(dict),
    velocity_(dict.get<vector>("velocity")),
    startTime_(dict.getOrDefault<scalar>("startTime", 0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::immersedBodyMotions::uniformTranslation::points
(
    const pointField& points0,
    const scalar t
) const
{
    return points0 + velocity_*(t - startTime_);
}


Foam::tmp<Foam::vectorField>
Foam::immersedBodyMotions::uniformTranslation::velocity
(
    const pointField& x,
    const scalar t
) const
{
    return tmp<vectorField>::New(x.size(), velocity_);
}


// ************************************************************************* //
