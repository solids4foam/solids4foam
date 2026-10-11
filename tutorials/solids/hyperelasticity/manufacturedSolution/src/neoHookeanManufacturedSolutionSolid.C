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

\*----------------------------------------------------------------------------*/

#include "neoHookeanManufacturedSolutionSolid.H"

#ifndef OPENFOAM_COM

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{
    defineTypeNameAndDebug(neoHookeanManufacturedSolutionSolid, 0);
    addToRunTimeSelectionTable
    (
        solidModel,
        neoHookeanManufacturedSolutionSolid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModels::neoHookeanManufacturedSolutionSolid::
neoHookeanManufacturedSolutionSolid
(
    Time& runTime,
    const word& region
)
:
    nonLinGeomTotalLagTotalDispSolid(runTime, region),
    mms_(mesh())
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::solidModels::neoHookeanManufacturedSolutionSolid::fvOptionsSource() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            "neoHookeanManufacturedSolutionSource",
            mms_.bodyForces()
        )
    );
}

#endif

// ************************************************************************* //
