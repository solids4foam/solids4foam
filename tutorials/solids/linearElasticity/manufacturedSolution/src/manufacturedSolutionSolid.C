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

#include "manufacturedSolutionSolid.H"

#ifndef OPENFOAM_COM

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{
    defineTypeNameAndDebug(manufacturedSolutionSolid, 0);
    addToRunTimeSelectionTable
    (
        solidModel,
        manufacturedSolutionSolid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModels::manufacturedSolutionSolid::manufacturedSolutionSolid
(
    Time& runTime,
    const word& region
)
:
    linGeomTotalDispSolid(runTime, region),
    sourceDict_
    (
        IOobject
        (
            "fvOptions",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mms_(mesh(), sourceDict_.subDict("momentumSource"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::solidModels::manufacturedSolutionSolid::fvOptionsSource() const
{
    return -mms_.bodyForces();
}

#endif

// ************************************************************************* //
