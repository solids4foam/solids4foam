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

Function
    Foam::bdfFsiCheck

Description
    Stops with a fatal error if the fluid ddt scheme is BDF.

    The Robin interface updates of pimpleFluid, fluidModel and interFluid,
    elasticWallVelocity, elasticSlipWallVelocity, newMovingWallVelocity and
    the Robin coefficient of elasticWallPressure choose their form by the
    name of the ddt scheme, Euler or backward. With BDF they would take their
    Euler form or skip the update without saying so, and the FSI coupling
    would not be consistent with the fluid ddt. These sites call this check
    until they take their weights from bdfDdtScheme::fieldWeights instead.

    Fluid-only cases, on static or moving meshes, do not reach these sites
    unless they use newMovingWallVelocity; OpenFOAM's movingWallVelocity,
    whose normal wall velocity comes from meshPhi, can be used with BDF.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef bdfFsiCheck_H
#define bdfFsiCheck_H

#include "error.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Stop if ddtScheme, the scheme name read for ddtName, is BDF. The caller
//  is named in the message by functionName.
inline void bdfFsiCheck
(
    const word& ddtScheme,
    const word& ddtName,
    const char* functionName
)
{
    if (ddtScheme == "BDF")
    {
        FatalErrorIn(functionName)
            << "The ddt scheme of " << ddtName << " is BDF, which is not yet"
            << " consistent with the fluid-solid interaction coupling: this"
            << " part of it supports only the Euler and backward schemes."
            << nl << "Use Euler or backward for " << ddtName
            << " in fluid-solid interaction cases."
            << exit(FatalError);
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
