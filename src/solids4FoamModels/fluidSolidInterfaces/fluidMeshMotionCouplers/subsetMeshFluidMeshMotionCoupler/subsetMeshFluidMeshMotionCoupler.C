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

#ifndef OPENFOAM_NOT_EXTEND

#include "subsetMeshFluidMeshMotionCoupler.H"
#include "newSubsetMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(subsetMeshFluidMeshMotionCoupler, 0);
    addToRunTimeSelectionTable
    (
        fluidMeshMotionCoupler,
        subsetMeshFluidMeshMotionCoupler,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsetMeshFluidMeshMotionCoupler::subsetMeshFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    velocityFluidMeshMotionCoupler
    (
        mesh,
        refCast<newSubsetMotionSolverFvMesh>(mesh).subsetMesh().subMesh(),
        patchIDs,
        dict
    )
{}

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
