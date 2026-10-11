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

#include "velocityFluidMeshMotionCoupler.H"
#include "fixedValuePointPatchFields.H"
#include "compatibilityFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityFluidMeshMotionCoupler, 0);
    addToRunTimeSelectionTable
    (
        fluidMeshMotionCoupler,
        velocityFluidMeshMotionCoupler,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityFluidMeshMotionCoupler::velocityFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const fvMesh& motionMesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    fluidMeshMotionCoupler(mesh, patchIDs, dict),
    pointMotionU_
    (
        const_cast<pointVectorField&>
        (
            motionMesh.objectRegistry::lookupObject<pointVectorField>
            (
                "pointMotionU"
            )
        )
    )
{
    forAll(patchIDs, i)
    {
        if
        (
            !isA<fixedValuePointPatchVectorField>
            (
                pointMotionU_.boundaryField()[patchIDs[i]]
            )
        )
        {
            FatalErrorInFunction
                << "The pointMotionU boundary condition on the fluid "
                << "interface patch "
                << mesh.boundaryMesh()[patchIDs[i]].name()
                << " must be of type fixedValue" << exit(FatalError);
        }
    }
}


Foam::velocityFluidMeshMotionCoupler::velocityFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    velocityFluidMeshMotionCoupler(mesh, mesh, patchIDs, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityFluidMeshMotionCoupler::setPatchMotion
(
    const label patchID,
    const vectorField& iterDisplacement
)
{
    fixedValuePointPatchVectorField& motionUPatch =
        refCast<fixedValuePointPatchVectorField>
        (
            boundaryFieldRef(pointMotionU_)[patchID]
        );

    motionUPatch == iterDisplacement/mesh().time().deltaT().value();
}


// ************************************************************************* //
