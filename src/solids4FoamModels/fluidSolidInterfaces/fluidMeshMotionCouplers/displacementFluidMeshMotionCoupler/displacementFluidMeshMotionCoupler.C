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

#ifdef OPENFOAM_NOT_EXTEND

#include "displacementFluidMeshMotionCoupler.H"
#include "fixedValuePointPatchFields.H"
#include "compatibilityFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementFluidMeshMotionCoupler, 0);
    addToRunTimeSelectionTable
    (
        fluidMeshMotionCoupler,
        displacementFluidMeshMotionCoupler,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementFluidMeshMotionCoupler::displacementFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    fluidMeshMotionCoupler(mesh, patchIDs, dict),
    motionSolver_(refCast<displacementMotionSolver>(*lookupMotionSolver(mesh)))
{
    forAll(patchIDs, i)
    {
        if
        (
            !isA<fixedValuePointPatchVectorField>
            (
                motionSolver_.pointDisplacement().boundaryField()[patchIDs[i]]
            )
        )
        {
            FatalErrorInFunction
                << "The pointDisplacement boundary condition on the fluid "
                << "interface patch "
                << mesh.boundaryMesh()[patchIDs[i]].name()
                << " must be of type fixedValue" << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::displacementFluidMeshMotionCoupler::totalPatchDisplacement
(
    const polyPatch& patch,
    const pointField& points0,
    const vectorField& iterDisplacement
)
{
    return
        patch.localPoints()
      + iterDisplacement
      - vectorField(points0, patch.meshPoints());
}


void Foam::displacementFluidMeshMotionCoupler::setPatchMotion
(
    const label patchID,
    const vectorField& iterDisplacement
)
{
    fixedValuePointPatchVectorField& pointDisplacementPatch =
        refCast<fixedValuePointPatchVectorField>
        (
            boundaryFieldRef(motionSolver_.pointDisplacement())[patchID]
        );

    pointDisplacementPatch ==
        totalPatchDisplacement
        (
            mesh().boundaryMesh()[patchID],
            motionSolver_.points0(),
            iterDisplacement
        );
}

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
