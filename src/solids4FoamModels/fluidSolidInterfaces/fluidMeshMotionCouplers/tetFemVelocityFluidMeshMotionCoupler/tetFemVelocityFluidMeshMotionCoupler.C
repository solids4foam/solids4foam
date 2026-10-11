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

#include "tetFemVelocityFluidMeshMotionCoupler.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetPolyPatchInterpolation.H"
#include "tetFemMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetFemVelocityFluidMeshMotionCoupler, 0);
    addToRunTimeSelectionTable
    (
        fluidMeshMotionCoupler,
        tetFemVelocityFluidMeshMotionCoupler,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetFemVelocityFluidMeshMotionCoupler::
tetFemVelocityFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    fluidMeshMotionCoupler(mesh, patchIDs, dict),
    motionU_
    (
        const_cast<tetPointVectorField&>
        (
            mesh.objectRegistry::lookupObject<tetPointVectorField>("motionU")
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetFemVelocityFluidMeshMotionCoupler::setPatchMotion
(
    const label patchID,
    const vectorField& iterDisplacement
)
{
    fixedValueTetPolyPatchVectorField& motionUPatch =
        refCast<fixedValueTetPolyPatchVectorField>
        (
            motionU_.boundaryField()[patchID]
        );

    tetPolyPatchInterpolation tppi
    (
        refCast<const faceTetPolyPatch>(motionUPatch.patch())
    );

    motionUPatch ==
        tppi.pointToPointInterpolate(iterDisplacement)
       /mesh().time().deltaT().value();
}

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
