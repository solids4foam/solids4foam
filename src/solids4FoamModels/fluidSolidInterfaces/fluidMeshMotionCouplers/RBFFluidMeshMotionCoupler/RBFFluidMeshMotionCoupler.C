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

#if !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)

#include "RBFFluidMeshMotionCoupler.H"
#include "primitivePatchInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBFFluidMeshMotionCoupler, 0);
    addToRunTimeSelectionTable
    (
        fluidMeshMotionCoupler,
        RBFFluidMeshMotionCoupler,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFFluidMeshMotionCoupler::RBFFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    fluidMeshMotionCoupler(mesh, patchIDs, dict),
    motionSolver_
    (
        refCast<RBFMeshMotionSolver>
        (
            *lookupMotionSolver(mesh)
        )
    ),
    motion_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBFFluidMeshMotionCoupler::setPatchMotion
(
    const label patchID,
    const vectorField& iterDisplacement
)
{
    // The patches without a prescribed motion do not move
    if (motion_.empty())
    {
        motion_.setSize(mesh().boundaryMesh().size());

        forAll(mesh().boundaryMesh(), patchI)
        {
            motion_[patchI] = vectorField
            (
                mesh().boundaryMesh()[patchI].size(), vector::zero
            );
        }
    }

    primitivePatchInterpolation interp(mesh().boundaryMesh()[patchID]);

    motion_[patchID] = interp.pointToFaceInterpolate(iterDisplacement);
}


void Foam::RBFFluidMeshMotionCoupler::commit()
{
    // Note: the RBF solver takes the displacement as opposed to the velocity
    motionSolver_.setMotion(motion_);

    motion_.clear();
}

#endif // !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)

// ************************************************************************* //
