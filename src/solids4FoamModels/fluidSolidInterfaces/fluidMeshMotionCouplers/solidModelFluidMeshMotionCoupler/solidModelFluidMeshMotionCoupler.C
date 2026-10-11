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

#ifdef OPENFOAM_COM

#include "solidModelFluidMeshMotionCoupler.H"
#include "displacementFluidMeshMotionCoupler.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"
#include "compatibilityFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidModelFluidMeshMotionCoupler, 0);
    addToRunTimeSelectionTable
    (
        fluidMeshMotionCoupler,
        solidModelFluidMeshMotionCoupler,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModelFluidMeshMotionCoupler::solidModelFluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    fluidMeshMotionCoupler(mesh, patchIDs, dict),
    motionSolver_
    (
        refCast<meshMotionSolidModelFvMotionSolver>(*lookupMotionSolver(mesh))
    ),
    referencePatchIDs_(mesh.boundaryMesh().size(), -1)
{
    if (!motionSolver_.referenceMesh())
    {
        FatalErrorInFunction
            << "The solids4foamSolidModel mesh motion solver must solve the "
            << "solid model on a separate reference copy of the fluid mesh: "
            << "set regionName to a region other than "
            << mesh.name() << exit(FatalError);
    }

    const solidModel& model = motionSolver_.model();
    const polyBoundaryMesh& referenceBoundary = model.mesh().boundaryMesh();

    forAll(patchIDs, i)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchIDs[i]];
        const label referencePatchID = referenceBoundary.findPatchID
        (
            patch.name()
        );

        if
        (
            referencePatchID < 0
         || referenceBoundary[referencePatchID].size() != patch.size()
        )
        {
            FatalErrorInFunction
                << "The fluid interface patch " << patch.name()
                << " is not a patch of the same size in the mesh motion "
                << "reference mesh " << model.mesh().name()
                << exit(FatalError);
        }

        referencePatchIDs_[patchIDs[i]] = referencePatchID;

        const fvPatchVectorField& DPatch =
            model.D().boundaryField()[referencePatchID];

        if
        (
            !isA<fixedDisplacementFvPatchVectorField>(DPatch)
         || refCast<const fixedDisplacementFvPatchVectorField>
            (
                DPatch
            ).enforcePointPatchField()
         || refCast<const fixedDisplacementFvPatchVectorField>
            (
                DPatch
            ).timeVarying()
        )
        {
            FatalErrorInFunction
                << "The " << model.D().name() << " boundary condition on "
                << "the patch " << patch.name() << " of the mesh motion "
                << "reference mesh must be of type fixedDisplacement, with "
                << "enforcePointPatchField no and no displacementSeries"
                << exit(FatalError);
        }

        if
        (
            !isA<fixedValuePointPatchVectorField>
            (
                model.pointD().boundaryField()[referencePatchID]
            )
        )
        {
            FatalErrorInFunction
                << "The " << model.pointD().name() << " boundary condition "
                << "on the patch " << patch.name() << " of the mesh motion "
                << "reference mesh must be of type fixedValue"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidModelFluidMeshMotionCoupler::setPatchMotion
(
    const label patchID,
    const vectorField& iterDisplacement
)
{
    solidModel& model = motionSolver_.model();
    const label referencePatchID = referencePatchIDs_[patchID];

    const vectorField totalDisp
    (
        displacementFluidMeshMotionCoupler::totalPatchDisplacement
        (
            mesh().boundaryMesh()[patchID],
            motionSolver_.points0(),
            iterDisplacement
        )
    );

    // The reference mesh is a copy of the fluid mesh, so the patch points
    // are in the same order
    fixedValuePointPatchVectorField& pointDPatch =
        refCast<fixedValuePointPatchVectorField>
        (
            boundaryFieldRef(model.pointD())[referencePatchID]
        );

    pointDPatch == totalDisp;

    fixedDisplacementFvPatchVectorField& DPatch =
        refCast<fixedDisplacementFvPatchVectorField>
        (
            boundaryFieldRef(model.D())[referencePatchID]
        );

    primitivePatchInterpolation interp
    (
        model.mesh().boundaryMesh()[referencePatchID]
    );

    DPatch.totalDisp() = interp.pointToFaceInterpolate(totalDisp);
}

#endif // OPENFOAM_COM

// ************************************************************************* //
