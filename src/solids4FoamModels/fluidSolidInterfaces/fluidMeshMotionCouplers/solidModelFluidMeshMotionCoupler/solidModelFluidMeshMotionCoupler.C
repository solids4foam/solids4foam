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
#include "meshMotionTotalDispSolid.H"
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
    )
{
    const solidModel& model = motionSolver_.model();

    // The partitioned coupling sets the total interface displacement from
    // the reference points, which needs a solid model of the total
    // displacement on the current mesh, solved with either solution
    // algorithm. The quasi-monolithic coupling, which does not use this
    // coupler, checks its own mesh motion solid model, and needs the PETSc
    // SNES algorithm
    if (!isA<solidModels::meshMotionTotalDispSolid>(model))
    {
        FatalErrorInFunction
            << "The " << meshMotionSolidModelFvMotionSolver::typeName
            << " mesh motion solver needs the "
            << solidModels::meshMotionTotalDispSolid::typeName
            << " solid model with this coupling, not " << model.type()
            << exit(FatalError);
    }

    forAll(patchIDs, i)
    {
        const label patchID = patchIDs[i];
        const word& patchName = mesh.boundaryMesh()[patchID].name();

        const fvPatchVectorField& DPatch = model.D().boundaryField()[patchID];

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
                << "The " << model.D().name() << " boundary condition of the "
                << "mesh motion solid model on the patch " << patchName
                << " must be of type fixedDisplacement, with "
                << "enforcePointPatchField no and no displacementSeries"
                << exit(FatalError);
        }

        if
        (
            !isA<fixedValuePointPatchVectorField>
            (
                model.pointD().boundaryField()[patchID]
            )
        )
        {
            FatalErrorInFunction
                << "The " << model.pointD().name() << " boundary condition "
                << "of the mesh motion solid model on the patch " << patchName
                << " must be of type fixedValue" << exit(FatalError);
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

    const vectorField totalDisp
    (
        displacementFluidMeshMotionCoupler::totalPatchDisplacement
        (
            mesh().boundaryMesh()[patchID],
            motionSolver_.points0(),
            iterDisplacement
        )
    );

    fixedValuePointPatchVectorField& pointDPatch =
        refCast<fixedValuePointPatchVectorField>
        (
            boundaryFieldRef(model.pointD())[patchID]
        );

    pointDPatch == totalDisp;

    fixedDisplacementFvPatchVectorField& DPatch =
        refCast<fixedDisplacementFvPatchVectorField>
        (
            boundaryFieldRef(model.D())[patchID]
        );

    primitivePatchInterpolation interp(mesh().boundaryMesh()[patchID]);

    DPatch.totalDisp() = interp.pointToFaceInterpolate(totalDisp);
}

#endif // OPENFOAM_COM

// ************************************************************************* //
