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

#include "fluidMeshMotionCoupler.H"
#include "pointFields.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "dynamicMotionSolverFvMesh.H"
    #include "displacementMotionSolver.H"
#else
    #include "tetPointFields.H"
    #include "newSubsetMotionSolverFvMesh.H"
#endif
#if !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)
    #include "RBFMeshMotionSolver.H"
#endif
#ifdef OPENFOAM_COM
    #include "meshMotionSolidModelFvMotionSolver.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidMeshMotionCoupler, 0);
    defineRunTimeSelectionTable(fluidMeshMotionCoupler, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidMeshMotionCoupler::fluidMeshMotionCoupler
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
:
    mesh_(mesh),
    patchIDs_(patchIDs),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionSolver* Foam::fluidMeshMotionCoupler::lookupMotionSolver
(
    const dynamicFvMesh& mesh
)
{
#ifdef OPENFOAM_NOT_EXTEND
    if (isA<dynamicMotionSolverFvMesh>(mesh))
    {
        return
            &const_cast<motionSolver&>
            (
                refCast<const dynamicMotionSolverFvMesh>(mesh).motion()
            );
    }
#else
    if (mesh.foundObject<motionSolver>("dynamicMeshDict"))
    {
        return
            &const_cast<motionSolver&>
            (
                mesh.lookupObject<motionSolver>("dynamicMeshDict")
            );
    }
#endif

    return nullptr;
}


Foam::word Foam::fluidMeshMotionCoupler::autoSelectType
(
    const dynamicFvMesh& mesh
)
{
    // The order of the checks is the order of precedence

#ifndef OPENFOAM_NOT_EXTEND
    if (mesh.foundObject<tetPointVectorField>("motionU"))
    {
        return "tetFemVelocity";
    }
#endif

    if (mesh.foundObject<pointVectorField>("pointMotionU"))
    {
        return "velocity";
    }

#ifndef OPENFOAM_NOT_EXTEND
    if (isA<newSubsetMotionSolverFvMesh>(mesh))
    {
        return "subsetMesh";
    }
#endif

    const motionSolver* motionPtr = lookupMotionSolver(mesh);

    if (motionPtr)
    {
#if !defined(OPENFOAM_ORG) && !defined(S4F_NO_RBF)
        if (isA<RBFMeshMotionSolver>(*motionPtr))
        {
            return "RBF";
        }
#endif

#ifdef OPENFOAM_COM
        if (isA<meshMotionSolidModelFvMotionSolver>(*motionPtr))
        {
            return "solidModel";
        }
#endif

#ifdef OPENFOAM_NOT_EXTEND
        if (isA<displacementMotionSolver>(*motionPtr))
        {
            return "displacement";
        }
#endif
    }

    FatalErrorInFunction
        << "No fluid mesh motion coupler is available for the fluid mesh "
        << "of type " << mesh.type();

    if (motionPtr)
    {
        FatalError
            << " with mesh motion solver of type " << motionPtr->type();
    }

    FatalError
        << nl << "The automatically selectable mesh motion solvers are "
        << "velocity-based solvers with a pointMotionU field (e.g. "
        << "velocityLaplacian), displacement-based solvers with a "
        << "pointDisplacement field (e.g. displacementLaplacian), the "
        << "RBF solver and solids4foamSolidModel; component solvers such as "
        << "velocityComponentLaplacian are not supported." << nl
        << "Valid fluidMeshMotion types are:" << nl
        << dictionaryConstructorTablePtr_->toc()
        << exit(FatalError);

    return word::null;
}


// ************************************************************************* //
