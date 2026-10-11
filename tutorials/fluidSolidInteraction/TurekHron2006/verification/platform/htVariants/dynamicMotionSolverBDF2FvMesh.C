#include "dynamicMotionSolverBDF2FvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverBDF2FvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverBDF2FvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverBDF2FvMesh,
        doInit
    );
}

Foam::dynamicMotionSolverBDF2FvMesh::dynamicMotionSolverBDF2FvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicMotionSolverFvMesh(io, doInit),
    updateTimeIndex_(-1),
    phiEulerPtr_(),
    phiEulerOldPtr_()
{}


bool Foam::dynamicMotionSolverBDF2FvMesh::update()
{
    // Same as dynamicMotionSolverFvMesh::update() up to the mesh flux
    fvMesh::movePoints
    (
        const_cast<motionSolver&>(motion()).newPoints()
    );

    refPtr<surfaceScalarField> tphi = setPhi();
    if (tphi.valid())
    {
        surfaceScalarField& mphi = tphi.ref();

        if (updateTimeIndex_ != time().timeIndex())
        {
            // New time step: the last Euler flux becomes the old one
            if (phiEulerPtr_.valid())
            {
                phiEulerOldPtr_ = std::move(phiEulerPtr_);
            }
            updateTimeIndex_ = time().timeIndex();
        }

        // Euler (swept volume) flux of this step, as set by movePoints
        phiEulerPtr_.reset
        (
            new surfaceScalarField("meshPhiEuler", mphi)
        );

        if (phiEulerOldPtr_.valid())
        {
            mphi == 1.5*phiEulerPtr_() - 0.5*phiEulerOldPtr_();
        }
    }

    volVectorField* Uptr = getObjectPtr<volVectorField>("U");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();
    }

    return true;
}
