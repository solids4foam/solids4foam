// Diagnostic study only (verification/hronturek-replay-variants).
//
// dynamicMotionSolverFvMesh whose mesh-motion flux is made consistent with
// the backward (BDF2) ddt scheme: with the Euler swept-volume flux
// phiE^n = (V^n - V^(n-1))/dt (exactly the geometric conservation law for
// Euler), the flux
//
//     meshPhi^n = 1.5 phiE^n - 0.5 phiE^(n-1)
//
// has div(meshPhi^n) = (3 V^n - 4 V^(n-1) + V^(n-2))/(2 dt), which is the
// volume derivative in the backward ddt term. The standard mesh flux is
// phiE^n, so the backward scheme violates the GCL by 0.5 (V^n - 2V^(n-1)
// + V^(n-2))/dt = O(dt). At the first update after a restart there is no
// old flux and the Euler flux is kept.
//
// Use with dynamicFvMesh dynamicMotionSolverBDF2FvMesh; and
// libs ("libhtVariants.so"); in controlDict. Constant time step only.

#ifndef dynamicMotionSolverBDF2FvMesh_H
#define dynamicMotionSolverBDF2FvMesh_H

#include "dynamicMotionSolverFvMesh.H"
#include "surfaceFields.H"

namespace Foam
{

class dynamicMotionSolverBDF2FvMesh
:
    public dynamicMotionSolverFvMesh
{
    //- Time index of the last update
    label updateTimeIndex_;

    //- Euler flux of the current time step (latest update)
    autoPtr<surfaceScalarField> phiEulerPtr_;

    //- Euler flux of the previous time step
    autoPtr<surfaceScalarField> phiEulerOldPtr_;

    dynamicMotionSolverBDF2FvMesh(const dynamicMotionSolverBDF2FvMesh&) = delete;
    void operator=(const dynamicMotionSolverBDF2FvMesh&) = delete;

public:

    TypeName("dynamicMotionSolverBDF2FvMesh");

    dynamicMotionSolverBDF2FvMesh(const IOobject& io, const bool doInit = true);

    virtual ~dynamicMotionSolverBDF2FvMesh() = default;

    virtual bool update();
};

} // End namespace Foam

#endif
