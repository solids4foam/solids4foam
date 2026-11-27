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

Class
    Foam::movingVentricleFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for a ventricle whose cavity volume follows a
    prescribed curve.

    The target volume at time t is the volume of the undeformed mesh plus the
    integral of a prescribed wall flux, timeVsFlowSeries (m3/s against s).
    Each point of the wall patch moves along its initial unit normal by

        d = a max(0, y0/yMin) n0

    where y0 is the initial y coordinate of the point and yMin the smallest
    on the patch, so the motion is largest at the apex (the most negative y)
    and vanishes for y0 >= 0. The scale factor a is found by a Newton loop
    (forward-difference slope, damped) that meets the target volume to relTol
    times the target. The interior points follow by a Laplacian of the cell
    displacement, as in displacementLaplacian.

    With smoothFlow the wall flux is a natural cubic spline through the
    samples and the volume its exact integral, so the volume is smooth in
    time. Otherwise the flux samples are integrated by the trapezoidal rule
    and the volume is interpolated linearly, which makes the wall flux
    piecewise constant and limits a time-step convergence study to first
    order.

Usage
    In constant/dynamicMeshDict:
    \verbatim
    dynamicFvMesh       dynamicMotionSolverFvMesh;
    motionSolverLibs    ("libfvMotionSolvers.so" "libcirculationMeshMotion.so");
    motionSolver        movingVentricle;

    movingVentricleCoeffs
    {
        patchName       WALL;       // the moving wall
        relTol          1e-12;      // volume tolerance / target (1e-6)
        smoothFlow      true;       // spline wall flux (false)
        timeVsFlowSeries
        {
            file        "$FOAM_CASE/constant/timeVsFlow";
            outOfBounds clamp;
        }
        diffusivity     quadratic inverseDistance (WALL);
    }
    \endverbatim

    The volume reached at each step is printed as "Converged: volume = ".

SourceFiles
    movingVentricleFvMotionSolver.C

Author
    Philip Cardiff, UCD.
    Chanikya Valeti, UCD.
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef movingVentricleFvMotionSolver_H
#define movingVentricleFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class movingVentricleFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class movingVentricleFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;

        //- Index of the moving wall patch
        const label patchID_;

        //- Volume of the undeformed mesh, points0
        const scalar volumeMin_;

        //- Smallest initial y coordinate on the wall patch
        scalar minYPoint_;

        //- Tolerance on the volume, as a fraction of the target volume
        const scalar relTol_;

        //- Interpolate the wall flux with a natural cubic spline and take the
        //  volume as its exact integral
        const Switch smoothFlow_;

        //- Times of the wall flux samples
        scalarField times_;

        //- Wall flux samples
        scalarField flowRates_;

        //- Target volumes at the sample times
        scalarField volumes_;

        //- Second derivatives of the wall flux spline at the samples
        scalarField splineM_;

        //- Initial wall patch points
        vectorField initialPatchPoints_;

        //- Initial wall patch point unit normals
        vectorField initialPatchPointNormals_;


    // Private Member Functions

        //- Read the wall flux series and set the target volumes and the
        //  initial wall geometry
        void initialise();

        //- Volume of the undeformed mesh, points0, by the divergence theorem
        //  over the boundary faces, without moving the mesh
        scalar points0Volume() const;

        //- Prescribed ventricle volume at time t
        scalar targetVolume(const scalar t) const;

        //- Linear interpolation of y(x) at xq, clamped at the ends
        scalar linearInterp
        (
            const scalarField& x,
            const scalarField& y,
            const scalar xq
        ) const;

        //- Wall displacement for the scale factor a
        tmp<vectorField> wallDisplacement(const scalar a) const;

        //- Set the wall displacement for the scale factor a on the
        //  pointDisplacement boundary
        void setWallDisplacement(const scalar a);

        //- Set the wall displacement for the scale factor a and return the
        //  mesh volume it gives
        scalar wallVolume(const scalar a);

        //- Mesh volume with the wall points displaced from points0
        scalar calculateVolume(const vectorField& displacement) const;

        //- No copy construct
        movingVentricleFvMotionSolver
        (
            const movingVentricleFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const movingVentricleFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("movingVentricle");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        movingVentricleFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );

        //- Construct from components
        movingVentricleFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    virtual ~movingVentricleFvMotionSolver() = default;


    // Member Functions

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
