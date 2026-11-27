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

#include "movingVentricleFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "interpolationTable.H"
#include "PatchTools.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingVentricleFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        movingVentricleFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        movingVentricleFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::movingVentricleFvMotionSolver::initialise()
{
    if (patchID_ == -1)
    {
        FatalErrorInFunction
            << "Did not find the patch "
            << coeffDict().get<word>("patchName") << exit(FatalError);
    }

    Info<< typeName << ": minimum volume = " << volumeMin_
        << " (undeformed mesh; current mesh "
        << gSum(fvMesh_.cellVolumes()) << ")" << endl;

    // Read the wall flux series
    const interpolationTable<scalar> flowTable
    (
        coeffDict().subDict("timeVsFlowSeries")
    );

    times_.setSize(flowTable.size());
    flowRates_.setSize(flowTable.size());

    forAll(flowTable, i)
    {
        times_[i] = flowTable[i].first();
        flowRates_[i] = flowTable[i].second();
    }

    // Integrate the flow rate with the trapezoidal rule to get the volume
    // against time
    volumes_.setSize(times_.size(), 0.0);
    for (label i = 1; i < times_.size(); ++i)
    {
        const scalar dt = times_[i] - times_[i - 1];
        const scalar avgFlow = 0.5*(flowRates_[i] + flowRates_[i - 1]);
        volumes_[i] = volumes_[i - 1] + avgFlow*dt;
    }

    if (smoothFlow_)
    {
        // Natural cubic spline through the flow samples: solve for the
        // second derivatives M (M = 0 at both ends) with the Thomas algorithm
        const label n = times_.size();
        splineM_.setSize(n, 0.0);

        if (n > 2)
        {
            scalarField c(n, 0.0);
            scalarField d(n, 0.0);

            for (label i = 1; i < n - 1; ++i)
            {
                const scalar h0 = times_[i] - times_[i - 1];
                const scalar h1 = times_[i + 1] - times_[i];
                const scalar rhs =
                    6
                   *(
                        (flowRates_[i + 1] - flowRates_[i])/h1
                      - (flowRates_[i] - flowRates_[i - 1])/h0
                    );

                // Row: h0 M[i-1] + 2(h0 + h1) M[i] + h1 M[i+1] = rhs
                const scalar diag = 2*(h0 + h1) - h0*c[i - 1];
                c[i] = h1/diag;
                d[i] = (rhs - h0*d[i - 1])/diag;
            }

            for (label i = n - 2; i > 0; --i)
            {
                splineM_[i] = d[i] - c[i]*splineM_[i + 1];
            }
        }

        // The exact integral of the spline over each interval
        for (label i = 1; i < n; ++i)
        {
            const scalar h = times_[i] - times_[i - 1];
            volumes_[i] =
                volumes_[i - 1]
              + 0.5*h*(flowRates_[i - 1] + flowRates_[i])
              - h*h*h*(splineM_[i - 1] + splineM_[i])/24;
        }
    }

    // Offset by the undeformed volume
    volumes_ += volumeMin_;

    Info<< typeName << ": using timeVsFlowSeries with "
        << times_.size() << " samples, volumeMin = " << volumeMin_
        << (smoothFlow_ ? ", flow as a cubic spline (smoothFlow)" : "")
        << endl;

    const polyPatch& wall = fvMesh_.boundaryMesh()[patchID_];

    minYPoint_ = gMin(wall.localPoints().component(vector::Y));

    initialPatchPoints_ = wall.localPoints();

    // Parallel-consistent point normals
    initialPatchPointNormals_ = PatchTools::pointNormals(fvMesh_, wall);
}


Foam::scalar Foam::movingVentricleFvMotionSolver::points0Volume() const
{
    // The prescribed volume is the volume of the undeformed mesh plus the
    // integral of the flow, so the reference must be the undeformed volume:
    // on a restart the mesh the solver is constructed on is already deformed.
    //
    // V = (1/3) sum over boundary faces of (x_f . S_f): with OpenFOAM's
    // face-centre pyramid decomposition this equals the sum of the cell
    // volumes exactly. Processor faces are seen from both sides with
    // opposite normals and cancel in the global sum.
    const pointField& p0 = points0();
    const faceList& faces = fvMesh_.faces();

    scalar sum = 0;
    for
    (
        label facei = fvMesh_.nInternalFaces();
        facei < fvMesh_.nFaces();
        ++facei
    )
    {
        const face& f = faces[facei];
        sum += (f.centre(p0) & f.areaNormal(p0));
    }

    return returnReduce(sum, sumOp<scalar>())/3;
}


Foam::scalar Foam::movingVentricleFvMotionSolver::targetVolume
(
    const scalar t
) const
{
    if (!smoothFlow_)
    {
        return linearInterp(times_, volumes_, t);
    }

    // Clamp outside the table, as the linear form does
    if (t <= times_.first())
    {
        return volumes_.first();
    }

    if (t >= times_.last())
    {
        return volumes_.last();
    }

    // Interval [times_[i-1], times_[i]] holding t
    const label i = label
    (
        std::upper_bound(times_.begin(), times_.end(), t) - times_.begin()
    );

    const scalar h = times_[i] - times_[i - 1];
    const scalar s = t - times_[i - 1];
    const scalar q0 = flowRates_[i - 1];
    const scalar q1 = flowRates_[i];
    const scalar M0 = splineM_[i - 1];
    const scalar M1 = splineM_[i];

    // Spline on the interval: q0 + b s + M0 s^2/2 + (M1 - M0) s^3/(6h)
    const scalar b = (q1 - q0)/h - h*(2*M0 + M1)/6;

    return
        volumes_[i - 1]
      + q0*s + 0.5*b*s*s + M0*s*s*s/6 + (M1 - M0)*s*s*s*s/(24*h);
}


Foam::scalar Foam::movingVentricleFvMotionSolver::linearInterp
(
    const scalarField& x,
    const scalarField& y,
    const scalar xq
) const
{
    if (xq <= x.first())
    {
        return y.first();
    }

    if (xq >= x.last())
    {
        return y.last();
    }

    label i = 1;
    while (i < x.size() - 1 && x[i] < xq)
    {
        ++i;
    }

    const scalar w = (xq - x[i - 1])/max(VSMALL, x[i] - x[i - 1]);

    return y[i - 1] + w*(y[i] - y[i - 1]);
}


Foam::tmp<Foam::vectorField>
Foam::movingVentricleFvMotionSolver::wallDisplacement(const scalar a) const
{
    return
        initialPatchPointNormals_
       *max(0.0, initialPatchPoints_.component(vector::Y)/minYPoint_)
       *a;
}


void Foam::movingVentricleFvMotionSolver::setWallDisplacement(const scalar a)
{
    // Push to the boundary and enforce continuity across processors
    const vectorField disp(wallDisplacement(a));
    pointDisplacement_.boundaryFieldRef()[patchID_] == disp;
    pointDisplacement_.correctBoundaryConditions();
}


Foam::scalar Foam::movingVentricleFvMotionSolver::wallVolume(const scalar a)
{
    setWallDisplacement(a);

    return calculateVolume
    (
        pointDisplacement_.boundaryField()[patchID_].patchInternalField()
    );
}


Foam::scalar Foam::movingVentricleFvMotionSolver::calculateVolume
(
    const vectorField& displacement
) const
{
    // Move the mesh to points0 with the wall points displaced, take its
    // volume and move it back. Only the boundary points affect the total
    // volume, so the interior points are left at points0.
    fvMesh& mesh = const_cast<fvMesh&>(fvMesh_);

    const pointField oldPoints(mesh.points());
    const labelList& wallPoints = mesh.boundaryMesh()[patchID_].meshPoints();

    pointField newPoints(points0());
    forAll(wallPoints, pointi)
    {
        newPoints[wallPoints[pointi]] += displacement[pointi];
    }

    mesh.movePoints(newPoints);

    const scalar totalVolume = gSum(mesh.cellVolumes());

    mesh.movePoints(oldPoints);

    return totalVolume;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingVentricleFvMotionSolver::movingVentricleFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    patchID_
    (
        mesh.boundaryMesh().findPatchID(coeffDict().get<word>("patchName"))
    ),
    volumeMin_(points0Volume()),
    minYPoint_(0),
    relTol_(coeffDict().getOrDefault<scalar>("relTol", 1e-6)),
    smoothFlow_(coeffDict().getOrDefault<Switch>("smoothFlow", false))
{
    initialise();
}


Foam::movingVentricleFvMotionSolver::movingVentricleFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    patchID_
    (
        mesh.boundaryMesh().findPatchID(coeffDict().get<word>("patchName"))
    ),
    volumeMin_(points0Volume()),
    minYPoint_(0),
    relTol_(coeffDict().getOrDefault<scalar>("relTol", 1e-6)),
    smoothFlow_(coeffDict().getOrDefault<Switch>("smoothFlow", false))
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::movingVentricleFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_)
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }

    return *diffusivityPtr_;
}


Foam::tmp<Foam::pointField>
Foam::movingVentricleFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );
    pointField& curPoints = tcurPoints.ref();

    // Implement frozen points
    if (frozenPointsZone_ != -1)
    {
        const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

        forAll(pz, i)
        {
            curPoints[pz[i]] = points0()[pz[i]];
        }
    }

    twoDCorrectPoints(curPoints);

    return tcurPoints;
}


void Foam::movingVentricleFvMotionSolver::solve()
{
    // The points have moved so before interpolation update the motionSolver
    // accordingly
    movePoints(fvMesh_.points());

    const scalar targetVolume = this->targetVolume(fvMesh_.time().value());

    // Newton loop on the wall displacement scale factor a, with a
    // forward-difference slope
    scalar a = 0;
    scalar V = wallVolume(a);

    const scalar tol = relTol_*targetVolume;
    const label maxIter = 20;
    label iter = 0;
    scalar res = mag(V - targetVolume);

    Info<< "Find the displacement scale factor using the Newton-Raphson"
        << " method" << nl
        << "    Target volume = " << targetVolume << endl;

    while (res > tol && ++iter < maxIter)
    {
        const scalar h = max(1e-8*mag(a), 1e-8);
        const scalar dVda = (wallVolume(a + h) - V)/h;

        if (mag(dVda) < VSMALL)
        {
            Info<< "    Newton: tiny slope, stopping" << endl;
            break;
        }

        scalar aNew = a - (V - targetVolume)/dVda;

        // Damping: halve the step, at most five times, until the residual
        // falls
        scalar Vnew = wallVolume(aNew);
        label nHalvings = 0;
        while
        (
            mag(Vnew - targetVolume) > 0.9*mag(V - targetVolume)
         && nHalvings++ < 5
        )
        {
            aNew = 0.5*(a + aNew);
            Vnew = wallVolume(aNew);
        }

        a = aNew;
        V = Vnew;
        res = mag(V - targetVolume);

        Info<< "    " << iter
            << ": scale factor = " << a
            << ", volume = " << V
            << ", |res| = " << res << endl;
    }

    if (iter == maxIter)
    {
        FatalErrorInFunction
            << "Max iterations reached in the Newton-Raphson loop"
            << exit(FatalError);
    }

    // Keep the converged wall displacement
    setWallDisplacement(a);

    Info<< "    Converged: volume = " << V << endl;

    // Update the diffusivity field
    diffusivity().correct();

    // Update the pointDisplacement boundary conditions: these will be used
    // by the cellDisplacement boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Solve for cellDisplacement
    Info<< "Solving the mesh motion for cellDisplacement" << endl;
    fvVectorMatrix DEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
    );

    DEqn.solveSegregatedOrCoupled();
}


void Foam::movingVentricleFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
