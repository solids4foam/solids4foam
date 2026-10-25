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

#include "meshMotionSolidModelFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "fvcAverage.H"
#include "mapPolyMesh.H"
#include "fvOptions.H"
#include "motionInterpolation.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshMotionSolidModelFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        meshMotionSolidModelFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        meshMotionSolidModelFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::volScalarField>
Foam::meshMotionSolidModelFvMotionSolver::makeStiffnessScaleFactor() const
{
    if (!coeffDict().found("diffusivity"))
    {
        return autoPtr<volScalarField>();
    }

    const scalar maxFactor =
        coeffDict().lookupOrDefault<scalar>("maxFactor", 100);
    const scalar minFactor =
        coeffDict().lookupOrDefault<scalar>("minFactor", 0.1);

    if (maxFactor < minFactor)
    {
        FatalErrorInFunction
            << "maxFactor cannot be less than minFactor!" << exit(FatalError);
    }

    const IOobject::writeOption writeOpt =
        Switch
        (
            coeffDict().lookupOrDefault<Switch>("writeStiffScaleFactor", false)
        )
      ? IOobject::AUTO_WRITE
      : IOobject::NO_WRITE;

    // The factor is calculated once, from the mesh at the start of the run.
    // On a restart the mesh is the deformed mesh of the start time, so the
    // factor of the original run is read if it was written (it is written
    // when the mesh motion solid has 'restart yes'); otherwise the mesh motion
    // after the restart is not that of a continuous run
    if (fvMesh_.time().startTimeIndex() > 0)
    {
        IOobject io
        (
            "stiffnessScaleFactor",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            writeOpt
        );

        if (io.typeHeaderOk<volScalarField>(true))
        {
            autoPtr<volScalarField> dPtr(new volScalarField(io, fvMesh_));

            Info<< "Stiffness scale factor read from "
                << fvMesh_.time().timeName() << ": min "
                << gMin(dPtr->primitiveField()) << ", max "
                << gMax(dPtr->primitiveField()) << endl;

            return dPtr;
        }

        WarningInFunction
            << "The stiffness scale factor was not found at the start time "
            << fvMesh_.time().timeName() << ", so it is calculated from "
            << "the current (deformed) mesh rather than the mesh at the start "
            << "of the original run: the mesh motion is not that of a "
            << "continuous run. Set 'restart yes' in the mesh motion solid "
            << "model coefficients to write it" << endl;
    }

    autoPtr<motionDiffusivity> diffusivityPtr
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    );

    // Limit the range of the face diffusivity and normalise it by its average
    surfaceScalarField df("stiffnessScaleFactorf", diffusivityPtr()());
    const dimensionedScalar avDf = average(df);

    if (avDf.value() < SMALL)
    {
        FatalErrorInFunction
            << "The average mesh motion diffusivity must be positive"
            << exit(FatalError);
    }

    df = min(maxFactor*avDf, max(minFactor*avDf, df));
    df /= avDf;

    autoPtr<volScalarField> dPtr
    (
        new volScalarField
        (
            IOobject
            (
                "stiffnessScaleFactor",
                fvMesh_.time().timeName(),
                fvMesh_,
                IOobject::NO_READ,
                writeOpt
            ),
            fvc::average(df),
            "zeroGradient"
        )
    );

    dPtr->correctBoundaryConditions();

    Info<< "Stiffness scale factor: min " << gMin(dPtr->primitiveField())
        << ", max " << gMax(dPtr->primitiveField()) << endl;

    return dPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshMotionSolidModelFvMotionSolver::meshMotionSolidModelFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    stiffnessScaleFactorPtr_(makeStiffnessScaleFactor()),
    modelPtr_
    (
        solidModel::New
        (
            const_cast<Time&>(fvMesh_.time()),
            word(coeffDict().lookup("regionName"))
        )
    ),
    solveEquations_(true)
{
    // A restart needs the stiffness scale factor of the original run
    if (stiffnessScaleFactorPtr_.valid() && modelPtr_->restart())
    {
        stiffnessScaleFactorPtr_->writeOpt() = IOobject::AUTO_WRITE;
    }
}


Foam::meshMotionSolidModelFvMotionSolver::meshMotionSolidModelFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    stiffnessScaleFactorPtr_(makeStiffnessScaleFactor()),
    modelPtr_
    (
        solidModel::New
        (
            const_cast<Time&>(fvMesh_.time()),
            word(coeffDict().lookup("regionName"))
        )
    ),
    solveEquations_(true)
{
    // A restart needs the stiffness scale factor of the original run
    if (stiffnessScaleFactorPtr_.valid() && modelPtr_->restart())
    {
        stiffnessScaleFactorPtr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    // Update pointD
    modelPtr_->pointD().primitiveFieldRef() =
        pointDisplacement.primitiveField();
    modelPtr_->pointD().correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshMotionSolidModelFvMotionSolver::
~meshMotionSolidModelFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::meshMotionSolidModelFvMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints
    (
        points0() + modelPtr_->pointD().primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    if (debug && Pstream::parRun())
    {
        // Check that the processors agree on the position of every point they
        // share; otherwise the moved mesh does not match across processor
        // boundaries
        pointField maxPoints(tcurPoints());
        pointField minPoints(tcurPoints());
        pointConstraints::syncUntransformedData
        (
            fvMesh_, maxPoints, maxEqOp<vector>()
        );
        pointConstraints::syncUntransformedData
        (
            fvMesh_, minPoints, minEqOp<vector>()
        );

        const scalarField mismatch(mag(maxPoints - minPoints));

        label nMismatch = 0;
        forAll(mismatch, pointI)
        {
            if (mismatch[pointI] > 1e-12)
            {
                ++nMismatch;
            }
        }

        Info<< "Processor point mismatch of the moved mesh: max "
            << gMax(mismatch) << " m, number of point copies above 1e-12 m "
            << returnReduce(nMismatch, sumOp<label>()) << endl;
    }

    pointDisplacement_.primitiveFieldRef() =
        modelPtr_->pointD().primitiveField();

    return tcurPoints;
}


void Foam::meshMotionSolidModelFvMotionSolver::solve()
{
    Info<< "Updating the mesh" << endl;

    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    // Evolve the solid model
    if (solveEquations_)
    {
        modelPtr_->evolve();
    }
}


void Foam::meshMotionSolidModelFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
