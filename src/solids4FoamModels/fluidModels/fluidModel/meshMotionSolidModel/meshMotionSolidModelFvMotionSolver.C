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

Foam::autoPtr<Foam::dynamicFvMesh>
Foam::meshMotionSolidModelFvMotionSolver::makeReferenceMesh() const
{
    const word regionName(coeffDict().lookup("regionName"));

    if (regionName == fvMesh_.name())
    {
        return autoPtr<dynamicFvMesh>();
    }

    Info<< "Reading the mesh motion reference mesh for region "
        << regionName << endl;

    return dynamicFvMesh::New
    (
        IOobject
        (
            regionName,
            fvMesh_.time().timeName(),
            fvMesh_.time(),
            IOobject::MUST_READ
        )
    );
}


Foam::autoPtr<Foam::volScalarField>
Foam::meshMotionSolidModelFvMotionSolver::makeStiffnessScaleFactor() const
{
    if (!referenceMeshPtr_.valid() || !coeffDict().found("diffusivity"))
    {
        return autoPtr<volScalarField>();
    }

    const fvMesh& mesh = referenceMeshPtr_();

    const scalar maxFactor =
        coeffDict().lookupOrDefault<scalar>("maxFactor", 100);
    const scalar minFactor =
        coeffDict().lookupOrDefault<scalar>("minFactor", 0.1);

    if (maxFactor < minFactor)
    {
        FatalErrorInFunction
            << "maxFactor cannot be less than minFactor!" << exit(FatalError);
    }

    autoPtr<motionDiffusivity> diffusivityPtr
    (
        motionDiffusivity::New(mesh, coeffDict().lookup("diffusivity"))
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
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                Switch(coeffDict().lookupOrDefault<Switch>
                (
                    "writeStiffScaleFactor", false
                ))
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
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


void Foam::meshMotionSolidModelFvMotionSolver::checkReferenceMesh() const
{
    if (!referenceMeshPtr_.valid())
    {
        return;
    }

    const fvMesh& mesh = referenceMeshPtr_();

    bool ok = (mesh.nPoints() == points0().size());
    reduce(ok, andOp<bool>());

    if (ok)
    {
        const scalar maxDiff = gMax(mag(mesh.points() - points0())());

        ok = (maxDiff <= 1e-6*mag(mesh.bounds().span()));
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "The points of the mesh motion reference mesh "
            << mesh.name() << " are not the reference points of the mesh "
            << fvMesh_.name() << ": the reference mesh must be a copy of "
            << "the initial mesh, decomposed in the same way in parallel"
            << exit(FatalError);
    }
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
    referenceMeshPtr_(makeReferenceMesh()),
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
    checkReferenceMesh();
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
    referenceMeshPtr_(makeReferenceMesh()),
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
    checkReferenceMesh();

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
