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

Foam::word Foam::meshMotionSolidModelFvMotionSolver::regionName() const
{
    const word regionName
    (
        coeffDict().lookupOrDefault<word>("regionName", fvMesh_.name())
    );

    if (regionName != fvMesh_.name())
    {
        FatalErrorInFunction
            << "The solid model is solved on the moving mesh " << fvMesh_.name()
            << ", so regionName, if given, must be " << fvMesh_.name()
            << ", not " << regionName << exit(FatalError);
    }

    return regionName;
}


Foam::autoPtr<Foam::volScalarField>
Foam::meshMotionSolidModelFvMotionSolver::makeStiffnessScaleFactor()
{
    if (!coeffDict().found("diffusivity"))
    {
        return autoPtr<volScalarField>();
    }

    diffusivityPtr_ =
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"));

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
                Switch(coeffDict().lookupOrDefault<Switch>
                (
                    "writeStiffScaleFactor", false
                ))
              ? IOobject::AUTO_WRITE
              : IOobject::NO_WRITE
            ),
            fvMesh_,
            dimensionedScalar("one", dimless, 1),
            "zeroGradient"
        )
    );

    return dPtr;
}


void Foam::meshMotionSolidModelFvMotionSolver::updateStiffnessScaleFactor()
{
    if (!stiffnessScaleFactorPtr_.valid())
    {
        return;
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

    diffusivityPtr_->correct();

    // Limit the range of the face diffusivity and normalise it by its average
    surfaceScalarField df("stiffnessScaleFactorf", diffusivityPtr_()());
    const dimensionedScalar avDf = average(df);

    if (avDf.value() < SMALL)
    {
        FatalErrorInFunction
            << "The average mesh motion diffusivity must be positive"
            << exit(FatalError);
    }

    df = min(maxFactor*avDf, max(minFactor*avDf, df));
    df /= avDf;

    volScalarField& d = stiffnessScaleFactorPtr_();
    d = fvc::average(df);
    d.correctBoundaryConditions();

    Info<< "Stiffness scale factor: min " << gMin(d.primitiveField())
        << ", max " << gMax(d.primitiveField()) << endl;
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
    diffusivityPtr_(),
    stiffnessScaleFactorPtr_(makeStiffnessScaleFactor()),
    modelPtr_(),
    solveEquations_(true)
{
    // The stiffness scale factor is needed when the solid model is created
    updateStiffnessScaleFactor();

    modelPtr_ = solidModel::New
    (
        const_cast<Time&>(fvMesh_.time()), regionName()
    );
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
    diffusivityPtr_(),
    stiffnessScaleFactorPtr_(makeStiffnessScaleFactor()),
    modelPtr_(),
    solveEquations_(true)
{
    // The stiffness scale factor is needed when the solid model is created
    updateStiffnessScaleFactor();

    modelPtr_ = solidModel::New
    (
        const_cast<Time&>(fvMesh_.time()), regionName()
    );

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
        updateStiffnessScaleFactor();

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
