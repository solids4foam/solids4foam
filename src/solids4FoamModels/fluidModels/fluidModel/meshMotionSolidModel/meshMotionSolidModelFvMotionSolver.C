/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meshMotionSolidModelFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"
#include "fvOptions.H"

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshMotionSolidModelFvMotionSolver::meshMotionSolidModelFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    modelPtr_
    (
        solidModel::New
        (
            const_cast<Time&>(fvMesh_.time()),
            word(coeffDict().lookup("regionName"))
        )
    )
{}


Foam::meshMotionSolidModelFvMotionSolver::
meshMotionSolidModelFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    modelPtr_
    (
        solidModel::New
        (
            const_cast<Time&>(fvMesh_.time()),
            word(coeffDict().lookup("regionName"))
        )
    )
{
    // Update pointD
    modelPtr_->pointD().primitiveFieldRef() =
        pointDisplacement.primitiveField();
    modelPtr_->pointD().correctBoundaryConditions();

    FatalError
        << "stop " << __LINE__ << abort(FatalError);
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
    modelPtr_->evolve();
}


void Foam::meshMotionSolidModelFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
