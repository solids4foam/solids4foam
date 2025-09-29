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

#include "decayingTaylorGreenVortexMeshMotion.H"
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
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decayingTaylorGreenVortexMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        decayingTaylorGreenVortexMeshMotion,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        decayingTaylorGreenVortexMeshMotion,
        displacement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decayingTaylorGreenVortexMeshMotion::decayingTaylorGreenVortexMeshMotion
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    A_
    (
        readScalar
        (
            dict.subDict
            (
                "decayingTaylorGreenVortexMeshMotionCoeffs"
            ).lookup("scaleFactor")
        )
    ),
    orthogonalMeshMotion_
    (
        dict.subDict
        (
            "decayingTaylorGreenVortexMeshMotionCoeffs"
        ).lookup("orthogonalMeshMotion")
    )
{}


Foam::decayingTaylorGreenVortexMeshMotion::decayingTaylorGreenVortexMeshMotion
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    A_
    (
        readScalar
        (
            dict.subDict
            (
                "decayingTaylorGreenVortexMeshMotionCoeffs"
            ).lookup("scaleFactor")
        )
    ),
    orthogonalMeshMotion_
    (
        dict.subDict
        (
            "decayingTaylorGreenVortexMeshMotionCoeffs"
        ).lookup("orthogonalMeshMotion")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decayingTaylorGreenVortexMeshMotion::
~decayingTaylorGreenVortexMeshMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::decayingTaylorGreenVortexMeshMotion::curPoints() const
{
    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::decayingTaylorGreenVortexMeshMotion::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    // Take a reference to the initial points
    const vectorField& points0 = this->points0();

    // Calculate Prepare x,y and t
    const scalar pi = constant::mathematical::pi;
    const scalarField x(points0.component(vector::X));
    const scalarField y(points0.component(vector::Y));
    const scalar t = time().value();

    // Calculate the displacement
    if (orthogonalMeshMotion_)
    {
        const scalar timeFactor = A_*Foam::sin(pi*t/0.8);

        pointDisplacement_.primitiveFieldRef() =
            timeFactor*Foam::sin(pi*x)*vector(1, 0, 0)
          + timeFactor*Foam::sin(pi*y)*vector(0, 1, 0);
    }
    else
    {
        pointDisplacement_.primitiveFieldRef() =
            A_*Foam::sin(pi*t/0.8)*Foam::sin(pi*x)
           *Foam::sin(pi*y)*vector(1, 1, 0);
    }
}


void Foam::decayingTaylorGreenVortexMeshMotion::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
