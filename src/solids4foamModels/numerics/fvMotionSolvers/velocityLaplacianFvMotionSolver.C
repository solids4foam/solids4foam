/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "velocityLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityLaplacianFvMotionSolver::curPoints() const
{
    leastSquaresVolPointInterpolation::New(fvMesh_).interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    tmp<pointField> tcurPoints(new pointField(fvMesh_.allPoints()));
    pointField& cp = tcurPoints();
    const pointField& pointMotionUI = pointMotionU_.internalField();

    forAll(pointMotionUI, pointI)
    {
        cp[pointI] += 
            pointMotionUI[pointI]*fvMesh_.time().deltaT().value();
    }
    
//     tmp<pointField> tcurPoints
//     (
//         fvMesh_.points()
//       + fvMesh_.time().deltaT().value()*pointMotionU_.internalField()
//     );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::velocityLaplacianFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.allPoints());
//     movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryField().updateCoeffs();

    for (label i=0; i<2; i++)
    {
        Foam::solve
        (
            fvm::laplacian
            (
                diffusivityPtr_->operator()(),
                cellMotionU_,
                "laplacian(diffusivity,cellMotionU)"
            )
        );
    }
}


// ************************************************************************* //
