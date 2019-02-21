/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "totalVelocityLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
// #include "volPointInterpolation.H"
#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(totalVelocityLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        totalVelocityLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalVelocityLaplacianFvMotionSolver::
totalVelocityLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    fvMotionSolver(mesh),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(fvMesh_)
    ),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            vector::zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(*this, lookup("diffusivity"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::totalVelocityLaplacianFvMotionSolver::
~totalVelocityLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::totalVelocityLaplacianFvMotionSolver::curPoints() const
{
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


void Foam::totalVelocityLaplacianFvMotionSolver::solve()
{
    pointIOField refAllPoints
    (
        IOobject
        (
            "points",
            fvMesh_.time().constant(),
            polyMesh::meshSubDir,
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fvMesh& mesh = const_cast<fvMesh&>(fvMesh_);

    vectorField oldAllPoints = fvMesh_.allPoints();

    mesh.movePoints(refAllPoints);

    diffusivityPtr_->correct();
//     pointMotionU_.boundaryField().updateCoeffs(); // Problem at symmetry ZT

    label nNonOrthCorr = 1;
    if (found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(lookup("nNonOrthogonalCorrectors"));
    }

    for (label i=0; i<nNonOrthCorr; i++)
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

    {
        leastSquaresVolPointInterpolation::New(fvMesh_).interpolate
        (
            cellMotionU_,
            pointMotionU_
        );
    }

    mesh.movePoints(oldAllPoints);
}


void Foam::totalVelocityLaplacianFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    fvMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(NULL);
    diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
