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

#include "triSurfaceMotion.H"
#include "RodriguesRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMotion, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMotion::triSurfaceMotion
(
    const Time& runTime,
    const dictionary& dict,
    const standAlonePatch& patch
)
:
    runTime_(runTime),
    dict_(dict),
    patch_(patch),
    displacement_(dict_.lookupOrDefault<vector>("displacement", vector::zero)),
    dispRampTime_(dict_.lookupOrDefault<scalar>("dispRampTime", 0.0)),
    rpm_(dict_.lookupOrDefault<scalar>("rpm", 0.0)),
    rotationRampTime_(dict_.lookupOrDefault<scalar>("rotationRampTime", 0.0)),
    rotationAxis_(dict_.lookupOrDefault<vector>("rotationAxis", vector(1,0,0))),
    currentAxisDisplacement_(vector::zero),
    currentAngle_(0.0),
    currentRpm_(0.0),
    origPoints_(patch.points()),
    origFaceCentres_(patch.faceCentres()),
    oldPointDisp_(patch.nPoints(), vector::zero),
    oldFaceDisp_(patch.size(), vector::zero),
    initialRotationOrigin_(vector::zero),
    faceDisplacementIncrement_(patch.size(), vector::zero)
{
    // Calculate initial rotation origin
    // standAlonePatch does not store face areas so we will calculate them
    scalarField faceAreas(patch_.size(), 0.0);
    forAll(faceAreas, faceI)
    {
        faceAreas[faceI] = mag(patch_[faceI].normal(patch_.points()));
    }

    initialRotationOrigin_ =
        gSum(patch_.faceCentres()*mag(faceAreas))
       /(gSum(mag(faceAreas)) + SMALL);

    // CHECK RESTART!!!!
    WarningIn("triSurfaceMotion::triSurfaceMotion(...)")
        << "Restart has yet to be implemented!" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::triSurfaceMotion::~triSurfaceMotion()
{}


// * * * * * * * * * * * * * * * Public Member Functions  * * *  * * * * * * //


Foam::tmp<Foam::pointField> Foam::triSurfaceMotion::newPoints()
{
    tmp<pointField> tnewPoints(new pointField(patch_.points()));
    pointField& newPoints = tnewPoints();

    // Current time
    const scalar curTime = runTime_.value();

    // Calculate  point and face displacement
    vectorField pointDisp(newPoints.size(), vector::zero);
    faceDisplacementIncrement_ = vector::zero;

    // Translation

    if (curTime <= (dispRampTime_ + SMALL))
    {
        Info<< "triSurfaceMotion: surface translation phase" << endl;

        // Linearly ramp displacement
        pointDisp = displacement_*curTime/dispRampTime_;
        faceDisplacementIncrement_ = displacement_*curTime/dispRampTime_;
        currentAxisDisplacement_ = displacement_*curTime/dispRampTime_;
    }
    else
    {
        pointDisp = displacement_;
        faceDisplacementIncrement_ = displacement_;
        currentAxisDisplacement_ = displacement_;
    }

    // Rotation

    // Analytically calculate current angle
    if (curTime <= (rotationRampTime_ + SMALL))
    {
        Info<< "triSurfaceMotion: rotational acceleration phase" << endl;

        // Roller acceleration
        currentAngle_ = Foam::pow(curTime, 2)*rpm_*3.0/(rotationRampTime_);

        // Current RPM
        currentRpm_ = (curTime/rotationRampTime_)*rpm_;
    }
    else
    {
        Info<< "triSurfaceMotion: constant rotation phase" << endl;

        // Constant angular velocity
        currentAngle_ =
            rotationRampTime_*rpm_*3.0 // acceleration
            + (curTime - rotationRampTime_)*rpm_*6.0; // constant velocity

        // Current RPM
        currentRpm_ = rpm_;
    }

    // Create rotation matrix
    const tensor rotMat = RodriguesRotation(rotationAxis_, currentAngle_);

    // Calculate position of rotated points
    const pointField rotatedPoints =
        (rotMat & (origPoints_ - initialRotationOrigin_))
        + initialRotationOrigin_;

    // Calculate position of rotated faces
    const pointField rotatedFaces =
        (rotMat & (origFaceCentres_ - initialRotationOrigin_))
        + initialRotationOrigin_;

    // Add total displacement due to rotation to the displacement due to
    // translation
    pointDisp += rotatedPoints - origPoints_;
    faceDisplacementIncrement_ += rotatedFaces - origFaceCentres_;

    // Subtract old displacement as we are setting the displacement increment
    pointDisp -= oldPointDisp_;
    faceDisplacementIncrement_ -= oldFaceDisp_;

    // Update old displacement
    oldPointDisp_ += pointDisp;
    oldFaceDisp_ += faceDisplacementIncrement_;

    // Update the new points
    newPoints += pointDisp;

    return tnewPoints;
}


// ************************************************************************* //
