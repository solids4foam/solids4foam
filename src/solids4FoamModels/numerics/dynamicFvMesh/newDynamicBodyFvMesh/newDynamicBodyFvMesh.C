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

#include "newDynamicBodyFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "pointFields.H"
#include "mathematicalConstants.H"
#include "tetMotionSolver.H"
#include "laplaceTetMotionSolver.H"
#include "fixedValueTetPolyPatchFields.H"
#include "transformField.H"
#include "fixedValuePointPatchFields.H"
#include "tetPolyPatchInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newDynamicBodyFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, newDynamicBodyFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newDynamicBodyFvMesh::newDynamicBodyFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false  // Do not register
            )
        ).subDict(typeName + "Coeffs")
    ),
    motionPtr_(motionSolver::New(*this)),
    bodyPatchName_
    (
        dynamicMeshCoeffs_.lookup("bodyPatchName")
    ),
    bodyPatchID_(-1),
    translationDirection_
    (
        dynamicMeshCoeffs_.lookup("translationDirection")
    ),
    translationAmplitude_
    (
        readScalar(dynamicMeshCoeffs_.lookup("translationAmplitude"))
    ),
    translationFrequency_
    (
        readScalar(dynamicMeshCoeffs_.lookup("translationFrequency"))
    ),
    initialRotationOrigin_
    (
        dynamicMeshCoeffs_.lookup("initialRotationOrigin")
    ),
    rotationAxis_
    (
        dynamicMeshCoeffs_.lookup("rotationAxis")
    ),
    rotationAmplitude_
    (
        readScalar(dynamicMeshCoeffs_.lookup("rotationAmplitude"))
    ),
    rotationFrequency_
    (
        readScalar(dynamicMeshCoeffs_.lookup("rotationFrequency"))
    )
{
    bodyPatchID_ = boundaryMesh().findPatchID(bodyPatchName_);

    if (bodyPatchID_<0)
    {
        FatalErrorIn
        (
            "newDynamicBodyFvMesh::newDynamicBodyFvMesh(const IOobject& io)"
        )
            << "Can't find patch: " << bodyPatchName_
                << exit(FatalError);
    }

    translationDirection_ /= mag(translationDirection_) + SMALL;

    rotationAxis_ /= mag(rotationAxis_) + SMALL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newDynamicBodyFvMesh::~newDynamicBodyFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::newDynamicBodyFvMesh::update()
{
    scalar curTime = time().value();
    scalar oldTime = curTime - time().deltaT().value();

    //Info<< "newDynamicBodyFvMesh::update()" << endl;

    {
        vector trans =
            translationAmplitude_
           *(
                sin(2*mathematicalConstant::pi*translationFrequency_*curTime)
              - sin(2*mathematicalConstant::pi*translationFrequency_*oldTime)
            )
           *translationDirection_;


        scalar rotAngle =
            rotationAmplitude_
           *(
                sin(2*mathematicalConstant::pi*rotationFrequency_*curTime)
              - sin(2*mathematicalConstant::pi*rotationFrequency_*oldTime)
            );

        vector curRotationOrigin =
            initialRotationOrigin_
          + translationDirection_
           *translationAmplitude_
           *sin(2*mathematicalConstant::pi*translationFrequency_*oldTime);

        const pointField& oldPoints =
            boundaryMesh()[bodyPatchID_].localPoints();

        vector r0(1, 1, 1);
        r0 -= rotationAxis_*(rotationAxis_ & r0);
        r0 /= mag(r0);

        // http://mathworld.wolfram.com/RotationFormula.html
        vector r1 =
            r0*cos(rotAngle)
          + rotationAxis_*(rotationAxis_ & r0)*(1 - cos(rotAngle))
          + (r0 ^ rotationAxis_)*sin(rotAngle);

        tensor T = rotationTensor(r0, r1);

        vectorField rot =
            transform(T, oldPoints - curRotationOrigin)
          + curRotationOrigin
          - oldPoints;

        // Check mesh motion solver type
        bool feMotionSolver =
            this->objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );

        bool fvMotionSolver =
            this->objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    this->objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bodyPatchID_]
                );

            tetPolyPatchInterpolation tppi
            (
                refCast<const faceTetPolyPatch>(motionUPatch.patch())
            );

            motionUPatch ==
                tppi.pointToPointInterpolate
                (
                    (trans + rot)/time().deltaT().value()
                );
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[bodyPatchID_]
                );

            motionUPatch ==
                (trans + rot)/time().deltaT().value();;
        }
        else
        {
            FatalErrorIn("newDynamicBodyFvMesh::update()")
                << "Problem with mesh motion solver selection"
                    << abort(FatalError);
        }
    }

    fvMesh::movePoints(motionPtr_->newPoints());

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //
