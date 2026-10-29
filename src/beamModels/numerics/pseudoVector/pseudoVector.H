/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    This class contains the mathematical formulae to:
    1. Construct a 3 x 3 rotation matrix from a 3-D rotation vector using the Rodrigues' formula (member functions 'rotationMatrix')
    2. Construct a (pseudo-)vector from a specified rotation matrix using the Spurrier algorithm, wherein
       quaternions are calculated from the rotation matix and then using the quaternions, the corresponding
       3-D (pseudo)-rotation vector is calculated (member functions 'pseudoVector')
    3. Interpolate rotation matrices between cell-centre and face-centre values of a
       beam cell-segment (member functions 'interpolateRotationMatrix')
    4. Interpolate rotational strains between cell-centre and face-centre values of a
       beam cell-segment (member function 'meanLineCurvature')
    5. Calculation of spatial derivative of the rotation matrix (member function 'rotationAngleDerivative')

    Currently used by the solver: 'pseudoVector', 'rotationMatrix' and the
    face-to-cell 'interpolateRotationMatrix(bm, Rf, R)', which calculates the
    cell-centre reference rotation matrix refLambda from refLambdaf.

    Currently unused (kept for reference; compiled so that they stay valid):
    - 'interpolateRotationMatrix(R)' and 'meanLineCurvature(R)': called only
      by the former 'objectiveInterpolation' option, see the notes in
      pseudoVector.C
    - 'interpolateRotationMatrix(R, Rf, snGradRotAngle)' and
      'rotationAngleDerivative(R)': no caller in the history of beamFoam

Author
    Zeljko Tukovic, FSB Zagreb
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pseudoVector_H
#define pseudoVector_H

#include "volFields.H"
#include "pointFields.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<volVectorField> pseudoVector(const volTensorField& rotationMatrix);

tmp<vectorField> pseudoVector(const tensorField& rotationMatrix);

vector pseudoVector(const tensor& R);

tmp<tensorField> rotationMatrix(const vectorField& rotatioAngle);
tmp<volTensorField> rotationMatrix(const volVectorField& rotatioAngle);
tmp<surfaceTensorField> rotationMatrix(const surfaceVectorField& rotatioAngle);

tensor rotationMatrix(const vector& rotationAngle);

//- Cell-to-face objective interpolation of rotation matrices
//  NOTE: currently unused (former objectiveInterpolation option)
tmp<surfaceTensorField> interpolateRotationMatrix
(
    const volTensorField& R
);

//- Face curvature from neighbouring cell rotation matrices
//  NOTE: currently unused (former objectiveInterpolation option)
tmp<surfaceVectorField> meanLineCurvature(const volTensorField& R);

//- Cell-to-face interpolation of rotation matrices and face rotation
//  angle gradient
//  NOTE: currently unused
void interpolateRotationMatrix
(
    const volTensorField& R,
    surfaceTensorField& Rf,
    surfaceVectorField& snGradRotAngle
);

//- Face-to-cell interpolation of rotation matrices
void interpolateRotationMatrix
(
    const beamModel& bm,
    const surfaceTensorField& Rf,
    volTensorField& R
);

//- Face rotation angle gradient from neighbouring cell rotation matrices
//  NOTE: currently unused
tmp<surfaceVectorField> rotationAngleDerivative
(
    const volTensorField& rotatioMatrix
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
