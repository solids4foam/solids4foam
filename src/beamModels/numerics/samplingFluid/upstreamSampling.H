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
    Helper function for Actuator line model's upstream sampling

Author
    Amirhossein Taran, UCD

\*---------------------------------------------------------------------------*/
#ifndef upstreamSampling_H
#define upstreamSampling_H

#include "fvMesh.H"
#include "vector.H"
#include "pointField.H"
#include "volFields.H"
#include "IOmanip.H"

namespace Foam
{
//     inline vector upstreamDirPerpToTangent(const vector& vf, const vector& beamTangent)
// {
//     vector t = beamTangent/(mag(beamTangent) + VSMALL);

//     // 1. Force the direction to be based on the true global freestream (-X)
//     vector globalUpstream(-1.0, 0.0, 0.0);

//     // 2. Make it perpendicular to the local beam tangent
//     vector samplingDir = globalUpstream - (globalUpstream & t)*t;

//     // 3. Prevent Z-escape for the 3D thin mesh
//     samplingDir.z() = 0.0;

//     scalar m = mag(samplingDir);
//     // 4. Fallback in case the beam ever bends to be completely horizontal
//     if (m < VSMALL)
//     {
//         samplingDir = vector(0.0, 1.0, 0.0);
//         m = 1.0;
//     }

//     return samplingDir / m;
// }
inline vector upstreamDirPerpToTangent
(
    const vector& vf,
    const vector& beamTangent,
    const vector& globalUpstreamDir
)
{
    const vector upHat = globalUpstreamDir/(mag(globalUpstreamDir) + VSMALL);
    const vector t = beamTangent/(mag(beamTangent) + VSMALL);
    const vector vPerp = vf - (vf & t)*t;

    vector samplingDir = -vPerp;

    // If samplingDir is not in the upstream half-space, project globalUpstreamDir
    // perpendicular to the beam tangent and use that instead.
    if ((samplingDir & upHat) <= 0)
    {
        samplingDir = upHat - (upHat & t)*t;
    }

    scalar m = mag(samplingDir);
    if (m < VSMALL)
    {
        // Beam nearly parallel to globalUpstreamDir: cannot form a meaningful
        // perpendicular direction. Fall back to globalUpstreamDir directly.
        WarningInFunction
            << "Sampling direction magnitude is below VSMALL." << nl
            << "The beam tangent is nearly parallel to almUpstreamDir." << nl
            << "Using almUpstreamDir = " << globalUpstreamDir
            << " directly as the sampling direction." << nl
            << "vf = " << vf << nl
            << "beamTangent = " << beamTangent << nl
            << endl;
        return upHat;
    }
    return samplingDir/(m + VSMALL);
}

// Compute sampling points Ps (step 1 of ALM)
inline void computeUpstreamSamplingPoints
(
    const pointField& beamCoords,
    const vectorField& vfAtBeam,
    const vectorField& beamTangents,
    const scalar sampleDist,
    const vector& samplingPlaneNormal,
    const vector& globalUpstreamDir,
    pointField& Ps,
    vectorField& upstreamDir
)
{
    Ps.setSize(beamCoords.size(), vector::zero);
    upstreamDir.setSize(beamCoords.size(), vector::zero);

    const vector nHat = samplingPlaneNormal/(mag(samplingPlaneNormal) + VSMALL);

    forAll(beamCoords, i)
    {
        const vector vf = vfAtBeam[i];
        const vector dir =
            upstreamDirPerpToTangent(vf, beamTangents[i], globalUpstreamDir);
        upstreamDir[i] = dir;
        Ps[i] = beamCoords[i] + sampleDist*dir;
        Ps[i] -= ((Ps[i] - beamCoords[i]) & nHat)*nHat;
    }
}

} // End namespace Foam

#endif
