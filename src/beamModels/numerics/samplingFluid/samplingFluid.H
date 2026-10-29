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
    Helper function for sampling fluid fields in two way coupled beam fluid interaction problems.
Author
    Amirhossein Taran, UCD
    Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#ifndef samplingFluid_H
#define samplingFluid_H

#include "fvMesh.H"
#include "volFields.H"
#include "meshSearch.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    std::tuple
        <
            tmp<volVectorField>,
            tmp<volScalarField>,
            tmp<volVectorField>,
            labelList
        >
        getFluidVelocity
        (
            const fvMesh& fluidMesh,
            const fvMesh& mesh,
            const vectorField& beamCellCenterCoord,
            labelList& seedCellIDs,
            const plane& groundPlane,
            const bool groundContactActive,
            const vectorField& dRdScell,
            const scalar samplingReferenceLength,
            const scalar samplingRadius,
            const vector& samplingPlaneNormal,
            const vector& globalUpstreamDir,
            const meshSearch& searchEngine
        );
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
