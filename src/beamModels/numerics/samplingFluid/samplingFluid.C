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

\*---------------------------------------------------------------------------*/

#include "samplingFluid.H"
#include <tuple>
#include "surfaceMesh.H"
#include "surfaceFields.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "FieldSumOp.H"
#include "upstreamSampling.H"
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
    )
    {
    tmp<volVectorField> tresult // this is the fluid velocity stored on beam mesh
    (
        new volVectorField
        (
            IOobject
            (
                "fluidVelocity",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zero", dimVelocity, vector::zero
            )
        )
    );
    tmp<volScalarField> tmarkerResults // this is the cellMarker stored on fluid mesh
    (
        new volScalarField
        (
            IOobject
            (
                "fluidCellMarker",
                fluidMesh.time().timeName(),
                fluidMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidMesh,
            dimensionedScalar(dimless, 0)
        )
    );
    tmp<volVectorField> tmarkerVelocity // this is the beam velocity stored on fluid mesh
    (
        new volVectorField
        (
            IOobject
            (
                "markerVelocity",
                fluidMesh.time().timeName(),
                fluidMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidMesh,
            dimensionedVector
            (
                "zero", dimVelocity, vector::zero
            )
        )
    );
    volVectorField& result = tresult.ref();
    volScalarField& markerResults = tmarkerResults.ref();
    volVectorField& markerVelocity = tmarkerVelocity.ref();

    vectorField beamCoords = beamCellCenterCoord;
    vectorField beamTangents =  dRdScell;
    const volVectorField& fluidVelocity = fluidMesh.lookupObject<volVectorField>("U");

    //- this is the beam cell velocities stored on the beam mesh
    // const volVectorField& beamU = mesh.lookupObject<volVectorField>("U");

    vectorField beamUVec(beamCoords.size(), vector::zero);

    if (mesh.nCells() != 0)
    {
        forAll(beamUVec, i)
        {
            beamUVec[i] = mesh.lookupObject<volVectorField>("U")[i];   // copy raw vectors
        }
    }
    Pstream::broadcast(beamUVec);
    Pstream::broadcast(beamCoords);
    Pstream::broadcast(beamTangents);


    seedCellIDs.setSize(beamCoords.size(), -1);
    // resultVect is the fluid velocity on the beam mesh
    vectorField resultVect(beamCoords.size(), vector::zero);

    labelList fluidCellIDs(beamCoords.size(), -1);
    labelList sampleCellIDs(beamCoords.size(), -1);
    // Beam cells outside the fluid domain or below the ground, which get no
    // fluid velocity
    boolList noFluid(beamCoords.size(), false);
    boundBox bb = fluidMesh.bounds();
    bb.inflate(0.01);
    forAll(beamCoords, beamCellI)
    {
        if (!bb.contains(beamCoords[beamCellI]))
        {
            // if its not in the bounding box, return zero for fluid's U
            noFluid[beamCellI] = true;
            resultVect[beamCellI] = vector::zero;
            seedCellIDs[beamCellI] = -1;
            // if its not within the bounding box, set the fluidCellID to -1
            fluidCellIDs[beamCellI] = -1;
            continue;
        }
        // This part only should activate when we have a mooring case!
        if
        (
            groundContactActive
         && groundPlane.signedDistance(beamCoords[beamCellI]) <= 0
        )
        {
            noFluid[beamCellI] = true;
            seedCellIDs[beamCellI] = -1;
            resultVect[beamCellI] = vector::zero;
            continue;
        }
        // Seed the search from the fluid cell found for the previous beam
        // cell; the first beam cell has no previous one, so it falls through
        // to the full search below
        if (seedCellIDs[beamCellI] == -1 && beamCellI > 0)
        {
            seedCellIDs[beamCellI] = seedCellIDs[beamCellI - 1];
        }
        label fluidCellID = -1;
        if (seedCellIDs[beamCellI] != -1)
        {
            fluidCellID = searchEngine.findCell(beamCoords[beamCellI], seedCellIDs[beamCellI]);
        }

        if (fluidCellID == -1)
        {
            fluidCellID = searchEngine.findCell(beamCoords[beamCellI], -1, true);
        }
        // commented out since we are not getting the velocity here!
        if (fluidCellID != -1)
        {
            fluidCellIDs[beamCellI] = fluidCellID;

            markerVelocity[fluidCellID] = beamUVec[beamCellI];
        }
        // else
        // {
        // }

        // fluidcellID may or may not be -1
        // if -1, we could keep the old seed ...
        seedCellIDs[beamCellI] = fluidCellID;
    }
    // @todo check that each beam cell is found by only one procs.

    pointField Ps;
    vectorField upstreamDir;

    const scalar sampleDist = samplingRadius*samplingReferenceLength;

    vectorField vfAtBeam(beamCoords.size(), vector::zero);
    forAll(beamCoords, i)
    {
        if (fluidCellIDs[i] >= 0)
        {
            // local rank contributes the velocity from its own cell
            vector U_fluid = fluidVelocity[fluidCellIDs[i]];
            // We should check the beam UVec!
            vector U_beam = beamUVec[i];
            // get the relative velocity, and use it for getting the correct sampling directions
            vfAtBeam[i] = U_fluid - U_beam; // check here @todo
        }
    }
    // Reduce so that *all ranks* now have the SAME vfAtBeam
    // this assumes each beam cell maps to only one procs.
    reduce(vfAtBeam, FieldSumOp<vector>());


    //---- parallel code up until here
    if (Pstream::master())
    {
        computeUpstreamSamplingPoints
        (
            beamCoords,
            vfAtBeam,
            beamTangents,
            sampleDist,
            samplingPlaneNormal,
            globalUpstreamDir,
            Ps,
            upstreamDir
        );
    }

    // Broadcast upstream sampling points and directions so all ranks use the same Ps / upstreamDir
    Pstream::broadcast(Ps);
    Pstream::broadcast(upstreamDir);


    forAll(Ps, pI)
    {
        // A cell outside the fluid domain or below the ground keeps its zero
        // fluid velocity, rather than the velocity at its upstream sampling
        // point
        if (noFluid[pI])
        {
            continue;
        }

        label sampleCellID = -1;
        if (!fluidMesh.bounds().contains(Ps[pI]))
        {
            resultVect[pI] = vector::zero;
            continue;
        }
        // try to also add searching with seed
        sampleCellID = searchEngine.findCell(Ps[pI],-1,true);
        if (sampleCellID == -1)
        {
            const vector dirHat = upstreamDir[pI]/(mag(upstreamDir[pI]) + VSMALL);
            const point pPert = Ps[pI] + 1e-9*dirHat;   // small directional nudge
            sampleCellID = searchEngine.findCell(pPert, -1, true);
        }
        if (sampleCellID == -1)
        {
            resultVect[pI] = vector::zero;
        }
        else
        {
            resultVect[pI] = fluidVelocity[sampleCellID];
            markerResults[sampleCellID] = 1.0;
        }
    }


    reduce(resultVect, FieldSumOp<vector>());

    // resultVect spans the whole beam on every rank, but the result field
    // must keep the local beam mesh size: all cells on the master, none
    // elsewhere (see checkBeamOnMasterProcessor)
    result.primitiveFieldRef() = SubField<vector>(resultVect, mesh.nCells());

    return std::make_tuple(std::move(tresult), std::move(tmarkerResults), std::move(tmarkerVelocity), std::move(fluidCellIDs));
    }
} // End namespace Foam


// ************************************************************************* //
