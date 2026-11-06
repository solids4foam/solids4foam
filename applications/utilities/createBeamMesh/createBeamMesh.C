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
    Create beam meshes for the beamFoam solver or solid model.

Authors
    Zeljko Tukovic, FSB
    Seevani Bali, UCD
    Philip Cardiff, UCD, All rights reserved.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "faceList.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "IOdictionary.H"
#include "cellZone.H"
#include "crossSectionModel.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    // argList::validOptions.insert("decomposePatches", "");
#   include "setRootCase.H"
#   include "createTime.H"

    word regionName;
    fileName polyMeshDir;

    if (args.found("region"))
    {
        // constant/<region>/polyMesh/blockMeshDict
        regionName  = args.get("region");
        polyMeshDir = regionName/polyMesh::meshSubDir;

        Info<< nl << "Generating mesh for region " << regionName << endl;
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        regionName  = polyMesh::defaultRegion;
        polyMeshDir = polyMesh::meshSubDir;
    }

    // const bool decomposePatches = args.found("decomposePatches");

    IOdictionary beamProperties
    (
        IOobject
        (
            "beamProperties",
            bool(regionName == polyMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/regionName),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    int nBeams = -1;
    scalarField R;
    scalarField L;
    labelList nZ;
    labelList nTheta;

    // Read cross-sections of beams
    PtrList<crossSectionModel> crossSections;
    if (beamProperties.found("beams"))
    {
        const PtrList<entry> entries(beamProperties.lookup("beams"));

        nBeams = entries.size();

        crossSections.setSize(nBeams);

        forAll(entries, beamI)
        {
            crossSections.set
            (
                beamI,
                crossSectionModel::New
                (
                 word(entries[beamI].dict().lookup("crossSectionModel")),
                 entries[beamI].dict()
                )
            );
        }

        R.setSize(nBeams);
        L.setSize(nBeams);
        nZ.setSize(nBeams);
        nTheta.setSize(nBeams);
        forAll(crossSections, beamI)
        {
            R[beamI] = crossSections[beamI].R();

            L[beamI] = readScalar(entries[beamI].dict().lookup("length"));

            nZ[beamI] = readInt(entries[beamI].dict().lookup("nSegments"));

            nTheta[beamI] = crossSections[beamI].nCircumferentialPoints();
        }
    }

    // if (beamProperties.found("nBeams"))
    // {
    //     nBeams = readInt(beamProperties.lookup("nBeams"));
    //     R.setSize
    //     (
    //         nBeams,
    //         dimensionedScalar(beamProperties.lookup("R")).value()
    //     );
    //     L.setSize
    //     (
    //         nBeams,
    //         dimensionedScalar(beamProperties.lookup("L")).value()
    //     );
    //     nZ.setSize
    //     (
    //         nBeams,
    //         readInt(beamProperties.lookup("nSegments"))
    //     );
    // }
    // else
    // {
    //     R = scalarField(beamProperties.lookup("R"));
    //     L = scalarField(beamProperties.lookup("L"));
    //     nZ = labelList(beamProperties.lookup("nSegments"));
    //     if (R.size() != L.size())
    //     {
    //         FatalErrorIn(args.executable())
    //             << "Sizes of R and L fields are not equal."
    //             << abort(FatalError);
    //     }

    //     if (R.size() != nZ.size())
    //     {
    //         FatalErrorIn(args.executable())
    //             << "Sizes of R and nZ fields are not equal."
    //             << abort(FatalError);
    //     }

    //     nBeams = R.size();
    // }

    if (nBeams < 1)
    {
        FatalErrorIn(args.executable())
            << "Number of beams is not defined."
            << abort(FatalError);
    }

    // Scalar R = dimensionedScalar(beamProperties.lookup("R")).value();
    // scalar L = dimensionedScalar(beamProperties.lookup("L")).value();
    // scalar R = 0.1;
    // scalar L = 20;

    // label nZ(readInt(beamProperties.lookup("nSegments")));
    // label nZ = 50;
    // label nTheta = 32;

    // const int nBeams
    // (
    //     beamProperties.lookupOrDefault<int>("nBeams", 1)
    // );

    labelList nOneBeamCells = nZ;

    labelList cellStart(nBeams, -1);
    cellStart[0] = 0;
    for (label i = 1; i < nBeams; i++)
    {
        cellStart[i] = cellStart[i - 1] + nOneBeamCells[i - 1];
    }

    labelList nOneBeamPoints = nTheta*(nZ + 1);
    label nPoints = sum(nOneBeamPoints); //*nBeams;

    labelList pointStart(nBeams, -1);
    pointStart[0] = 0;
    for (label i = 1; i < nBeams; i++)
    {
        pointStart[i] = pointStart[i - 1] + nOneBeamPoints[i - 1];
    }

    labelList nOneBeamFaces = ((nZ + 1) + nZ*nTheta);
    label nFaces = sum(nOneBeamFaces); //*nBeams;

    labelList nOneBeamInternalFaces = (nZ - 1);
    label nInternalFaces = sum(nOneBeamInternalFaces); //*nBeams;

    // Points
    vectorField points(nPoints, vector::zero);

    label pI = 0;
    for (label k = 0; k < nBeams; k++)
    {
        // scalar dTheta = -2*M_PI/nTheta[k];

        // vectorField circumferencePoints(nTheta[k], vector::zero);
        // forAll(circumferencePoints, pointI)
        // {
        //     circumferencePoints[pointI].x() = R[k]*::cos(pointI*dTheta);
        //     circumferencePoints[pointI].y() = R[k]*::sin(pointI*dTheta);
        //     circumferencePoints[pointI].z() = 0;
        // }

        const vectorField circumferencePoints
        (
            crossSections[k].circumferentialPoints()
        );

        scalar dZ = L[k]/nZ[k];
        for (label i = 0; i <= nZ[k]; i++)
        {
            for (label j = 0; j < nTheta[k]; j++)
            {
                vector curPoint = circumferencePoints[j];
                curPoint.x() = i*dZ;

                points[pI++] = curPoint;
            }
        }
    }

    // Faces
    faceList faces(nFaces, face());

    //-Internal faces

    label fI = 0;
    for (label k = 0; k < nBeams; k++)
    {
        for (label i = 0; i <= nZ[k]; i++)
        {
            labelList curPointLabels(nTheta[k], -1);
            for (label j = 0; j < nTheta[k]; j++)
            {
                curPointLabels[j] = j + i*nTheta[k] + pointStart[k];
                // curPointLabels[j] = j + i*nTheta[k] + k*nOneBeamPoints;
            }

            if (i > 0 && i < nZ[k])
            {
                face curFace(curPointLabels);
                curFace = curFace.reverseFace();
                faces[fI++] = curFace;
            }
        }
    }

    //- Inlet faces

    for (label k = 0; k < nBeams; k++)
    {
        labelList curPointLabels(nTheta[k], -1);
        for (label j = 0; j < nTheta[k]; j++)
        {
            curPointLabels[j] = j + pointStart[k]; //k*nOneBeamPoints;
        }

        face curFace(curPointLabels);
        // curFace = curFace.reverseFace();
        faces[fI++] = curFace;
    }

    //- Outlet faces

    for (label k = 0; k < nBeams; k++)
    {
        labelList curPointLabels = labelList(nTheta[k], -1);
        for (label j = 0; j < nTheta[k]; j++)
        {
            curPointLabels[j] = j + nZ[k]*nTheta[k] + pointStart[k]; //k*nOneBeamPoints;
        }

        face curFace = face(curPointLabels);
        curFace = curFace.reverseFace();
        faces[fI++] = curFace;
    }

    //- beam faces

    for (label k = 0; k < nBeams; k++)
    {
        for (label i = 0; i < nZ[k]; i++)
        {
            for (label j = 1; j < nTheta[k]; j++)
            {
                labelList curPointLabels(4, -1);

                curPointLabels[0] = (j - 1) + i*nTheta[k];
                curPointLabels[1] = (j - 1) + (i + 1)*nTheta[k];
                curPointLabels[2] = j + (i + 1)*nTheta[k];
                curPointLabels[3] = j + i*nTheta[k];

                forAll(curPointLabels, pI)
                {
                    curPointLabels[pI] += pointStart[k]; //k*nOneBeamPoints;
                }

                face curFace(curPointLabels);
                // curFace = curFace.reverseFace();
                faces[fI++] = curFace;
            }

            labelList curPointLabels(4, -1);

            curPointLabels[0] = (nTheta[k] - 1) + i*nTheta[k];
            curPointLabels[1] = (nTheta[k] - 1) + (i + 1)*nTheta[k];
            curPointLabels[2] = 0 + (i + 1)*nTheta[k];
            curPointLabels[3] = 0 + i*nTheta[k];

            forAll(curPointLabels, pI)
            {
                curPointLabels[pI] += pointStart[k]; //k*nOneBeamPoints;
            }

            face curFace(curPointLabels);
            // curFace = curFace.reverseFace();
            faces[fI++] = curFace;
        }
    }

    // Owners
    labelList owners(nFaces, -1);

    fI = 0;
    for (label k = 0; k < nBeams; k++)
    {
        for (label i = 0; i <= nZ[k]; i++)
        {
            if (i > 0 && i < nZ[k])
            {
                owners[fI++] = (i - 1) + cellStart[k]; //k*nOneBeamCells;
            }
        }
    }

    // Inlet faces
    for (label k = 0; k < nBeams; k++)
    {
        owners[fI++] = cellStart[k]; //k*nOneBeamCells;
    }

    // Outlet faces
    for (label k = 0; k < nBeams; k++)
    {
        owners[fI++] = (nZ[k] - 1) + cellStart[k]; //k*nOneBeamCells;
    }

    // Beam faces
    for (label k = 0; k < nBeams; k++)
    {
        for (label i = 0; i < nZ[k]; i++)
        {
            for (label j = 1; j < nTheta[k]; j++)
            {
                owners[fI++] = i + cellStart[k]; //k*nOneBeamCells;
            }

            owners[fI++] = i + cellStart[k]; //k*nOneBeamCells;
        }
    }

    // Neighbours
    labelList neighbours(nInternalFaces, -1);

    fI = 0;
    for (label k = 0; k < nBeams; k++)
    {
        for (label i = 1; i <= (nZ[k]-1); i++)
        {
            neighbours[fI++] = i + cellStart[k]; //k*nOneBeamCells;
        }
    }

    // Create mesh
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points),
        std::move(faces),
        std::move(owners),
        std::move(neighbours)
    );

    // if ((nBeams > 1) && decomposePatches)
    if (nBeams > 1)
    {
        // The boundary faces of multiple beams form multiple
        // patches but the circumeference faces of internal
        // region of all beams are combined into 'beam'
        // patch that is empty.
        List<polyPatch*> patches(nBeams + nBeams + 1);

        label pI = 0;
        label start = nInternalFaces;

        // Inlet patches

        for (label k = 0; k < nBeams; k++)
        {
            OStringStream PatchName;
            PatchName() << "left_" << k;

            patches[pI] =
                polyPatch::New
                (
                    polyPatch::typeName,
                    word(PatchName.str()),
                    1,
                    start,
                    pI,
                    mesh.boundaryMesh()
                ).ptr();
            pI++;
            start += 1;
        }

        // Outlet patches
        for (label k = 0; k < nBeams; k++)
        {
            OStringStream PatchName;
            PatchName() << "right_" << k;

            patches[pI] =
                polyPatch::New
                (
                    polyPatch::typeName,
                    word(PatchName.str()),
                    1,
                    start,
                    pI,
                    mesh.boundaryMesh()
                ).ptr();
            pI++;
            start += 1;
        }

        // Beam patch
        patches[pI] =
            polyPatch::New
            (
                emptyPolyPatch::typeName,
                "beam",
                sum(nTheta*nZ),
                // (nTheta[k]*nZ)*nBeams,
                start,
                pI,
                mesh.boundaryMesh()
            ).ptr();

        mesh.addPatches(patches);
    }
    else
    {
        List<polyPatch*> patches(3);

        // Inlet patch
        patches[0] =
            polyPatch::New
            (
                polyPatch::typeName,
                "left",
                nBeams,
                nInternalFaces,
                0,
                mesh.boundaryMesh()
            ).ptr();

        // Outlet patch
        patches[1] =
            polyPatch::New
            (
                polyPatch::typeName,
                "right",
                nBeams,
                nInternalFaces + nBeams,
                1,
                mesh.boundaryMesh()
            ).ptr();

        // Beam patch
        patches[2] =
            polyPatch::New
            (
                emptyPolyPatch::typeName,
                "beam",
                sum(nTheta*nZ),
                // (nTheta[k]*nZ)*nBeams,
                nInternalFaces + nBeams + nBeams,
                2,
                mesh.boundaryMesh()
            ).ptr();

        mesh.addPatches(patches);
    }

    //mesh.removeZones();

    List<pointZone*> pz(nBeams);
    List<faceZone*> fz(0);
    List<cellZone*> cz(nBeams);

    for (label k = 0; k < nBeams; k++)
    {
        labelList curBeamCells(nOneBeamCells[k], -1);
        labelList curBeamPoints(nOneBeamPoints[k], -1);
        forAll(curBeamCells, cellI)
        {
            curBeamCells[cellI] = cellI + cellStart[k]; //k*nOneBeamCells;
        }
        forAll(curBeamPoints, pointI)
        {
            curBeamPoints[pointI] = pointI + pointStart[k]; //k*nOneBeamCells;
        }

        OStringStream RegionName;
        RegionName() << "beam_" << k;

        pz[k] =
            new pointZone
            (
                word(RegionName.str()),
                curBeamPoints,
                k,
                mesh.pointZones()
            );

        cz[k] =
            new cellZone
            (
                word(RegionName.str()),
                curBeamCells,
                k,
                mesh.cellZones()
            );
    }

    mesh.addZones(pz, fz, cz);

    mesh.write();

    #include "printMeshSummary.H"

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
