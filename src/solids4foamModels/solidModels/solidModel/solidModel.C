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

#include "solidModel.H"
#include "volFields.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidModel, 0);
    defineRunTimeSelectionTable(solidModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidModel::calcGlobalFaceZones() const
{
    // Find global face zones
    if (globalFaceZonesPtr_)
    {
        FatalErrorIn
        (
            "void solidModel::calcGlobalFaceZones() const"
        )
            << "Global face zones already found"
            << abort(FatalError);
    }

    if (Pstream::parRun())
    {
        SLList<label> globalFaceZonesSet;

        // Previous method
        // const faceZoneMesh& faceZones = mesh().faceZones();
        // forAll(faceZones, zoneI)
        // {
        //     const faceZone& curFaceZone = faceZones[zoneI];
        //     forAll(curFaceZone, faceI)
        //     {
        //         // If unused face exist
        //         if (curFaceZone[faceI] >= mesh().nFaces())
        //         {
        //             globalFaceZonesSet.insert(zoneI);
        //             break;
        //         }
        //     }
        // }

        // New method: directly lookup globalFaceZones from decomposeParDict


        // For FSI cases, we need to look in a different location for the dict

        word decompDictName = "system/decomposeParDict";

        if
        (
            isDir
            (
                mesh().time().rootPath()/mesh().time().caseName()
                /"../system/solid"
            )
        )
        {
            decompDictName = "../system/solid/decomposeParDict";
        }

        Info<< "Reading decomposeParDict " << decompDictName << endl;

        IOdictionary decompDict
        (
            IOobject
            (
                decompDictName,
                mesh().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (decompDict.found("globalFaceZones"))
        {
            wordList globalFaceZoneNames(decompDict.lookup("globalFaceZones"));

            const faceZoneMesh& faceZones = mesh().faceZones();

            forAll(globalFaceZoneNames, nameI)
            {
                const label zoneID =
                    faceZones.findZoneID(globalFaceZoneNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn(type() + "::findGlobalFaceZones")
                        << "Cannot find globalFaceZone:"
                        << " " << globalFaceZoneNames[nameI]
                        << abort(FatalError);
                }

                globalFaceZonesSet.insert(zoneID);
            }

            globalFaceZonesPtr_ = new labelList(globalFaceZonesSet);
        }
        else
        {
            globalFaceZonesPtr_ = new labelList(0);
        }
    }
    else
    {
        globalFaceZonesPtr_ = new labelList(0);
    }
}


void Foam::solidModel::calcGlobalToLocalFaceZonePointMap() const
{
    // Find global face zones
    if (globalToLocalFaceZonePointMapPtr_)
    {
        FatalErrorIn
        (
            "void solidModel::calcGlobalToLocalFaceZonePointMap() const"
        )   << "Global to local face zones point map already exists"
            << abort(FatalError);
    }

    globalToLocalFaceZonePointMapPtr_ =
        new labelListList(globalFaceZones().size());

    labelListList& globalToLocalFaceZonePointMap =
        *globalToLocalFaceZonePointMapPtr_;

    const labelList& globalFaceZones = this->globalFaceZones();

    forAll(globalFaceZones, zoneI)
    {
        const label curZoneID = globalFaceZones[zoneI];

        Info<< "Creating faceMap for globalFaceZones "
            << mesh().faceZones()[curZoneID].name()<< endl;

        labelList curMap(mesh().faceZones()[curZoneID]().nPoints(), -1);

        vectorField fzGlobalPoints =
            mesh().faceZones()[curZoneID]().localPoints();

        // Set all slave points to zero because only the master order is used
        if(!Pstream::master())
        {
            fzGlobalPoints = vector::zero;
        }

        // Pass points to all procs
        reduce(fzGlobalPoints, sumOp<vectorField>());

        // Now every proc has the master's list of FZ points
        // every proc must now find the mapping from their local FZ points to
        // the global FZ points

        const vectorField& fzLocalPoints =
            mesh().faceZones()[curZoneID]().localPoints();

        const edgeList& fzLocalEdges =
            mesh().faceZones()[curZoneID]().edges();

        const labelListList& fzPointEdges =
            mesh().faceZones()[curZoneID]().pointEdges();

        scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

        forAll(minEdgeLength, pI)
        {
            const labelList& curPointEdges = fzPointEdges[pI];

            forAll(curPointEdges, eI)
            {
                const scalar Le =
                    fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);

                if (Le < minEdgeLength[pI])
                {
                    minEdgeLength[pI] = Le;
                }
            }
        }

        forAll(fzGlobalPoints, globalPointI)
        {
            boolList visited(fzLocalPoints.size(), false);

            forAll(fzLocalPoints, procPointI)
            {
                if (!visited[procPointI])
                {
                    visited[procPointI] = true;

                    label nextPoint = procPointI;

                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[nextPoint]
                          - fzGlobalPoints[globalPointI]
                        );

                    if (curDist < 1e-4*minEdgeLength[nextPoint])
                    {
                        curMap[globalPointI] = nextPoint;
                        break;
                    }

                    label found = false;

                    while (nextPoint != -1)
                    {
                        const labelList& nextPointEdges =
                            fzPointEdges[nextPoint];

                        scalar minDist = GREAT;
                        label index = -1;
                        forAll(nextPointEdges, edgeI)
                        {
                            label curNgbPoint =
                                fzLocalEdges[nextPointEdges[edgeI]]
                               .otherVertex(nextPoint);

                            if (!visited[curNgbPoint])
                            {
                                visited[curNgbPoint] = true;

                                scalar curDist =
                                    mag
                                    (
                                        fzLocalPoints[curNgbPoint]
                                      - fzGlobalPoints[globalPointI]
                                    );

                                if (curDist < 1e-4*minEdgeLength[curNgbPoint])
                                {
                                    curMap[globalPointI] = curNgbPoint;
                                    found = true;
                                    break;
                                }
                                else if (curDist < minDist)
                                {
                                    minDist = curDist;
                                    index = curNgbPoint;
                                }
                            }
                        }

                        nextPoint = index;
                    }

                    if (found)
                    {
                        break;
                    }
                }
            }
        }

        forAll(curMap, globalPointI)
        {
            if (curMap[globalPointI] == -1)
            {
                FatalErrorIn
                (
                    "solidModel::calcGlobalToLocalFaceZonePointMap()"
                )   << "local to global face zone point map is not correct"
                    << abort(FatalError);
            }
        }

        globalToLocalFaceZonePointMap[zoneI] = curMap;
    }
}


void Foam::solidModel::updateGlobalFaceZoneNewPoints
(
    const pointField& pointDDI,
    pointField& newPoints
)
{
    const labelList& globalFaceZones = this->globalFaceZones();
    const labelListList& globalToLocalFaceZonePointMap =
        this->globalToLocalFaceZonePointMap();

    forAll(globalFaceZones, zoneI)
    {
        const label curZoneID = globalFaceZones[zoneI];

        const labelList& curMap = globalToLocalFaceZonePointMap[zoneI];

        const labelList& curZoneMeshPoints =
            mesh().faceZones()[curZoneID]().meshPoints();

        vectorField curGlobalZonePointDispl
        (
            curZoneMeshPoints.size(),
            vector::zero
        );

        // Inter-proc points are shared by multiple procs
        // pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

        forAll(curGlobalZonePointDispl, globalPointI)
        {
            label localPoint = curMap[globalPointI];

            if(curZoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = curZoneMeshPoints[localPoint];

                curGlobalZonePointDispl[globalPointI] = pointDDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(curGlobalZonePointDispl, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            // Now average the displacement between all procs
            curGlobalZonePointDispl /= pointNumProcs;
        }

        // The curZonePointsDisplGlobal now contains the correct face zone
        // displacement in a global master processor order, now convert them
        // back into the local proc order

        vectorField curZonePointDispl(curZoneMeshPoints.size(), vector::zero);

        forAll(curGlobalZonePointDispl, globalPointI)
        {
            label localPoint = curMap[globalPointI];

            curZonePointDispl[localPoint] =
                curGlobalZonePointDispl[globalPointI];
        }

        forAll(curZonePointDispl, pointI)
        {
            // Unused points
            if (curZoneMeshPoints[pointI] >= mesh().nPoints())
            {
                newPoints[curZoneMeshPoints[pointI]] +=
                    curZonePointDispl[pointI];
            }
        }
    }
}


void Foam::solidModel::makeMechanicalModel() const
{
    if (!mechanicalPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeMechanicalModel() const")
            << "pointer alrady set!" << abort(FatalError);
    }

    mechanicalPtr_.set
    (
        new mechanicalModel(mesh(), nonLinGeom())
    );
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

Foam::mechanicalModel& Foam::solidModel::mechanical()
{
    if (mechanicalPtr_.empty())
    {
        makeMechanicalModel();
    }

    return mechanicalPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    dynamicFvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "solidProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    solidProperties_(subDict(type + "Coeffs")),
    mechanicalPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{
    deleteDemandDrivenData(globalFaceZonesPtr_);
    deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::mechanicalModel& Foam::solidModel::mechanical() const
{
    if (mechanicalPtr_.empty())
    {
        makeMechanicalModel();
    }

    return mechanicalPtr_();
}


const Foam::labelList& Foam::solidModel::globalFaceZones() const
{
    if (!globalFaceZonesPtr_)
    {
        calcGlobalFaceZones();
    }

    return *globalFaceZonesPtr_;
}


const Foam::labelListList&
Foam::solidModel::globalToLocalFaceZonePointMap() const
{
    if (!globalToLocalFaceZonePointMapPtr_)
    {
        calcGlobalToLocalFaceZonePointMap();
    }

    return *globalToLocalFaceZonePointMapPtr_;
}


void Foam::solidModel::setTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setTraction(patchID, patchTraction);
}


void Foam::solidModel::setPressure
(
    const label patchID,
    const label zoneID,
    const scalarField& faceZonePressure
)
{
    scalarField patchPressure(mesh().boundary()[patchID].size(), 0.0);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchPressure, i)
    {
        patchPressure[i] =
            faceZonePressure
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setPressure(patchID, patchPressure);
}


void Foam::solidModel::writeFields(const Time& runTime)
{
    runTime.write();
}


Foam::scalar Foam::solidModel::newDeltaT()
{
    return mechanical().newDeltaT();
}


void Foam::solidModel::moveMesh
(
    const pointField& oldPoints,
    const volVectorField& DD,
    pointVectorField& pointDD
)
{
    Info<< "Moving the mesh to the deformed configuration" << nl << endl;

    //- Move mesh by interpolating displacement field to vertices
    // To be checked: sync boundary and global points across procs to make sure
    // numiercal error does not build up and when end up with the error
    // "face area does not match neighbour..."
    // We could sync points as a pointVectorField just as we sync pointDD

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD, pointDD);

    // Ensure continuous displacement across processor boundary
    // Something strange is happening here
    pointDD.correctBoundaryConditions();

    vectorField& pointDDI = pointDD.internalField();

    vectorField newPoints = oldPoints;

    // Correct symmetryPlane points

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());

            const vector i(1, 0, 0);
            const vector j(0, 1, 0);
            const vector k(0, 0, 1);

            if (mag(avgN & i) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].x() = 0;
                }
            }
            else if (mag(avgN & j) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].y() = 0;
                }
            }
            else if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
        else if (isA<emptyPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());
            const vector k(0, 0, 1);

            if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
    }

    // Note: allPoints will have more points than pointDD if there are
    // globalFaceZones
    forAll (pointDDI, pointI)
    {
        newPoints[pointI] += pointDDI[pointI];
    }

    // Move unused globalFaceZone points
    updateGlobalFaceZoneNewPoints(pointDDI, newPoints);

    twoDPointCorrector twoDCorrector(mesh());
    twoDCorrector.correctPoints(newPoints);
    twoDCorrector.correctPoints(pointDDI);
    mesh().movePoints(newPoints);
    mesh().V00();
    mesh().moving(false);
    mesh().changing(false);
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

    // Tell the mechanical model to move the subMeshes, if they exist
    mechanical().moveSubMeshes();
}


bool Foam::solidModel::read()
{
    if (regIOobject::read())
    {
        solidProperties_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
