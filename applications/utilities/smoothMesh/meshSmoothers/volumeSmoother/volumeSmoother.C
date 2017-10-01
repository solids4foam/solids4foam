/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "volumeSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "twoDPointCorrector.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volumeSmoother, 0);
    addToRunTimeSelectionTable(meshSmoother, volumeSmoother, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
volumeSmoother::volumeSmoother
(
    const word& name,
    fvMesh& mesh
)
:
    meshSmoother(name, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

volumeSmoother::~volumeSmoother()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Smoothing function
scalar volumeSmoother::smooth(bool writeIters)
{
    const scalar lambda(readScalar(dict().lookup("lambda")));
    const scalar mu(readScalar(dict().lookup("mu")));
    const int nCorr(readInt(dict().lookup("nCorrectors")));

    fvMesh& mesh = this->mesh();

    boolList cellZonePoints(mesh.nPoints(), false);

    const scalar edgeAngle
    (
        dict().lookupOrDefault<scalar>
        (
            "edgeAngle",
            45
        )
    );
    const scalar edgeAngleDot =
        Foam::cos(edgeAngle*mathematicalConstant::pi/180.0);

    const word cellZoneName
    (
        dict().lookupOrDefault<word>
        (
            "cellZone",
            "All"
        )
    );

    label cellZoneID = -1;

    if (cellZoneName != "All")
    {

        cellZoneID = mesh.cellZones().findZoneID(cellZoneName);

        if (cellZoneID == -1)
        {
            FatalError
                << "cellZone " << cellZoneName << " not found!"
                << abort(FatalError);
        }

        Info<< "Smoothing points in cellZone " << cellZoneName
            << nl << endl;

        // Mark all points that are part of any cell in the cellZone
        const labelListList& cellPoints = mesh.cellPoints();
        forAll(mesh.cellZones()[cellZoneID], cellI)
        {
            const label cellID = mesh.cellZones()[cellZoneID][cellI];
            forAll(cellPoints[cellID], cpI)
            {
                cellZonePoints[cellPoints[cellID][cpI]] = true;
            }
        }
    }
    else
    {
        Info<< "Smoothing all points" << endl;
        cellZonePoints = true;
    }

    // Designate boundary points
    boolList boundaryPoints(mesh.nPoints(), false);

    forAll(mesh.boundaryMesh(), patchI)
    {
        const labelList& meshPointsOfBoundaryPatch =
            mesh.boundaryMesh()[patchI].meshPoints();

        forAll(meshPointsOfBoundaryPatch, pointI)
        {
            boundaryPoints[meshPointsOfBoundaryPatch[pointI]] = true;
        }
    }

    const labelListList& pointCells = mesh.pointCells();
    const labelListList& pointFaces = mesh.pointFaces();

    scalarField residual(mesh.nPoints(), 0);
    label counter = 0;
    bool shrink = true;
    pointField oldPoints = mesh.points();
    volVectorField oldC = mesh.C();
    vectorField oldCI = mesh.C().internalField();
    surfaceScalarField oldMagSf = mesh.magSf();
    scalarField oldV = mesh.V();
    pointField newPoints = oldPoints;
    do
    {
        oldPoints = newPoints;
        oldCI = mesh.C().internalField();
        oldV = mesh.V();

        if (!shrink)
        {
            counter++;
        }

        forAll(newPoints, pointI)
        {
            // Smooth points in cellZone
            if (cellZonePoints[pointI])
            {
                vector curNewPoint = vector::zero;

                scalar sumW = 0;

                const labelList& neiCells = pointCells[pointI];

                // for all neighbour cell centres
                forAll(neiCells, pcI)
                {
                    const label neiCellID = neiCells[pcI];

                    vector d = oldCI[neiCellID] - oldPoints[pointI];

                    // Weight based on volumes of adjacent cells
                    scalar w = oldV[neiCellID];

                    curNewPoint += w*d;

                    sumW += w;
                }

                curNewPoint /= sumW;

                curNewPoint += oldPoints[pointI];

                vector disp = curNewPoint - oldPoints[pointI];

                if (shrink)
                {
                    // Relaxed motion
                    newPoints[pointI] = oldPoints[pointI] + lambda*disp;

                    residual[pointI] = lambda*mag(disp);
                }
                else
                {
                    // inflate
                    newPoints[pointI] = oldPoints[pointI] - mu*disp;

                    residual[pointI] = mu*mag(disp);
                }
            }
        }

        // Correct patch motion to be slip beahviour
        forAll(mesh.boundaryMesh(), patchI)
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            const pointField& pointNormals =
                mesh.boundaryMesh()[patchI].pointNormals();

            // Smooth boundary points independently based on
            // neighbouring boundary face areas
            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                if (cellZonePoints[pointID])
                {
                    // Calculate boundary points based on neighbouring boundary
                    // face areas
                    vector curNewPoint = vector::zero;

                    scalar sumW = 0;

                    const labelList& neiFaces = pointFaces[pointID];

                    // for all neighbour cell centres
                    forAll(neiFaces, pfI)
                    {
                        const label neiFaceID = neiFaces[pfI];

                        if (!mesh.isInternalFace(neiFaceID))
                        {
                            const label opID =
                                mesh.boundaryMesh().whichPatch(neiFaceID);

                            if (mesh.boundaryMesh()[opID].type() != "empty")
                            {
                                const label start =
                                    mesh.boundaryMesh()[opID].start();

                                vector d =
                                    oldC.boundaryField()
                                    [
                                        opID
                                    ][neiFaceID - start]
                                  - oldPoints[pointID];

                                // Weight based on volumes of adjacent cells
                                scalar w =
                                    oldMagSf.boundaryField()
                                    [
                                        opID
                                    ][neiFaceID - start];

                                curNewPoint += w*d;

                                sumW += w;
                            }
                        }
                    }

                    vector disp = vector::zero;

                    if (sumW > SMALL)
                    {
                        curNewPoint /= sumW;

                        curNewPoint += oldPoints[pointID];

                        disp = curNewPoint - oldPoints[pointID];

                        if (shrink)
                        {
                            // Relaxed motion
                            newPoints[pointID] =
                                oldPoints[pointID] + lambda*disp;

                            residual[pointID] = lambda*mag(disp);
                        }
                        else
                        {
                            // inflate
                            newPoints[pointID] = oldPoints[pointID] - mu*disp;

                            residual[pointID] = mu*mag(disp);
                        }
                    }


                    // Keep only motion in the patch plane i.e. remove motion
                    // normal to the patch

                    disp = newPoints[pointID] - oldPoints[pointID];

                    newPoints[pointID] =
                        oldPoints[pointID]
                        + (
                            (I - sqr(pointNormals[pI])) & disp
                        );

                    residual[pointID] = mag((I - sqr(pointNormals[pI])) & disp);
                }
            }

            // Correct patch boundary edge motion

            const labelListList& edgeLoops =
                mesh.boundaryMesh()[patchI].edgeLoops();

            forAll(edgeLoops, elI)
            {
                const labelList edgeLoopI = edgeLoops[elI];

                forAll(edgeLoopI, i)
                {
                    label prevPointID = -1;
                    const label pointID = meshPoints[edgeLoopI[i]];
                    label nextPointID = -1;

                    if (i == 0)
                    {
                        prevPointID =
                            meshPoints[edgeLoopI[edgeLoopI.size() - 1]];
                        nextPointID =
                            meshPoints[edgeLoopI[1]];
                    }
                    else if (i == edgeLoopI.size() - 1)
                    {
                        prevPointID =
                            meshPoints[edgeLoopI[i - 1]];
                        nextPointID =
                            meshPoints[edgeLoopI[0]];
                    }
                    else
                    {
                        prevPointID =
                            meshPoints[edgeLoopI[i - 1]];
                        nextPointID =
                            meshPoints[edgeLoopI[i + 1]];
                    }

                    const bool independentEdges = true;
                    if (independentEdges && cellZonePoints[pointID])
                    {
                        vector d0 = oldPoints[prevPointID] - oldPoints[pointID];

                        vector d1 = oldPoints[nextPointID] - oldPoints[pointID];

                        // Calculate inverse distance weights
                        // scalar w0 = 1.0/mag(d0);
                        // scalar w1 = 1.0/mag(d1);
                        // Uniform weights
                        scalar w0 = 1.0;
                        scalar w1 = 1.0;

                        vector disp = w0*d0 + w1*d1;

                        scalar sumW = w0 + w1;
                        disp /= sumW;

                        // Keep disp component only in edge direction
                        vector l =
                            oldPoints[nextPointID] - oldPoints[prevPointID];
                        l /= mag(l);
                        disp = sqr(l) & disp;

                        // No smoothing if the angle between the two
                        // adjoining edges is large
                        vector e0 = oldPoints[pointID] - oldPoints[prevPointID];
                        e0 /= mag(e0);
                        vector e1 = oldPoints[nextPointID] - oldPoints[pointID];
                        e1 /= mag(e1);

                        if ((e0 & e1) < edgeAngleDot)
                        {
                            disp = vector::zero;
                        }

                        if (shrink)
                        {
                            // Relaxed motion
                            newPoints[pointID] =
                                oldPoints[pointID] + lambda*disp;

                            residual[pointID] = lambda*mag(disp);
                        }
                        else
                        {
                            // Inflate
                            //newPoints[pointID] -= mu*disp;
                            // They may be smoothed more than once
                            newPoints[pointID] = oldPoints[pointID] - mu*disp;

                            residual[pointID] = mu*mag(disp);
                        }
                    }

                    // Just correct edge point motion
                    vector disp = newPoints[pointID] - oldPoints[pointID];

                    // No smoothing if the angle between the two
                    // adjoining edges is large
                    vector e0 = oldPoints[pointID] - oldPoints[prevPointID];
                    if (mag(e0) < SMALL)
                    {
                        FatalErrorIn("volumeSmoother")
                            << "Two boundary points are coincident"
                            << exit(FatalError);
                    }
                    e0 /= mag(e0);

                    vector e1 = oldPoints[nextPointID] - oldPoints[pointID];
                    if (mag(e1) < SMALL)
                    {
                        FatalErrorIn("volumeSmoother")
                            << "Two boundary points are coincident"
                            << exit(FatalError);
                    }
                    e1 /= mag(e1);

                    if ((e0 & e1) < edgeAngleDot)
                    {
                        disp = vector::zero;
                    }

                    // Keep only motion in the edge direction i.e. points on
                    // an edge can only slide along that edge
                    vector l =
                        oldPoints[nextPointID] - oldPoints[prevPointID];
                    l /= mag(l);

                    // We will increase the lambda parameter on edges to account
                    // for the lower number of neighbours: factor found by
                    // trial-and-error
                    const scalar sFac = 1.0;
                    if (shrink)
                    {
                        // Relaxed motion
                        newPoints[pointID] =
                            oldPoints[pointID] + (sqr(l) & sFac*lambda*disp);
                        residual[pointID] = mag(sqr(l) & sFac*lambda*disp);
                    }
                    else
                    {
                        // inflate
                        newPoints[pointID] =
                            oldPoints[pointID] - (sqr(l) & sFac*mu*disp);
                        residual[pointID] = mag(sqr(l) & sFac*mu*disp);
                    }
                }
            }
        }

        Info<< "Iteration: " << counter << ", max displacement: "
            << max(residual) << endl;

        twoDPointCorrector twoDCorrector(mesh);
        twoDCorrector.correctPoints(newPoints);

        // Fix: the points instance must stored prior to moving the points
        oldInstance() = mesh.pointsInstance();

        mesh.movePoints(newPoints);
        mesh.moving(false);
        mesh.changing(false);
        mesh.setPhi().writeOpt() = IOobject::NO_WRITE;

        if (shrink)
        {
            shrink = false;
        }
        else
        {
            shrink = true;
        }
    }
    while (counter < nCorr);

    return 0;

}

// ************************************************************************* //

} // end of namespace foam
