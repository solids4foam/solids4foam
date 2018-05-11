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

#include "taubinSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "twoDPointCorrector.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(taubinSmoother, 0);
    addToRunTimeSelectionTable(meshSmoother, taubinSmoother, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
taubinSmoother::taubinSmoother
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    meshSmoother(mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

taubinSmoother::~taubinSmoother()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate weighing for given distance between two points
scalar taubinSmoother::weight(const word& weighingFunctionName, const vector& d)
{
    scalar w = 1.0;

    if (weighingFunctionName == "arithmeticAverage")
    {
        w = 1.0;
    }
    else if (weighingFunctionName == "inverseDistance")
    {
        w = 1.0/mag(d);
    }
    else if (weighingFunctionName == "inverseDistanceSquared")
    {
        w = 1.0/Foam::sqrt(mag(d));
    }
    else if (weighingFunctionName == "inverseDistanceSquareRoot")
    {
        w = 1.0/Foam::sqrt(mag(d));
    }
    else
    {
        FatalErrorIn
        (
            "scalar taubinSmoother::weight"
            "(const word& weighingFunctionName, const vector& d)"
        )   << "weighingFunctionName " << weighingFunctionName
            << " is not known!" << abort(FatalError);
    }

    return w;
}


// Smoothing function
scalar taubinSmoother::smooth()
{

    const scalar lambda(readScalar(dict().lookup("lambda")));
    const scalar mu(readScalar(dict().lookup("mu")));
    const int nCorr(readInt(dict().lookup("nCorrectors")));

    boolList cellZonePoints(mesh().nPoints(), false);

    int edgeAngle
    (
        dict().lookupOrDefault<int>
        (
            "edgeAngle",
            45
        )
    );
    const scalar edgeAngleDot =
        Foam::cos(edgeAngle*mathematicalConstant::pi/180.0);


    bool independentPatches
    (
        dict().lookupOrDefault<Switch>
        (
            "independentPatches",
            false
        )
    );

    word cellZoneName
    (
        dict().lookupOrDefault<word>
        (
            "cellZone",
            "All"
        )
    );

    word weighingFunctionName
    (
        dict().lookupOrDefault<word>
        (
            "weighingFunction",
            "arithmeticAverage"
        )
    );
    Info<< "        weighingFunction: " << weighingFunctionName << endl;

    label cellZoneID = -1;

    if (cellZoneName != "All")
    {

        cellZoneID = mesh().cellZones().findZoneID(cellZoneName);

        if (cellZoneID == -1)
        {
            FatalError
                << "cellZone " << cellZoneName << " not found!"
                << abort(FatalError);
        }

        Info<< "Smoothing points in cellZone " << cellZoneName
            << nl << endl;

        // Mark all points that are part of any cell in the cellZone
        const labelListList& cellPoints = mesh().cellPoints();
        forAll(mesh().cellZones()[cellZoneID], cellI)
        {
            const label cellID = mesh().cellZones()[cellZoneID][cellI];
            forAll(cellPoints[cellID], cpI)
            {
                cellZonePoints[cellPoints[cellID][cpI]] = true;
            }
        }
    }
    else
    {
        cellZonePoints = true;
    }

    // Designate boundary points
    boolList boundaryPoints(mesh().nPoints(), false);

    forAll(mesh().boundaryMesh(), patchI)
    {
        const labelList& meshPointsOfBoundaryPatch =
            mesh().boundaryMesh()[patchI].meshPoints();
        forAll(meshPointsOfBoundaryPatch, pointI)
        {
            boundaryPoints[meshPointsOfBoundaryPatch[pointI]] = true;
        }
    }

    // Calculate point-point weights using inverse distances
    // Note: these weights are divided by the sum of the weights during
    // smoothing
    const labelListList& pointPoints = mesh().pointPoints();

    scalarField residual(mesh().nPoints(), 0);
    label counter = 0;
    bool shrink = true;
    pointField oldPoints = mesh().points();
    pointField newPoints = oldPoints;
    do
    {
        oldPoints = newPoints;

        if (!shrink)
        {
            counter++;
        }

        forAll(newPoints, pointI)
        {
            // Smooth internal points or all points if dependent patches
            if
            (
                (!boundaryPoints[pointI] || !independentPatches)
                && cellZonePoints[pointI]
            )
            {

                vector curNewPoint = vector::zero;

                scalar sumW = 0;

                forAll(pointPoints[pointI], ppI)
                {
                    vector d =
                        oldPoints[pointPoints[pointI][ppI]] - oldPoints[pointI];

                    scalar w = weight(weighingFunctionName, d);

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

        // Correct patch motion and optionally smooth independently
        forAll(mesh().boundaryMesh(), patchI)
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const pointField& pointNormals =
                mesh().boundaryMesh()[patchI].pointNormals();

            forAll(meshPoints, pointI)
            {
                label pointID = meshPoints[pointI];

                if (independentPatches && cellZonePoints[pointID])
                {
                    vector curNewPoint = vector::zero;

                    scalar sumW = 0;

                    forAll(pointPoints[pointID], ppI)
                    {
                        label otherPointID = pointPoints[pointID][ppI];

                        vector d = oldPoints[otherPointID] - oldPoints[pointID];

                        scalar w = weight(weighingFunctionName, d);

                        curNewPoint += w*d;
                        sumW += w;
                    }

                    curNewPoint /= sumW;

                    curNewPoint += oldPoints[pointID];

                    vector disp = curNewPoint - oldPoints[pointID];

                    // Remove point normal component, so point motion is along
                    // patch
                    disp = (I - sqr(pointNormals[pointI])) & disp;

                    if (shrink)
                    {
                        // Relaxed motion
                        newPoints[pointID] = oldPoints[pointID] + lambda*disp;

                        residual[pointID] = lambda*mag(disp);
                    }
                    else
                    {
                        // Inflate
                        newPoints[pointID] = oldPoints[pointID] - mu*disp;

                        residual[pointID] = mu*mag(disp);
                    }
                }
                else
                {
                    // Keep only motion in the patch plane i.e. remove motion
                    // normal to the patch
                    newPoints[pointID] =
                        oldPoints[pointID]
                        +
                        (
                            (I - sqr(pointNormals[pointI])) &
                            (newPoints[pointID] - oldPoints[pointID])
                        );

                    residual[pointID] =
                        mag
                        (
                            (I - sqr(pointNormals[pointI])) &
                            (newPoints[pointID] - oldPoints[pointID])
                        );
                }
            }

            // Correct patch boundary edge motion, and optionally smooth
            // independently

            const labelListList& edgeLoops =
                mesh().boundaryMesh()[patchI].edgeLoops();

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

                    if (independentPatches && cellZonePoints[pointID])
                    {
                        vector d0 = oldPoints[prevPointID] - oldPoints[pointID];

                        vector d1 = oldPoints[nextPointID] - oldPoints[pointID];

                        // Calculate inverse distance weights
                        scalar w0 = weight(weighingFunctionName, d0);
                        scalar w1 = weight(weighingFunctionName, d1);

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
                    else
                    {
                        // Just correct edge point motion
                        vector disp = newPoints[pointID] - oldPoints[pointID];

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

                        // Keep only motion in the edge direction i.e. points on
                        // an edge can only slide along that edge
                        vector l =
                            oldPoints[nextPointID] - oldPoints[prevPointID];
                        l /= mag(l);

                        // newPoints[pointID] =
                        //     oldPoints[pointID] + (sqr(l) & lambda*disp);

                        // Bug-fix 08-Dec-15 PC: we forgot to inflate the edges
                        // This only makes a small difference
                        if (shrink)
                        {
                            // Relaxed motion
                            newPoints[pointID] =
                                oldPoints[pointID] + (sqr(l) & lambda*disp);
                            residual[pointID] = mag(sqr(l) & lambda*disp);
                        }
                        else
                        {
                            // inflate
                            newPoints[pointID] =
                                oldPoints[pointID] - (sqr(l) & mu*disp);
                            residual[pointID] = mag(sqr(l) & mu*disp);
                        }
                    }
                }
            }
        }

/*        if (writeIters)
        {
            twoDPointCorrector twoDCorrector(mesh());
            twoDCorrector.correctPoints(newPoints);

            mesh().movePoints(newPoints);
        }
*/
        Info<< "Iteration: " << counter << ", max displacement: "
            << max(residual) << endl;

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

    twoDPointCorrector twoDCorrector(mesh());
    twoDCorrector.correctPoints(newPoints);

    Info<< "The points instance is: " << mesh().pointsInstance() << endl;

    // Fix: the points instance must stored prior to moving the points
    oldInstance() = mesh().pointsInstance();

    mesh().movePoints(newPoints);
    mesh().moving(false);
    mesh().changing(false);
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

    return 0;

}

// ************************************************************************* //

} // end of namespace foam
