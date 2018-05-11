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

#include "explicitVolumeSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "twoDPointCorrector.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "ZoneID.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcDiv.H"
#include "fixedGradientFvPatchFields.H"
#include "cellQuality.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(explicitVolumeSmoother, 0);
    addToRunTimeSelectionTable
    (
        meshSmoother, explicitVolumeSmoother, dictionary
    );


// * * * * * * * * * * * * Private Data Members  * * * * * * * * * * * * * * //

fvMesh& explicitVolumeSmoother::subMeshToSmooth()
{
    return subsetMeshToSmooth().subMesh();
}


const fvMesh& explicitVolumeSmoother::subMeshToSmooth() const
{
    return subsetMeshToSmooth().subMesh();
}


newFvMeshSubset& explicitVolumeSmoother::subsetMeshToSmooth()
{
    if (subMeshToSmooth_.empty())
    {
        makeSubMeshToSmooth();
    }

    return subMeshToSmooth_();
}


const newFvMeshSubset& explicitVolumeSmoother::subsetMeshToSmooth() const
{
    if (subMeshToSmooth_.empty())
    {
        makeSubMeshToSmooth();
    }

    return subMeshToSmooth_();
}


void explicitVolumeSmoother::makeSubMeshToSmooth() const
{
    if (debug)
    {
        InfoIn("void explicitVolumeSmoother::makeSubMeshToSmooth() const")
            << endl;
    }

    if (!subMeshToSmooth_.empty())
    {
        FatalErrorIn("void explicitVolumeSmoother::makeSubMeshToSmooth() const")
            << "pointer already set!" << abort(FatalError);
    }

    // Take a reference to the mesh for efficiency
    const fvMesh& mesh = this->mesh();

    // Read the list of cell zones to smooth from the dict; if none are
    // specified then we will smooth the entire mesh

    // List indicating cells to be smoothed (1) and cells to not smooth (0)
    labelList cellsToSmooth(mesh.nCells(), 0);

    if (dict().found("cellZones"))
    {
        // Read the cell zones to smooth
        const wordList cellZonesNames = wordList(dict().lookup("cellZones"));

        forAll(cellZonesNames, czI)
        {
            // Check the cell zone exist
            const ZoneID<cellZone> cellZoneID =
                ZoneID<cellZone>(cellZonesNames[czI], mesh.cellZones());

            if (!cellZoneID.active())
            {
                FatalErrorIn("scalar explicitVolumeSmoother::smooth()")
                    << "cellZone " << cellZonesNames[czI] << " not found!"
                    << abort(FatalError);
            }

            // Set cells in zone to be smoothed
            Info<< "Cell zone " << cellZonesNames[czI] << " will be smoothed"
                << endl;
            const labelList& curCellZone = mesh.cellZones()[cellZoneID.index()];
            forAll(curCellZone, cI)
            {
                cellsToSmooth[curCellZone[cI]] = 1;
            }
        }
    }
    else
    {
        // Smooth all the cells
        Info<< "All cells in the mesh will be smoothed"
            << endl;
        cellsToSmooth = 1.0;
    }


    // Create a subMesh containing only the cells to be smoothed
    subMeshToSmooth_.set
    (
        new newFvMeshSubset
        (
            IOobject
            (
                "cellsToSmooth",
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    if (debug)
    {
        InfoIn("void explicitVolumeSmoother::makeSubMeshToSmooth() const")
            << "Creating the fvMeshSubset" << endl;
    }

    // Select cells with cellsToSmooth set to 1
    subMeshToSmooth_().setLargeCellSubset(cellsToSmooth, 1);
}


void explicitVolumeSmoother::calculatePointMotionExplicitly
(
    pointVectorField& pointMotionD
)
{
    // Read settings specific to the explicit algorithm
    const scalar lambda(readScalar(dict().lookup("lambda")));
    const Switch simpleAverageWeights(dict().lookup("simpleAverageWeights"));

    // Take a reference to the interalFields for efficiency and convenience
    const fvMesh& subMesh = subMeshToSmooth();
    vectorField& pointMotionDI = pointMotionD.internalField();
    const scalarField& oldV = subMesh.V();
    const labelListList& pointCells = subMesh.pointCells();
    const labelListList& pointFaces = subMesh.pointFaces();
    const pointField& oldPoints = subMesh.points();
    const volVectorField& oldC = subMesh.C();
    const vectorField& oldCI = subMesh.C().internalField();


    // Calculate the weights for each point based on the adjacent cell volumes

    // Initialise point weights
    scalarListList weights(subMesh.nPoints());
    forAll(weights, pointI)
    {
        weights[pointI].setSize(pointCells[pointI].size());
    }

    // We need to create a field to store the sum of the weights so that points
    // on processor boundaries can be corrected
    pointScalarField sumWeights
    (
        IOobject
        (
            "volumeSmootherExplicitWeights",
            subMesh.polyMesh::instance(),
            subMesh
        ),
        pointMesh::New(subMesh),
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Set un-normalised weights to cell volume
    forAll(weights, pI)
    {
        scalarList& pw = weights[pI];
        const labelList& pcp = pointCells[pI];

        forAll(pcp, pointCellI)
        {
            // Point cell index
            const label curPointCellID = pcp[pointCellI];

            if (simpleAverageWeights)
            {
                // Set the un-normalised point weight
                pw[pointCellI] = 1.0;

                // Increment the weights sum
                sumWeights[pI] += 1.0;
            }
            else
            {
                // Set the un-normalised point weight
                pw[pointCellI] = oldV[curPointCellID];

                // Increment the weights sum
                sumWeights[pI] += oldV[curPointCellID];
            }
        }
    }

    // Sync weights in parallel
    forAll(sumWeights.boundaryField(), patchI)
    {
        if (sumWeights.boundaryField()[patchI].coupled())
        {
            sumWeights.boundaryField()[patchI].initAddField();
        }
    }

    forAll(sumWeights.boundaryField(), patchI)
    {
        if (sumWeights.boundaryField()[patchI].coupled())
        {
            sumWeights.boundaryField()[patchI].addField
            (
                sumWeights.internalField()
            );
        }
    }

    // Normalise the weights
    forAll(weights, pI)
    {
        scalarList& pw = weights[pI];

        forAll(pw, pointCellI)
        {
            pw[pointCellI] /= sumWeights[pI];
        }
    }


    // Reset pointMotionDI to zero
    pointMotionDI = vector::zero;


    // Calculate the point motion
    forAll(pointMotionDI, pointI)
    {
        // Take references
        vector& curPointMotionD = pointMotionDI[pointI];
        const scalarList& curPointWeights = weights[pointI];
        const labelList& curPointCells = pointCells[pointI];

        // Loop through all neighbour cell centres for this point
        forAll(curPointCells, pcI)
        {
            // Neighbour cell centre index
            const label neiCellID = curPointCells[pcI];

            // Vector from the previous cell centre to the previous point
            const vector d = oldCI[neiCellID] - oldPoints[pointI];

            // Add contribution to the point motion
            curPointMotionD += curPointWeights[pcI]*d;
        }

        // Relax the motion motion
        curPointMotionD *= lambda;
    }

    // Add contributions from other processors
    forAll(pointMotionD.boundaryField(), patchI)
    {
        if (pointMotionD.boundaryField()[patchI].coupled())
        {
            pointMotionD.boundaryField()[patchI].initAddField();
        }
    }

    forAll(pointMotionD.boundaryField(), patchI)
    {
        if (pointMotionD.boundaryField()[patchI].coupled())
        {
            pointMotionD.boundaryField()[patchI].addField
            (
                pointMotionD.internalField()
            );
        }
    }

    // Update coupled and constrained boundaries
    pointMotionD.correctBoundaryConditions();


    // Correct boundary motion
    correctBoundaryMotion(pointMotionD);

    // Correct point motions for points shared by more than one processor
    pointMotionD.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
explicitVolumeSmoother::explicitVolumeSmoother
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    meshSmoother(mesh, dict),
    subMeshToSmooth_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

explicitVolumeSmoother::~explicitVolumeSmoother()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void explicitVolumeSmoother::clearOut()
{
    subMeshToSmooth_.clear();
}


scalar explicitVolumeSmoother::smooth()
{
    if (debug)
    {
        InfoIn("scalar explicitVolumeSmoother::smooth()")
            << endl;
    }

    Info<< "Smoothing the mesh" << endl;

    // Get the sub-mesh containing the cells to be smoothed
    // This could be the whole mesh or specified cells zones
    fvMesh& subMesh = subMeshToSmooth();

    // Create the point mesh
    const pointMesh& pMesh = pointMesh::New(subMesh);

    // Create the point motion field
    pointVectorField pointMotionD
    (
        IOobject
        (
            "pointMotionD",
            subMesh.time().timeName(),
            subMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Take a reference to the interalFields for efficiency and convenience
    vectorField& pointMotionDI = pointMotionD.internalField();

    // Perform smoothing iterations
    int i = 0;
    const int nCorr(dict().lookupOrDefault<int>("nCorrectors", 1));
    blockLduMatrix::debug = 0;
    do
    {
        Info<< "    Iteration " << (i + 1);

        calculatePointMotionExplicitly(pointMotionD);

        // TODO: try adding a limiter to limit the point motion to a fraction of
        // edge length: this may stablise it
        // Also, what happens when a volume is negative? Check this on simple
        // cases


        Info<< ", maximum point motion is " << gMax(mag(pointMotionDI))
            << endl;

        // Sub-mesh new points
        pointField newPointsSubMesh = subMesh.points() + pointMotionDI;

        // Apply corrections for 2-D meshes
        twoDPointCorrector twoDCorrectorSubMesh(subMesh);
        twoDCorrectorSubMesh.correctPoints(newPointsSubMesh);

        const scalarField oldVSubMesh = subMesh.V();

        // It is important to set the old points otherwise the swept volumes
        // will be wrong
        subMesh.setOldPoints(subMesh.points());

        // Move the subMesh
        subMesh.movePoints(newPointsSubMesh);

        // Disable flags saying that the mesh is moved as we do not need to
        // store the mesh motion fluxes
        subMesh.moving(false);
        subMesh.changing(false);
        subMesh.setPhi().writeOpt() = IOobject::NO_WRITE;

        // Take a reference to the main mesh: be careful not to confuse the
        // sub-mesh and main mesh
        fvMesh& mesh = this->mesh();

        // Map the point mesh displacements from the sub-mesh to the main mesh
        const labelList& pointMap = subsetMeshToSmooth().pointMap();

        // Set newPoints field by mapping the sub-mesh point displacement to the
        // main mesh
        pointField newPoints = mesh.points();
        forAll(pointMap, subMeshPI)
        {
            const label pointID = pointMap[subMeshPI];
            newPoints[pointID] += pointMotionDI[subMeshPI];
        }

        // Apply corrections for 2-D meshes
        twoDPointCorrector twoDCorrector(mesh);
        twoDCorrector.correctPoints(newPoints);

        // Store the old cell volumes
        const scalarField volOld = scalarField(mesh.V());

        // It is important to set the old points otherwise the swept volumes
        // will be wrong
        mesh.setOldPoints(mesh.points());

        // Move the mesh and store the swept face volumes
        const scalarField sweptVol = mesh.movePoints(newPoints);

        if (Switch(dict().lookup("advectFields")))
        {
            // Store the new cell volumes
            const scalarField volNew = scalarField(mesh.V());

            // Advect the volFields
            advectFields(sweptVol, volOld, volNew);
        }

        // Disable flags saying that the mesh is moved as we do not need to
        // store the mesh motion fluxes
        mesh.moving(false);
        mesh.changing(false);
        mesh.setPhi().writeOpt() = IOobject::NO_WRITE;

        // Clear mesh data
        mesh.clearOut();
        subMesh.clearOut();
    }
    while (++i < nCorr);

    blockLduMatrix::debug = 0;

    return 0;
}

// ************************************************************************* //

} // end of namespace foam
