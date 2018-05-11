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

#include "implicitVolumeSmoother.H"
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
    defineTypeNameAndDebug(implicitVolumeSmoother, 0);
    addToRunTimeSelectionTable
    (
        meshSmoother, implicitVolumeSmoother, dictionary
    );


// * * * * * * * * * * * * Private Data Members  * * * * * * * * * * * * * * //

fvMesh& implicitVolumeSmoother::subMeshToSmooth()
{
    return subsetMeshToSmooth().subMesh();
}


const fvMesh& implicitVolumeSmoother::subMeshToSmooth() const
{
    return subsetMeshToSmooth().subMesh();
}


newFvMeshSubset& implicitVolumeSmoother::subsetMeshToSmooth()
{
    if (subMeshToSmooth_.empty())
    {
        makeSubMeshToSmooth();
    }

    return subMeshToSmooth_();
}


const newFvMeshSubset& implicitVolumeSmoother::subsetMeshToSmooth() const
{
    if (subMeshToSmooth_.empty())
    {
        makeSubMeshToSmooth();
    }

    return subMeshToSmooth_();
}


void implicitVolumeSmoother::makeSubMeshToSmooth() const
{
    if (debug)
    {
        InfoIn("void implicitVolumeSmoother::makeSubMeshToSmooth() const")
            << endl;
    }

    if (!subMeshToSmooth_.empty())
    {
        FatalErrorIn("void implicitVolumeSmoother::makeSubMeshToSmooth() const")
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
                FatalErrorIn("scalar implicitVolumeSmoother::smooth()")
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
        InfoIn("void implicitVolumeSmoother::makeSubMeshToSmooth() const")
            << "Creating the fvMeshSubset" << endl;
    }

    // Select cells with cellsToSmooth set to 1
    subMeshToSmooth_().setLargeCellSubset(cellsToSmooth, 1);
}


tmp<surfaceScalarField> implicitVolumeSmoother::diffusivity() const
{
    // Take a reference to the sub-mesh
    const fvMesh& subMesh = subMeshToSmooth();

    // Diffusivity field for mesh smoothing
    tmp<surfaceScalarField> tdiffusivity
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshDiffusivity",
                subMesh.time().timeName(),
                subMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            subMesh,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    // Calculate the cell volumes field
    const volScalarField V = cellVolumesField();

    // Interpolate the volumes to the faces
    tdiffusivity() = fvc::interpolate(V);

    return tdiffusivity;
}


tmp<volVectorField> implicitVolumeSmoother::cellCentresField() const
{
    // Take a reference to the sub-mesh
    const fvMesh& subMesh = subMeshToSmooth();

    // Cell centres field with fixedValue boundaries
    tmp<volVectorField> tcellCentres
    (
        new volVectorField
        (
            IOobject
            (
                "cellCentresMeshSmoother",
                subMesh.time().timeName(),
                subMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            subMesh,
            dimensionedVector("zero", dimLength, vector::zero),
            "fixedGradient"
        )
    );

    // Set the cell centres internal field
    tcellCentres().internalField() = subMesh.C().internalField();

    // Set the cell centres boundary field
    forAll(tcellCentres().boundaryField(), patchI)
    {
        vectorField& cellCentresP = tcellCentres().boundaryField()[patchI];
        const vectorField& CP = subMesh.C().boundaryField()[patchI];
        forAll(cellCentresP, faceI)
        {
            cellCentresP[faceI] = CP[faceI];
        }

        // if (tcellCentres().boundaryField()[patchI].type() == "fixedGradient")
        // {
        //     fixedGradientFvPatchVectorField& cellCentreP =
        //         refCast<fixedGradientFvPatchVectorField>
        //         (
        //             tcellCentres().boundaryField()[patchI]
        //         );

        //     cellCentreP.gradient() =
        //         subMesh.boundaryMesh()[patchI].faceNormals();
        // }
    }

    return tcellCentres;
}


tmp<volScalarField> implicitVolumeSmoother::cellVolumesField() const
{
    // Take a reference to the sub-mesh
    const fvMesh& subMesh = subMeshToSmooth();

    // Cell volumes field with zero gradient boundaries
    tmp<volScalarField> tcellVolumes
    (
        new volScalarField
        (
            IOobject
            (
                "cellVolumesMeshSmoother",
                subMesh.time().timeName(),
                subMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            subMesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        )
    );

    // Set the cell volumes internal field
    tcellVolumes().internalField() = subMesh.V();

    // Optional: scale the volume field by a cell quality metric
    // In this was, we can make it seem like "bad cells" are smaller than they
    // are, causing them to be expanded, and hopefully improving their quality
    if (dict().lookupOrDefault<Switch>("useCellQuality", false))
    {
        Info<< nl << "Using cell quality field" << endl;

        cellQuality quality(subMesh);

        const word cellMethod = word(dict().lookup("cellQualityMethod"));
        const scalar qualityThreshold =
            readScalar(dict().lookup("cellQualityThreshold"));

        if (cellMethod == "nonOrthogonality")
        {
            scalarField cellNonOrthogonality = quality.nonOrthogonality();
            cellNonOrthogonality /= gMax(cellNonOrthogonality);
            cellNonOrthogonality *= cellNonOrthogonality;
            cellNonOrthogonality = max(cellNonOrthogonality, qualityThreshold);
            cellNonOrthogonality /= qualityThreshold;

            Info<< "    Cell non-orthogonality: "
                << "max = " << max(cellNonOrthogonality)
                << ", average = " << average(cellNonOrthogonality) << endl;

            tcellVolumes().internalField() *= 1.0/cellNonOrthogonality;
        }
        else if (cellMethod == "skewness")
        {
            scalarField cellSkewness = quality.skewness();
            cellSkewness /= gMax(cellSkewness);
            cellSkewness = max(cellSkewness, qualityThreshold);

            Info<< "    Cell skewness: "
                << "max = " << max(cellSkewness)
                << ", average = " << average(cellSkewness) << endl;

            tcellVolumes().internalField() *= 1.0/cellSkewness;
        }
        else
        {
            FatalErrorIn
            (
                "tmp<volScalarField> implicitVolumeSmoother::"
                "cellVolumesField() const"
            )   << "cellQualityMethod " << cellMethod
                << " is unknown! Options are nonOrthogonality or skewness"
                << abort(FatalError);
        }

        // Reset the cells field values at the boundary as these cells are
        // typically fine
        if (dict().found("boundaryCellScaleFactor"))
        {
            Info<< "Overwriting boundary cell values" << endl;

            scalarField& cellVolumesI = tcellVolumes().internalField();
            const scalarField& VI = subMesh.V();
            const scalar boundaryCellScaleFactor =
                readScalar(dict().lookup("boundaryCellScaleFactor"));
            forAll(subMesh.boundaryMesh(), patchI)
            {
                if (!subMesh.boundaryMesh()[patchI].coupled())
                {
                    const labelList& faceCells =
                        subMesh.boundaryMesh()[patchI].faceCells();

                    forAll(faceCells, fI)
                    {
                        const label cellID = faceCells[fI];
                        cellVolumesI[cellID] =
                            boundaryCellScaleFactor*VI[cellID];
                    }
                }
            }
        }
    }

    // Update the boundaries
    tcellVolumes().correctBoundaryConditions();

    return tcellVolumes;
}


void implicitVolumeSmoother::mapVolToPoint
(
    const volVectorField& cellMotionD,
    pointVectorField& pointMotionD
)
{
    if (debug)
    {
        Info<< "    Mapping cellMotionD to pointMotionD" << endl;
    }

    // One option to perform an interpolation; however, this can lead to
    // oscillations occuring in the vertices e.g. cells become skewed and this
    // is not noticed by the cell centres. So a better option is to force the
    // vertices to lie at the average point of the cell centres. This enforces
    // the volume averaging but also detects and removed cell skewness

    // Option 1: interpolation - this can result in skewed cells
    // newLeastSquaresVolPointInterpolation::New(subMesh).interpolate
    // (
    //     cellMotionD,
    //     pointMotionD
    // );

    // Option 2: force vertices to lie at the average position of the cell
    // centres

    // Take references for efficiency
    const fvMesh& subMesh = cellMotionD.mesh();
    vectorField& pointMotionDI = pointMotionD.internalField();
    const vectorField& cellMotionDI = cellMotionD.internalField();
    const labelListList& pointCells = subMesh.pointCells();
    const pointField oldPoints = subMesh.points();
    const volVectorField cellCentres = cellCentresField();
    const vectorField& cellCentresI = cellCentres.internalField();


    //
    // Calculate point weights to interpolate from the cell centres to the
    // points: we will use a simple arithmetic average, but we need to be
    // careful in parallel
    //


    // Initialise point weights for interpolating from the cells to the points
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
            "MapVolToPointSumWeights",
            subMesh.polyMesh::instance(),
            subMesh
        ),
        pointMesh::New(subMesh),
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Set un-normalised weights to 1.0 (i.e. arithmetic average)
    forAll(weights, pI)
    {
        scalarList& pw = weights[pI];
        const labelList& pcp = pointCells[pI];

        forAll(pcp, pointCellI)
        {
            // We are taking an arithmetic average
            pw[pointCellI] = 1.0;
            sumWeights[pI] += 1.0;
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


    //
    // Calculate the new point positions using the point weights
    //

    scalar oscillCorr =
        dict().lookupOrDefault<scalar>("oscillationCorrection", 0.1);

    if (oscillCorr > 1.0)
    {
        WarningIn
        (
            "void implicitVolumeSmoother::mapVolToPoint\n"
            "(\n"
            "    const volVectorField& cellMotionD,\n"
            "    pointVectorField& pointMotionD\n"
            ")"
        )   << "Limiting oscillationCorrection to 1.0" << endl;

        oscillCorr = 1.0;
    }

    forAll(pointMotionDI, pI)
    {
        vector cuPointMotion = vector::zero;
        const labelList& curPointCells = pointCells[pI];
        const scalarList& curWeights = weights[pI];

        forAll(curPointCells, pcI)
        {
            const label cellID = curPointCells[pcI];

            // Simple interpolation
            // Susceptible to checkerboarding oscillations
            cuPointMotion += curWeights[pcI]*cellMotionDI[cellID];

            // Correction to force the point to lie at the average position of
            // the cell centres to remove checkerboard-like oscillations in the
            // mesh
            cuPointMotion +=
                oscillCorr*curWeights[pcI]
               *(
                   cellCentresI[cellID] - oldPoints[pI]
                );
        }

        pointMotionDI[pI] = cuPointMotion;
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


void implicitVolumeSmoother::calculatePointMotionImplicitly
(
    pointVectorField& pointMotionD
)
{
    // Take a reference to the subMesh
    fvMesh& subMesh = subMeshToSmooth();

    // Create the cell displacement field
    // We will set all the boundaries to slip
    volVectorField cellMotionD
    (
        IOobject
        (
            "cellMotionD",
            subMesh.time().timeName(),
            subMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        subMesh,
        dimensionedVector("zero", dimLength, vector::zero),
        "slip"
    );

    // Take a reference to the interalFields for efficiency and convenience
    vectorField& pointMotionDI = pointMotionD.internalField();
    vectorField& cellMotionDI = cellMotionD.internalField();

    // Calculate the diffusivity field for the smooth: this governs the
    // smoothed shape and size of the cells
    const surfaceScalarField alpha = diffusivity();

    // The old "not smoothed" cell centres provide the driving force source
    // term for smoothing
    const volScalarField V = cellVolumesField();

    // Non-orthogonal correction loop
    for (int iNonOrthoCorr = 0; iNonOrthoCorr < 3; iNonOrthoCorr++)
    {
        // Create the Laplacian smoothing equation where the weights are a
        // function of the cell volumes
        // The grad of the volumes is the driving force
        fvVectorMatrix cellMotionDEqn
        (
            fvm::laplacian(V, cellMotionD)
          + fvc::grad(V)
        );

        // Solve the linear system
        cellMotionDEqn.solve();
    }

    // Limit motion to a fraction of the Courant number
    if (dict().found("cellMotionLimit"))
    {
        const scalar cellMotionLimit
            (
                readScalar(dict().lookup("cellMotionLimit"))
            );

        const scalarField& VI = V.internalField();
        forAll(cellMotionDI, cellI)
        {
            vector& curCellMotionD = cellMotionDI[cellI];
            const scalar curMagCellMotionD = mag(curCellMotionD);
            const scalar delta = cellMotionLimit*Foam::cbrt(VI[cellI]);

            if (curMagCellMotionD > delta)
            {
                curCellMotionD = delta*curCellMotionD/curMagCellMotionD;
            }
        }

        forAll(cellMotionD.boundaryField(), patchI)
        {
            if (!cellMotionD.boundaryField()[patchI].coupled())
            {
                vectorField& cellMotionDP =
                    cellMotionD.boundaryField()[patchI];
                const scalarField& VP = V.boundaryField()[patchI];

                forAll(cellMotionDP, faceI)
                {
                    vector& curCellMotionD = cellMotionDP[faceI];
                    const scalar curMagCellMotionD = mag(curCellMotionD);
                    const scalar delta =
                        cellMotionLimit*Foam::cbrt(VP[faceI]);

                    if (curMagCellMotionD > delta)
                    {
                        curCellMotionD =
                            delta*curCellMotionD/curMagCellMotionD;
                    }
                }
            }
        }

        cellMotionD.correctBoundaryConditions();
    }

    if (debug)
    {
        Info<< "    max(mag(cellMotionD)): " << max(mag(cellMotionD))
            << endl;
    }

    // We should move the subMesh here so that the cell volumes are updated
    // Or we could equivalently include the Jacobian in the equation
    // Ah the subMesh only takes a reference to the base mesh so we need to
    // move the base mesh to update the volumes: also, remember that we need
    // to update the boundary normals etc. It is easier to move the mesh,
    // though this could be a place to improve efficiency

    // Interpolate the cellMotionD field to the mesh vertices

    // Map motionD from the cell centres to the vertices
    mapVolToPoint(cellMotionD, pointMotionD);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
implicitVolumeSmoother::implicitVolumeSmoother
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    meshSmoother(mesh, dict),
    subMeshToSmooth_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

implicitVolumeSmoother::~implicitVolumeSmoother()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void implicitVolumeSmoother::clearOut()
{
    subMeshToSmooth_.clear();
}


scalar implicitVolumeSmoother::smooth()
{
    if (debug)
    {
        InfoIn("scalar implicitVolumeSmoother::smooth()")
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

        calculatePointMotionImplicitly(pointMotionD);

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
