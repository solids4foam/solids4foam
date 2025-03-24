/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    perturbMeshPoints

Description
    Add a random perturbation to each mesh points.The perturbation of a point
    is calculated as a scale factor times the local minimum edge length.

    By default, boundary points are slide along the patch. Patches which should
    not move can be defined via the fixedPatches entry.

    In OpenFOAM.com, points on feature edges are not moved, where feature edges
    are defined by the minimum cosine angle (minCos). This feature is not
    currently implemented with OpenFOAM.org and foam-extend.

    The inputs are defined in $FOAM_CASE/system/perturbMeshPointsDict, and
    consist of a seed (for the random number generator) and a scaling factor
    to scale the perturbations. The scaling factor is a vector to allow
    different scalings in different directions; for example, for 2-D, the Z
    component should be set to 0.0.

    If there are negative or small cell volumes after moving the points, the
    local motion is reduced by the factor beta (defaults to 0.8) and the motion
    is performed again. The maximum number of corrections iterations is set with
    maxIter (defaults to 1000). A small volume is defined as the factor
    minCellVol (defaults to 0.1) times the original cel volume.

    This utility is useful for creating distorted grids for testing
    discretisations.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "Random.H"
#include "twoDPointCorrector.H"
#include "unitConversion.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "primitiveMeshTools.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef OPENFOAM_COM

// Copied from Foam::polyDualMesh in OpenFOAM-v2312
void calcFeatures
(
    const polyMesh& mesh,
    const scalar featureCos,
    labelList& featureEdges,
    labelList& featurePoints
)
{
    // Create big primitivePatch for all outside.
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // For ease of use store patch number per face in allBoundary.
    labelList allRegion(allBoundary.size());

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, i)
        {
            allRegion[i + pp.start() - mesh.nInternalFaces()] = patchi;
        }
    }


    // Calculate patch/feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelListList& edgeFaces = allBoundary.edgeFaces();
    const vectorField& faceNormals = allBoundary.faceNormals();
    const labelList& meshPoints = allBoundary.meshPoints();

    bitSet isFeatureEdge(edgeFaces.size(), false);

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() != 2)
        {
            // Non-manifold. Problem?
            const edge& e = allBoundary.edges()[edgeI];

            WarningInFunction
                << meshPoints[e[0]] << ' ' << meshPoints[e[1]]
                << "  coords:" << mesh.points()[meshPoints[e[0]]]
                << mesh.points()[meshPoints[e[1]]]
                << " has more than 2 faces connected to it:"
                << eFaces.size() << endl;

            isFeatureEdge.set(edgeI);
        }
        else if (allRegion[eFaces[0]] != allRegion[eFaces[1]])
        {
            isFeatureEdge.set(edgeI);
        }
        else if
        (
            (faceNormals[eFaces[0]] & faceNormals[eFaces[1]])
          < featureCos
        )
        {
            isFeatureEdge.set(edgeI);
        }
    }


    // Calculate feature points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    const labelListList& pointEdges = allBoundary.pointEdges();

    DynamicList<label> allFeaturePoints(pointEdges.size());

    forAll(pointEdges, pointi)
    {
        const labelList& pEdges = pointEdges[pointi];

        label nFeatEdges = 0;

        forAll(pEdges, i)
        {
            if (isFeatureEdge.test(pEdges[i]))
            {
                ++nFeatEdges;
            }
        }
        if (nFeatEdges > 2)
        {
            // Store in mesh vertex label
            allFeaturePoints.append(allBoundary.meshPoints()[pointi]);
        }
    }
    featurePoints.transfer(allFeaturePoints);

    // Get all feature edges.
    labelList meshEdges
    (
        allBoundary.meshEdges
        (
            mesh.edges(),
            mesh.cellEdges(),
            SubList<label>
            (
                mesh.faceOwner(),
                allBoundary.size(),
                mesh.nInternalFaces()
            )
        )
    );

    DynamicList<label> allFeatureEdges(isFeatureEdge.size());
    forAll(isFeatureEdge, edgeI)
    {
        if (isFeatureEdge.test(edgeI))
        {
            // Store in mesh edge label.
            allFeatureEdges.append(meshEdges[edgeI]);
        }
    }
    featureEdges.transfer(allFeatureEdges);
}

#endif // OPENFOAM_COM


#ifdef OPENFOAM_NOT_EXTEND

// Modified form OpenFOAM-v2312 primitiveMeshCheck.C
label numSevereNonOrthoFaces(const fvMesh& mesh)
{
    // Calculate the mesh orthogonality
    tmp<scalarField> tortho = primitiveMeshTools::faceOrthogonality
    (
        mesh,
        mesh.faceAreas(),
        mesh.cellCentres()
    );
    const scalarField& ortho = tortho();

    // Severe nonorthogonality threshold
    const scalar nonOrthThreshold = 70;
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(nonOrthThreshold));

    label severeNonOrth = 0;

    forAll(ortho, facei)
    {
        if (ortho[facei] < severeNonorthogonalityThreshold)
        {
            severeNonOrth++;
        }
    }

    return severeNonOrth;
}

#endif // OPENFOAM_NOT_EXTEND


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    argList::noParallel();

    // Read dictionary
    Info<< "Reading perturbMeshPointsDict dictionary" << nl << endl;
    IOdictionary perturbDict
    (
        IOobject
        (
            "perturbMeshPointsDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read inputs
    const int seed(readInt(perturbDict.lookup("seed")));
    const vector scaleFactor(perturbDict.lookup("scaleFactor"));
    const wordList fixedPatchesList(perturbDict.lookup("fixedPatches"));
    const scalar beta(perturbDict.lookupOrDefault("beta", 0.8));
    const int maxIter(perturbDict.lookupOrDefault("maxIter", 1000));
    const scalar minCellVol(perturbDict.lookupOrDefault("minCellVol", 0.1));
#ifdef OPENFOAM_COM
    const scalar minCos(readScalar(perturbDict.lookup("minCos")));
    const Switch Gaussian(perturbDict.lookup("Gaussian"));
#endif

    // Convert fixedPatches list to a set
    HashSet<word> fixedPatches;
    forAll(fixedPatchesList, pI)
    {
        fixedPatches.insert(fixedPatchesList[pI]);
    }

    // Create random number generator
    Random rnd(seed);

    // Store original points
    const pointField oldPoints = mesh.points();

    // Store the original cell volumes times the minCellVol factor (e.g. 0.1)
    // We will not allow the volume of each cell to become less than this
    const scalarField minV(minCellVol*mesh.V());

    // Calculate and store a copy of the mesh indexing as this will not change
    // and we don't want to recalculate it
    const labelListList cellCells(mesh.cellCells());
    const labelListList pointCells(mesh.pointCells());
    const labelListList pointPoints(mesh.pointPoints());

    // Calculate the minimum edge length connected to each point
    const labelListList& pointEdges = mesh.pointEdges();
    const edgeList& edges = mesh.edges();
    scalarField minEdgeLength(oldPoints.size(), GREAT);
    forAll(minEdgeLength, pI)
    {
        const labelList& curPointEdges = pointEdges[pI];
        scalar minLen = GREAT;
        forAll(curPointEdges, peI)
        {
            const label curEdgeID = curPointEdges[peI];
            const edge& curEdge = edges[curEdgeID];

            minLen = min(minLen, curEdge.mag(oldPoints));
        }

        minEdgeLength[pI] = minLen;
    }

    // Store the original minEdgeLength
    const scalarField oldMinEdgeLength(minEdgeLength);

    // Calculate a mask to identify fixed points
    boolList fixedPoint(oldPoints.size(), false);
    forAll(mesh.boundary(), patchI)
    {
        const word& patchName = mesh.boundary()[patchI].name();
        if (fixedPatches.found(patchName))
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, mpI)
            {
                const label pointID = meshPoints[mpI];
                fixedPoint[pointID] = true;
            }
        }
    }

#ifdef OPENFOAM_COM
    // Mark all feature points as fixed
    {
        // Calculate feature points and edges
        labelList featureEdges;
        labelList featurePoints;
        calcFeatures
        (
            mesh, minCos, featureEdges, featurePoints
        );

        forAll(featureEdges, eI)
        {
            const label curEdgeID = featureEdges[eI];
            const edge& curEdge = edges[curEdgeID];

            // Fix start and end points of the edge
            fixedPoint[curEdge.start()] = true;
            fixedPoint[curEdge.end()] = true;
        }
    }
#endif // OPENFOAM_COM

    // Perform loop
    //     1. Apply mesh motion based on the scaled local minimum edge length
    //     2. If checkMesh fails, reduce the local minimum edge length by the
    //        factor beta (e.g., beta = 0.8)
    //     3. Exit loop if checkMesh passes or the maximum number of iterations
    //        is reached
    bool validMesh = false;
    int iter = 0;
    int nRndReset = 0;
    do
    {
        Info<< "Iteration = " << iter << endl;

        // Calculate new points
        pointField newPoints(oldPoints);

        forAll(newPoints, pointI)
        {
            if (!fixedPoint[pointI])
            {
    #ifdef OPENFOAM_COM
                if (Gaussian)
                {
                    // Gaussian distribution
                    newPoints[pointI] +=
                        minEdgeLength[pointI]
                       *vector
                        (
                            scaleFactor.x()*rnd.GaussNormal<scalar>(),
                            scaleFactor.y()*rnd.GaussNormal<scalar>(),
                            scaleFactor.z()*rnd.GaussNormal<scalar>()
                        );
                }
                else
    #endif
                {
                    // Uniform distribution
                    newPoints[pointI] +=
                        minEdgeLength[pointI]
                       *vector
                        (
    #ifdef FOAMEXTEND
                            scaleFactor.x()*(2.0*rnd.scalar01() - 1.0),
                            scaleFactor.y()*(2.0*rnd.scalar01() - 1.0),
                            scaleFactor.z()*(2.0*rnd.scalar01() - 1.0)
    #else
                            scaleFactor.x()*(2.0*rnd.sample01<scalar>() - 1.0),
                            scaleFactor.y()*(2.0*rnd.sample01<scalar>() - 1.0),
                            scaleFactor.z()*(2.0*rnd.sample01<scalar>() - 1.0)
    #endif
                        );
                }
            }
        }

        // Remove the normal component on boundary patches
        forAll(mesh.boundary(), patchI)
        {
            const pointField& pointNormals =
                mesh.boundaryMesh()[patchI].pointNormals();
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(pointNormals, pI)
            {
                const vector& n = pointNormals[pI];
                const label pointID = meshPoints[pI];
                const vector disp = newPoints[pointID] - oldPoints[pointID];

                newPoints[pointID] = oldPoints[pointID] + ((I - sqr(n)) & disp);
            }
        }

        // Correct points for 2-D
        twoDPointCorrector twoD(mesh);
        twoD.correctPoints(newPoints);

        // Move the mesh
        Info<< "Applying the perturbation to the points" << endl;
        mesh.movePoints(newPoints);
#ifdef OPENFOAM_COM
    #if (OPENFOAM >= 2112)
        mesh.setPhi()->writeOpt() = IOobject::NO_WRITE;
    #else
        mesh.setPhi().writeOpt() = IOobject::NO_WRITE;
    #endif
#endif

        // Check for negative or small cell volumes
        const scalarField& VI = mesh.V();
        boolList negativeCellVol(VI.size(), false);
        int nNegCellVol = 0;
        int nSmallCellVol = 0;
        forAll(VI, cellI)
        {
            if (VI[cellI] < VSMALL)
            {
                negativeCellVol[cellI] = true;
                ++nNegCellVol;
            }
            else if (VI[cellI] < minV[cellI])
            {
                negativeCellVol[cellI] = true;
                ++nSmallCellVol;
            }
        }

        // Check if there are any severely non-orthogonal faces
#ifdef OPENFOAM_NOT_EXTEND
        const label nNonOrthoFaces = numSevereNonOrthoFaces(mesh);
#else
        const label nNonOrthoFaces = 0;
#endif

        // A valid mesh has no negative or small volumes and no severely
        // non-orthogonal faces
        validMesh = bool((nNegCellVol + nSmallCellVol + nNonOrthoFaces) == 0);

        if (validMesh)
        {
            Info<< "    There are no negative or small cell volumes or "
                << "severely non-orthogonal faces" << endl;
        }
        else
        {
            Info<< "    Number of cells with negative volumes: " << nNegCellVol
                << nl
                << "    Number of cells with small volumes: " << nSmallCellVol
                << nl
                << "    Number of severely non-orthogonal faces: "
                << nNonOrthoFaces
                << endl;

            // Expand the negativeCellVol marker field by layers of neighbours
            // This will more quickly fix bad cells, and may be required in some
            // case, e.g., a cell may be bad because neighbouring cells have
            // moved in a way that makes it bad
            // We will determine the number of layers based on the number of
            // iteration
            const label nFreq = 1;
            const label maxLayers = 100;
            const label nLayers = min(iter/nFreq + 1, maxLayers);
            Info<< "    Expanding the smoothing region by " << nLayers
                << " layers" << endl;
            for (label layerI = 0; layerI < nLayers; ++layerI)
            {
                // Take a copy of the negative volume field
                const boolList oldNegativeCellVol(negativeCellVol);
                forAll(oldNegativeCellVol, cI)
                {
                    if (oldNegativeCellVol[cI])
                    {
                        const labelList& curCellCells = cellCells[cI];
                        forAll(curCellCells, ccI)
                        {
                            const label neiCellID = curCellCells[ccI];
                            negativeCellVol[neiCellID] = true;
                        }
                    }
                }
            }

            // Reduce the minEdgeLength for all points, which are contained in
            // a negative volume cell
            boolList minEdgeLengthUpdated(minEdgeLength.size(), false);
            forAll(minEdgeLength, pI)
            {
                const labelList& curPointCells = pointCells[pI];
                bool negVol = false;
                forAll(curPointCells, pcI)
                {
                    const label cellID = curPointCells[pcI];

                    if (negativeCellVol[cellID])
                    {
                        negVol = true;
                        break;
                    }
                }

                if (negVol)
                {
                    // Update the minEdgeLength for this point
                    if (!minEdgeLengthUpdated[pI])
                    {
                        minEdgeLength[pI] *= beta;
                        minEdgeLengthUpdated[pI] = true;
                    }

                    // Update the minEdgeLength for the neighbouring points
                    // const labelList& curPointPoints = pointPoints[pI];
                    // forAll(curPointPoints, ppI)
                    // {
                    //     const label curPointID = curPointPoints[ppI];
                    //     if (!minEdgeLengthUpdated[curPointID])
                    //     {
                    //         minEdgeLength[curPointID] *= beta;
                    //         minEdgeLengthUpdated[curPointID] = true;
                    //     }

                    //     // Update the minEdgeLength for the second neighbours
                    //     // const labelList& secondPointPoints =
                    //     //     pointPoints[curPointID];
                    //     // forAll(secondPointPoints, sppI)
                    //     // {
                    //     //     const label secondPointID = secondPointPoints[sppI];
                    //     //     if (!minEdgeLengthUpdated[secondPointID])
                    //     //     {
                    //     //         minEdgeLength[secondPointID] *= beta;
                    //     //         minEdgeLengthUpdated[secondPointID] = true;
                    //     //     }
                    //     // }
                    // }
                }
            }

            if (iter == maxIter)
            {
                if (nRndReset++ == 10)
                {
                    FatalError
                        << "Maximum mesh correction and seed reset steps "
                        << "reached, but the mesh is still invalid."
                        << abort(FatalError);
                }
                else
                {
#ifdef OPENFOAM_COM
                    Warning
                        << "Maximum mesh correction steps reached, but the mesh "
                        << "is still invalid" << endl;

                    Info<< "Resetting the random number generator seed" << endl;

                    // Change the seed for the random number generator
                    rnd.reset(seed + 1);

                    // Reset minEdgeLength
                    minEdgeLength = oldMinEdgeLength;

                    // Reset iter
                    iter = 0;
#else
                    FatalError
                        << "Maximum mesh correction steps reached, but the mesh "
                        << "is still invalid" << abort(FatalError);
#endif
                }
            }
        }
    }
    while (!validMesh && iter++ < maxIter);

    // Write the mesh
    Info<< "Writing the mesh" << endl;
    mesh.setInstance(mesh.polyMesh::instance());
    mesh.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
