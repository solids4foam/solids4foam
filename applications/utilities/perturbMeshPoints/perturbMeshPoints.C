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

    This utility is useful for creating distorted grids for testing
    discretisations.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "Random.H"
#include "twoDPointCorrector.H"

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
    const scalar seed(readScalar(perturbDict.lookup("seed")));
    const vector scaleFactor(perturbDict.lookup("scaleFactor"));
    const wordList fixedPatchesList(perturbDict.lookup("fixedPatches"));
#ifdef OPENFOAM_COM
    const scalar minCos(readScalar(perturbDict.lookup("minCos")));
    const Switch Gaussian(perturbDict.lookup("Gaussian"));
#else
    const scalar minCos = 0;
#endif

    // Convert fixedPatches list to a set
    HashSet<word> fixedPatches;
    forAll(fixedPatchesList, pI)
    {
        fixedPatches.insert(fixedPatchesList[pI]);
    }

    // Create random number generator
    Random rnd(seed);

    // Calculate new points
    const pointField oldPoints = mesh.points();
    pointField newPoints(mesh.points());

    // Calculate the minimum edge length connected to each point
    const labelListList& pointEdges = mesh.pointEdges();
    const edgeList& edges = mesh.edges();
    scalarField minEdgeLength(newPoints.size(), GREAT);
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

    // Calculate a mask to identify fixed points
    boolList fixedPoint(newPoints.size(), false);
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
        const labelList& meshPoints = mesh.boundaryMesh()[patchI].meshPoints();

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
    mesh.movePoints(newPoints);

    // Write the mesh
    Info<< "Writing the mesh" << endl;
    mesh.setInstance(mesh.polyMesh::instance());
    mesh.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
