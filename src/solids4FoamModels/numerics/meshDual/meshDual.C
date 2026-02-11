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

\*---------------------------------------------------------------------------*/

#include "meshDual.H"
#include "meshDualiser.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "meshTools.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "polyTopoChange.H"
#else
    #include "directTopoChange.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshDual, 0);
}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::meshDual::makeDualMesh()
{
    if (dualMeshToMeshMapPtr_.valid())
    {
        FatalErrorIn("void Foam::solidModel::makeDualMesh() const")
            << "Pointers already set!" << abort(FatalError);
    }

    Info<< "Creating dualMesh" << endl;

    // Steps:
    // 1. Copy the current mesh: it is assumed that this already happen during
    //    the constructor, apart from the boundary
    // 2. Convert the mesh to its dual mesh
    // 3. Create the dual-mesh-to-mesh and mesh-to-dual-mesh maps

    // Copy boundary patches
    fvMesh& dualMesh = *this;
    List<polyPatch*> patches(mesh().boundaryMesh().size());
    forAll(patches, patchI)
    {
        // Clone the patch and reset the referenced boundary mesh
        patches[patchI] =
            mesh().boundaryMesh()[patchI].clone(dualMesh.boundaryMesh()).ptr();
    }

    // Add boundary patches
    dualMesh.addFvPatches(patches, true);


    // 2. Convert the mesh to its dual mesh

    // Create one dual face between cell centres and face centres
    const bool splitAllFaces = true;

    // Hard-coded settings
    const scalar featureAngle
    (
        dict_.lookupOrDefault<scalar>
        (
            word("featureAngle"), scalar(30)
        )
    );
    const bool doNotPreserveFaceZones
    (
        dict_.lookupOrDefault<bool>("doNotPreserveFaceZones", false)
    );
    const bool concaveMultiCells
    (
        dict_.lookupOrDefault<bool>("concaveMultiCells", false)
    );

    Info<< "    featureAngle: " << featureAngle << nl
        << "    doNotPreserveFaceZones: " << doNotPreserveFaceZones << nl
        << "    concaveMultiCells: " << concaveMultiCells << endl;

    // Mark boundary edges and points
    PackedBoolList isBoundaryEdge(dualMesh.nEdges());
    for
    (
        label faceI = dualMesh.nInternalFaces();
        faceI < dualMesh.nFaces();
        faceI++
    )
    {
        const labelList& fEdges = dualMesh.faceEdges()[faceI];

        forAll(fEdges, i)
        {
            isBoundaryEdge.set(fEdges[i], 1);
        }
    }

    // Face(centre)s that need inclusion in the dual mesh
    labelList featureFaces;
    // Edge(centre)s
    labelList featureEdges;
    // Points (that become a single cell) that need inclusion
    labelList singleCellFeaturePoints;
    // Points (that become a multiple cells)
    labelList multiCellFeaturePoints;

    // Sample implementation of feature detection.
    simpleMarkFeatures
    (
        dualMesh,
        isBoundaryEdge,
        featureAngle,
        concaveMultiCells,
        doNotPreserveFaceZones,
        featureFaces,
        featureEdges,
        singleCellFeaturePoints,
        multiCellFeaturePoints
    );

    // If we want to split all polyMesh faces into one dualface per cell
    // we are passing through we also need a point
    // at the polyMesh facecentre and edgemid of the faces we want to
    // split.
    if (splitAllFaces)
    {
        featureEdges = identity(dualMesh.nEdges());
        featureFaces = identity(dualMesh.nFaces());
    }

    // Topo change container
#ifdef OPENFOAM_NOT_EXTEND
    polyTopoChange meshMod(dualMesh.boundaryMesh().size());
#else
    directTopoChange meshMod(dualMesh.boundaryMesh().size());
#endif

    // Mesh dualiser engine
    meshDualiser dualMaker(dualMesh);

    // Insert all commands into directTopoChange to create dual of mesh.
    // This does all the hard work.
    dualMaker.setRefinement
    (
        splitAllFaces,
        featureFaces,
        featureEdges,
        singleCellFeaturePoints,
        multiCellFeaturePoints,
        meshMod
    );

    // Create mesh, return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(dualMesh, false);

    // Update fields
    dualMesh.updateMesh(map);

    if (map().hasMotionPoints())
    {
        dualMesh.movePoints(map().preMotionPoints());
    }

    // Find a way to add this default value to 'solidProperties' dict in the future
    if (dict_.lookupOrDefault<Switch>("writeDualMesh", false))
    {
        dualMesh.setInstance(runTime().constant());
        Info<< "Writing dualMesh to " << dualMesh.polyMesh::instance() << endl;
        dualMesh.write();

        // Write feature set fields
        const pointMesh& pMesh = pointMesh::New(mesh());
        pointScalarField cellFeaturePoints
        (
            IOobject
            (
                "cellFeaturePoints",
                runTime().timeName(),
                runTime(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", dimless, 0.0)
        );
        forAll(singleCellFeaturePoints, i)
        {
            cellFeaturePoints[singleCellFeaturePoints[i]] += 1;
        }
        forAll(multiCellFeaturePoints, i)
        {
            cellFeaturePoints[multiCellFeaturePoints[i]] += 2;
        }
        Info<< "Writing " << cellFeaturePoints.name() << endl;
        cellFeaturePoints.write();
    }


    // 3. Create the dual-mesh-to-mesh and mesh-to-dual-mesh maps
    dualMeshToMeshMapPtr_.set
    (
        new dualMeshToMeshMap(mesh(), dualMesh, dualMaker)
    );
}


// Copied from polyDualMeshApp.C
void Foam::meshDual::simpleMarkFeatures
(
    const polyMesh& mesh,
    const PackedBoolList& isBoundaryEdge,
    const scalar featureAngle,
    const bool concaveMultiCells,
    const bool doNotPreserveFaceZones,
    labelList& featureFaces,
    labelList& featureEdges,
    labelList& singleCellFeaturePoints,
    labelList& multiCellFeaturePoints
) const
{
#ifdef OPENFOAM_NOT_EXTEND
    scalar minCos = Foam::cos(featureAngle*constant::mathematical::pi/180.0);
#else
    scalar minCos = Foam::cos(featureAngle*mathematicalConstant::pi/180.0);
#endif

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Working sets
    labelHashSet featureEdgeSet;
    labelHashSet singleCellFeaturePointSet;
    labelHashSet multiCellFeaturePointSet;


    // 1. Mark all edges between patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& meshEdges = pp.meshEdges();

        // All patch corner edges. These need to be feature points & edges!
        for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
        {
            label meshEdgeI = meshEdges[edgeI];
            featureEdgeSet.insert(meshEdgeI);
            singleCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][0]);
            singleCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][1]);
        }
    }



    // 2. Mark all geometric feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Make distinction between convex features where the boundary point becomes
    // a single cell and concave features where the boundary point becomes
    // multiple 'half' cells.

    // Addressing for all outside faces
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces()-mesh.nInternalFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // Check for non-manifold points (surface pinched at point)
    allBoundary.checkPointManifold(false, &singleCellFeaturePointSet);

    // Check for non-manifold edges (surface pinched at edge)
    const labelListList& edgeFaces = allBoundary.edgeFaces();
    const labelList& meshPoints = allBoundary.meshPoints();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            const edge& e = allBoundary.edges()[edgeI];

            //Info<< "Detected non-manifold boundary edge:" << edgeI
            //    << " coords:"
            //    << allBoundary.points()[meshPoints[e[0]]]
            //    << allBoundary.points()[meshPoints[e[1]]] << endl;

            singleCellFeaturePointSet.insert(meshPoints[e[0]]);
            singleCellFeaturePointSet.insert(meshPoints[e[1]]);
        }
    }

    // Check for features.
    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 2)
        {
            label f0 = eFaces[0];
            label f1 = eFaces[1];

            // check angle
            const vector& n0 = allBoundary.faceNormals()[f0];
            const vector& n1 = allBoundary.faceNormals()[f1];

            if ((n0 & n1) < minCos)
            {
                const edge& e = allBoundary.edges()[edgeI];
                label v0 = meshPoints[e[0]];
                label v1 = meshPoints[e[1]];

                label meshEdgeI = meshTools::findEdge(mesh, v0, v1);
                featureEdgeSet.insert(meshEdgeI);

                // Check if convex or concave by looking at angle
                // between face centres and normal
                vector c1c0
                (
                    allBoundary[f1].centre(allBoundary.points())
                  - allBoundary[f0].centre(allBoundary.points())
                );

                if (concaveMultiCells && (c1c0 & n0) > SMALL)
                {
                    // Found concave edge. Make into multiCell features
                    Info<< "Detected concave feature edge:" << edgeI
                        << " cos:" << (c1c0 & n0)
                        << " coords:"
                        << allBoundary.points()[v0]
                        << allBoundary.points()[v1]
                        << endl;

                    singleCellFeaturePointSet.erase(v0);
                    multiCellFeaturePointSet.insert(v0);
                    singleCellFeaturePointSet.erase(v1);
                    multiCellFeaturePointSet.insert(v1);
                }
                else
                {
                    // Convex. singleCell feature.
                    if (!multiCellFeaturePointSet.found(v0))
                    {
                        singleCellFeaturePointSet.insert(v0);
                    }
                    if (!multiCellFeaturePointSet.found(v1))
                    {
                        singleCellFeaturePointSet.insert(v1);
                    }
                }
            }
        }
    }


    // 3. Mark all feature faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Face centres that need inclusion in the dual mesh
    labelHashSet featureFaceSet(mesh.nFaces()-mesh.nInternalFaces());
    // A. boundary faces.
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        featureFaceSet.insert(faceI);
    }

#ifdef OPENFOAM_ORG
    typedef meshFaceZones faceZoneMesh;
#endif

    // B. face zones.
    const faceZoneMesh& faceZones = mesh.faceZones();

    if (doNotPreserveFaceZones)
    {
        if (faceZones.size() > 0)
        {
            WarningIn("simpleMarkFeatures(..)")
                << "Detected " << faceZones.size()
                << " faceZones. These will not be preserved."
                << endl;
        }
    }
    else
    {
        if (faceZones.size() > 0)
        {
            Info<< "Detected " << faceZones.size()
                << " faceZones. Preserving these by marking their"
                << " points, edges and faces as features." << endl;
        }

        forAll(faceZones, zoneI)
        {
            const faceZone& fz = faceZones[zoneI];

            Info<< "Inserting all faces in faceZone " << fz.name()
                << " as features." << endl;

            forAll(fz, i)
            {
                label faceI = fz[i];
                const face& f = mesh.faces()[faceI];
                const labelList& fEdges = mesh.faceEdges()[faceI];

                featureFaceSet.insert(faceI);
                forAll(f, fp)
                {
                    // Mark point as multi cell point (since both sides of
                    // face should have different cells)
                    singleCellFeaturePointSet.erase(f[fp]);
                    multiCellFeaturePointSet.insert(f[fp]);

                    // Make sure there are points on the edges.
                    featureEdgeSet.insert(fEdges[fp]);
                }
            }
        }
    }

    // Transfer to arguments
    featureFaces = featureFaceSet.toc();
    featureEdges = featureEdgeSet.toc();
    singleCellFeaturePoints = singleCellFeaturePointSet.toc();
    multiCellFeaturePoints = multiCellFeaturePointSet.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshDual::meshDual(const fvMesh& mesh, const dictionary& dict)
:
    fvMesh
    (
        Foam::IOobject
        (
            "dualMesh",
            mesh.time().constant(),
            mesh.time(),
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE
        ),
#ifdef OPENFOAM_NOT_EXTEND
        pointField(mesh.points()),
        faceList(mesh.faces()),
        labelList(mesh.faceOwner()),
        labelList(mesh.faceNeighbour()),
#else
        xferCopy(mesh.points()),
        xferCopy(mesh.faces()),
        xferCopy(mesh.faceOwner()),
        xferCopy(mesh.faceNeighbour()),
#endif
        true
    ),
    mesh_(mesh),
    dict_(dict),
    dualMeshToMeshMapPtr_()
{
    // Create the dual mesh
    makeDualMesh();
}


// ************************************************************************* //
