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

#include "solidModel.H"
#include "volFields.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "solidTractionFvPatchVectorField.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "primitivePatchInterpolation.H"
    #include "polyTopoChange.H"
#else
    #include "blockSolidTractionFvPatchVectorField.H"
    #include "directTopoChange.H"
#endif
#include "fvcGradf.H"
#include "wedgePolyPatch.H"
#include "meshDualiser.H"
#include "meshTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidModel, 0);
    defineRunTimeSelectionTable(solidModel, dictionary);
    addToRunTimeSelectionTable(physicsModel, solidModel, physicsModel);

#ifdef OPENFOAMFOUNDATION
    typedef meshFaceZones faceZoneMesh;
#endif
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidModel::makeDualMesh() const
{
    if (dualMeshPtr_.valid() || dualMeshToMeshMapPtr_.valid())
    {
        FatalErrorIn("void Foam::solidModel::makeDualMesh() const")
            << "Pointers already set!" << abort(FatalError);
    }

    Info<< "Creating dualMesh" << endl;

    // Steps:
    // 1. Copy the current mesh
    // 2. Convert the mesh to its dual mesh
    // 3. Create the dual-mesh-to-mesh and mesh-to-dual-mesh maps

    // 1. Copy the current mesh
    dualMeshPtr_.set
    (
        new fvMesh
        (
            Foam::IOobject
            (
                "dualMesh",
                runTime().constant(),
                runTime(),
                Foam::IOobject::READ_IF_PRESENT
            ),
#ifdef OPENFOAMESIORFOUNDATION
            pointField(mesh().points()),
            faceList(mesh().faces()),
            labelList(mesh().faceOwner()),
            labelList(mesh().faceNeighbour()),
#else
            Xfer<pointField>(mesh().points()),
            Xfer<faceList>(mesh().faces()),
            Xfer<labelList>(mesh().faceOwner()),
            Xfer<labelList>(mesh().faceNeighbour()),
#endif
            true
        )
    );

    // Copy boundary patches
    fvMesh& dualMesh = dualMeshPtr_();
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
        solidModelDict().lookupOrDefault<scalar>("featureAngle", 30)
    );
    const bool doNotPreserveFaceZones
    (
        solidModelDict().lookupOrDefault<bool>("doNotPreserveFaceZones", false)
    );
    const bool concaveMultiCells
    (
        solidModelDict().lookupOrDefault<bool>("concaveMultiCells", false)
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
#ifdef OPENFOAMESIORFOUNDATION
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

    if (solidModelDict().lookupOrDefault<Switch>("writeDualMesh", false))
    {
        dualMesh.setInstance(runTime().constant());
        Info<< "Writing dualMesh to " << dualMesh.polyMesh::instance() << endl;
        dualMesh.write();

        // Write feature set fields
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
            pMesh_,
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
void Foam::solidModel::simpleMarkFeatures
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
#ifdef OPENFOAMESIORFOUNDATION
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


void Foam::solidModel::checkWedges() const
{
    const fvMesh& mesh = this->mesh();

    label nWedgePatches = 0;
    vector wedgeDirVec = vector::zero;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());

            // Make sure that solidWedge is used instead of wedge
            if
            (
                DD_.boundaryField()[patchI].type() == "wedge"
             && D_.boundaryField()[patchI].type() == "wedge"
            )
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "solidWedge should be used on displacement solution "
                    << "field wedge patches as non-orthogonal corrections "
                    << "are important!"
                    << abort(FatalError);
            }
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        if (nWedgePatches != 2)
        {
            FatalErrorIn("void Foam::solidModel::checkWedges() const")
                << "For axisymmetric cases, there should be exactly two wedge "
                << "patches!" << abort(FatalError);
        }

        Info<< nl << "Axisymmetric case: disabling the solution in the "
            << "out-of-plane direction" << endl;

        // We will use const_cast to disable the out-of-lane direction
        Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;

                wordList dirs(3);
                dirs[0] = "x";
                dirs[1] = "y";
                dirs[2] = "z";
                Info<< "    out-of-plane direction: " << dirs[cmpt] << nl
                    << endl;
            }
            else
            {
                solD[cmpt] = 1;
            }
        }
    }


    // Check all the face normals are in the same direction on the wedge patches
    // This is to avoid the case where a wedge patch is composed of two
    // disconnected regions with one on the front and one on the back
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // Unit face normals on processor
            const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

            if (nf.size() == 0)
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                    << nl << "Every processor should have at least one face on "
                    << "each wedge patch"
                    << abort(FatalError);
            }

            // Check that all the wedge face normals point in the same direction

            vector firstFaceNOnMasterProc = vector::zero;

            if (Pstream::master())
            {
                firstFaceNOnMasterProc = nf[0];
            }

            // Sync in parallel so that all processors have the master vector
            reduce(firstFaceNOnMasterProc, sumOp<vector>());

            forAll(nf, faceI)
            {
                if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
                {
                    FatalErrorIn("void Foam::solidModel::checkWedges() const")
                        << "On wedge patch "
                        << mesh.boundaryMesh()[patchI].name()
                        << " there are at "
                        << "least two faces with unit normals in the opposite "
                        << "directions" << nl
                        << "Please check that the wedge patches are correctly "
                        << "defined"
                        << abort(FatalError);
                }
            }
        }
    }
}


void Foam::solidModel::makeThermalModel() const
{
    if (!thermalPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeThermalModel() const")
            << "pointer already set!" << abort(FatalError);
    }

    thermalPtr_.set
    (
        new thermalModel(mesh())
    );
}


void Foam::solidModel::makeMechanicalModel() const
{
    if (!mechanicalPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeMechanicalModel() const")
            << "pointer already set!" << abort(FatalError);
    }

    mechanicalPtr_.set
    (
        new mechanicalModel(mesh(), nonLinGeom(), incremental())
    );
}


void Foam::solidModel::makeRho() const
{
    if (!rhoPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeRho() const")
            << "pointer already set!" << abort(FatalError);
    }

    rhoPtr_.set
    (
        new volScalarField(mechanical().rho())
    );
}


const Foam::pointVectorField& Foam::solidModel::pointDorPointDD() const
{
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Updated Lagrangian approaches move the mesh at the end of each
        // time-step so we use the increment of displacement field to calculate
        // the current deformed face zone points
        return pointDD();
    }
    else
    {
        // As linearGeometry and total Lagrangian approaches do not move the
        // mesh, we use the total displacement field to calculate the current
        // deformed face zone points
        return pointD();
    }
}


void Foam::solidModel::makeSetCellDisps() const
{
    if (setCellDispsPtr_.valid())
    {
        FatalErrorIn(type() + "::makeSetCellDisps() const")
            << "pointer already set!" << abort(FatalError);
    }

    if (solidModelDict().found("cellDisplacements"))
    {
        setCellDispsPtr_.set
        (
            new setCellDisplacements
            (
                mesh(), solidModelDict().subDict("cellDisplacements")
            )
        );
    }
    else
    {
        dictionary dict;
        setCellDispsPtr_.set(new setCellDisplacements(mesh(), dict));
    }
}


const Foam::setCellDisplacements& Foam::solidModel::setCellDisps() const
{
    if (setCellDispsPtr_.empty())
    {
        makeSetCellDisps();
    }

    return setCellDispsPtr_();
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::solidModel::dualMesh() const
{
    if (dualMeshPtr_.empty())
    {
        makeDualMesh();
    }

    return dualMeshPtr_();
}


Foam::fvMesh& Foam::solidModel::dualMesh()
{
    if (dualMeshPtr_.empty())
    {
        makeDualMesh();
    }

    return dualMeshPtr_();
}


const Foam::dualMeshToMeshMap& Foam::solidModel::dualMeshMap() const
{
    if (dualMeshToMeshMapPtr_.empty())
    {
        makeDualMesh();
    }

    return dualMeshToMeshMapPtr_();
}


Foam::thermalModel& Foam::solidModel::thermal()
{
    if (thermalPtr_.empty())
    {
        makeThermalModel();
    }

    return thermalPtr_();
}


Foam::mechanicalModel& Foam::solidModel::mechanical()
{
    if (mechanicalPtr_.empty())
    {
        makeMechanicalModel();
    }

    return mechanicalPtr_();
}


Foam::volScalarField& Foam::solidModel::rho()
{
    if (rhoPtr_.empty())
    {
        makeRho();
    }

    return rhoPtr_();
}


void Foam::solidModel::setCellDisps(fvVectorMatrix& DEqn)
{
    if (setCellDisps().cellIDs().size() == 0)
    {
        return;
    }

    if (incremental())
    {
        // Prepare the list of incremental displacements
        const vectorField& Dold = D().oldTime().internalField();
        const vectorField cellDisps = setCellDisps().cellDisps();
        vectorField cellIncrDisps(cellDisps.size(), vector::zero);
        const labelList cellIDs = setCellDisps().cellIDs();
        forAll(cellIncrDisps, cI)
        {
            cellIncrDisps[cI] = cellDisps[cI] - Dold[cellIDs[cI]];
        }

        DEqn.setValues(cellIDs, cellIncrDisps);
    }
    else
    {
        DEqn.setValues(setCellDisps().cellIDs(), setCellDisps().cellDisps());
    }
}


void Foam::solidModel::relaxField(volVectorField& D, int iCorr)
{
    if (relaxationMethod_ == "fixed")
    {
        // Fixed under-relaxation
        D.relax();
    }
    else if (relaxationMethod_ == "Aitken")
    {
        // See Aitken method at:
        // http://empire-multiphysics.com/projects/empire/wiki/Aitken_Relaxation
        // and
        // A partitioned solution approach for electro-thermo-
        // problems, Patrick Erbts, Stefan Hartmann, Alexander Duster.

        // Store aitkenResidual previous iteration
        aitkenResidual_.storePrevIter();

        // Calculate new aitkenResidual
        aitkenResidual_ = D.prevIter() - D;

        if (iCorr == 0)
        {
            // Fixed under-relaxation is applied on the first iteration
            aitkenAlpha_ = 1.0;

#ifdef OPENFOAMESIORFOUNDATION
            if (mesh().relaxField(D.name()))
            {
                aitkenAlpha_ =
                    mesh().fieldRelaxationFactor(D.name());
            }
#else
            if (mesh().solutionDict().relaxField(D.name()))
            {
                aitkenAlpha_ =
                    mesh().solutionDict().fieldRelaxationFactor(D.name());
            }
#endif
        }
        else
        {
            const volVectorField aitkenResidualDelta
            (
                aitkenResidual_.prevIter() - aitkenResidual_
            );

            // Update the relaxation factor field
            aitkenAlpha_ =
                aitkenAlpha_*(aitkenResidual_.prevIter() & aitkenResidualDelta)
               /(
                    magSqr(aitkenResidualDelta)
                  + dimensionedScalar("SMALL", dimLength*dimLength, SMALL)
                );

            // Bound alpha between 0.0 and 2.0
            // This may not be necessary but it seems to help convergence
            aitkenAlpha_ = max(0.0, min(2.0, aitkenAlpha_));
        }

        // Relax the field
        D -= aitkenAlpha_*aitkenResidual_;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::relaxField(volVectorField& D, int iCorr)"
        )   << "relaxationMethod '" << relaxationMethod_ << "' unknown!"
            << " Options are fixed or Aitken" << abort(FatalError);
    }
}


Foam::dictionary& Foam::solidModel::solidModelDict()
{
    return solidProperties_.subDict(type_ + "Coeffs");
}


void Foam::solidModel::makeRhoD2dt2D() const
{
    if (rhoD2dt2DPtr_.valid())
    {
        FatalErrorIn("void Foam::solidModel::makeRhoD2dt2D() const")
            << "Pointer already set" << abort(FatalError);
    }

    rhoD2dt2DPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "rhoD2dt2D",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector("zero", dimForce/dimVolume, vector::zero)
        )
    );
}


Foam::volVectorField& Foam::solidModel::rhoD2dt2D() const
{
    if (rhoD2dt2DPtr_.empty())
    {
        makeRhoD2dt2D();
    }

    return rhoD2dt2DPtr_();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime),
    regIOobject // ZT, Jul18: allow for multiple solid regions
    (
        IOobject
        (
            "solidModel_" + region,
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),
    dualMeshPtr_(),
    dualMeshToMeshMapPtr_(),
    solidProperties_
    (
        IOobject
        (
            "solidProperties",
            // If region == "region0" then read from the main case
            // Otherwise, read from the region/sub-mesh directory e.g.
            // constant/fluid or constant/solid
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    type_(type),
    thermalPtr_(),
    mechanicalPtr_(),
    Dheader_("D", runTime.timeName(), mesh(), IOobject::MUST_READ),
    DDheader_("DD", runTime.timeName(), mesh(), IOobject::MUST_READ),
    Uheader_("U", runTime.timeName(), mesh(), IOobject::MUST_READ),
    pointDheader_("pointD", runTime.timeName(), mesh(), IOobject::MUST_READ),
    D_
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    DD_
    (
        IOobject
        (
            "DD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    U_
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pMesh_(pointMesh::New(meshPtr_())),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradU_
    (
        IOobject
        (
            "grad(" + U_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless/dimTime, tensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    rhoPtr_(),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    dampingCoeff_
    (
        solidModelDict().lookupOrDefault<dimensionedScalar>
        (
            "dampingCoeff",
            dimensionedScalar("dampingCoeff", dimless/dimTime, 0)
        )
    ),
    stabilisationPtr_(),
    solutionTol_
    (
        solidModelDict().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
    ),
    alternativeTol_
    (
        solidModelDict().lookupOrDefault<scalar>("alternativeTolerance", 1e-07)
    ),
    materialTol_
    (
        solidModelDict().lookupOrDefault<scalar>("materialTolerance", 1e-05)
    ),
    infoFrequency_
    (
        solidModelDict().lookupOrDefault<int>("infoFrequency", 100)
    ),
    nCorr_(solidModelDict().lookupOrDefault<int>("nCorrectors", 10000)),
    minCorr_(solidModelDict().lookupOrDefault<int>("minCorrectors", 1)),
    maxIterReached_(0),
    residualFilePtr_(),
    writeResidualField_
    (
        solidModelDict().lookupOrDefault<Switch>("writeResidualField", false)
    ),
    enforceLinear_(false),
    relaxationMethod_
    (
        solidModelDict().lookupOrDefault<word>("relaxationMethod", "fixed")
    ),
    aitkenAlpha_
    (
        IOobject
        (
            "aitkenAlpha",
            runTime.constant(),
            meshPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    aitkenResidual_
    (
        IOobject
        (
            "aitkenResidual",
            runTime.constant(),
            meshPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    QuasiNewtonRestartFreq_
    (
        solidModelDict().lookupOrDefault<int>("QuasiNewtonRestartFrequency", 25)
    ),
    QuasiNewtonV_(QuasiNewtonRestartFreq_ + 2),
    QuasiNewtonW_(QuasiNewtonRestartFreq_ + 2),
    QuasiNewtonT_(QuasiNewtonRestartFreq_ + 2),
    DRef_
    (
        IOobject
        (
            "DRef",
            runTime.constant(),
            meshPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    unrelaxedDRef_
    (
        IOobject
        (
            "unrelaxedDRef",
            runTime.constant(),
            meshPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    globalPatchesPtrList_(),
    setCellDispsPtr_(),
    restart_
    (
        solidModelDict().lookupOrDefault<Switch>("restart", false)
    ),
    rhoD2dt2DPtr_()
{
    // Force old time fields to be stored
    D_.oldTime().oldTime();
    DD_.oldTime().oldTime();
    pointD_.oldTime();
    pointDD_.oldTime();
    gradD_.oldTime();
    gradDD_.oldTime();
    gradU_.oldTime();
    sigma_.oldTime();

    if (restart_)
    {
        // Enable writing of fields which are needed for restart
        D_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        D_.oldTime().oldTime().writeOpt() = IOobject::AUTO_WRITE;
        DD_.writeOpt() = IOobject::AUTO_WRITE;
        DD_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        DD_.oldTime().oldTime().writeOpt() = IOobject::AUTO_WRITE;
        U_.writeOpt() = IOobject::AUTO_WRITE;
        U_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        pointD_.writeOpt() = IOobject::AUTO_WRITE;
        pointD_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        pointDD_.writeOpt() = IOobject::AUTO_WRITE;
        gradD_.writeOpt() = IOobject::AUTO_WRITE;
        gradD_.oldTime().writeOpt() = IOobject::AUTO_WRITE;
        gradDD_.writeOpt() = IOobject::AUTO_WRITE;
        gradU_.writeOpt() = IOobject::AUTO_WRITE;
    }
    else
    {
        D_.oldTime().writeOpt() = IOobject::NO_WRITE;
        D_.oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
        DD_.writeOpt() = IOobject::NO_WRITE;
        DD_.oldTime().writeOpt() = IOobject::NO_WRITE;
        DD_.oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
        U_.oldTime().writeOpt() = IOobject::NO_WRITE;
        U_.oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
        pointD_.writeOpt() = IOobject::AUTO_WRITE;
        pointD_.oldTime().writeOpt() = IOobject::NO_WRITE;
        pointDD_.writeOpt() = IOobject::NO_WRITE;
        gradD_.writeOpt() = IOobject::NO_WRITE;
        gradD_.oldTime().writeOpt() = IOobject::NO_WRITE;
        gradDD_.writeOpt() = IOobject::NO_WRITE;
        gradU_.writeOpt() = IOobject::NO_WRITE;
    }

    // Print out the relaxation factor
    Info<< "    under-relaxation method: " << relaxationMethod_ << endl;
    if (relaxationMethod_ == "QuasiNewton")
    {
        Info<< "        restart frequency: " << QuasiNewtonRestartFreq_ << endl;
    }

    // If requested, create the residual file
    if (solidModelDict().lookupOrDefault<Switch>("residualFile", false))
    {
        if (Pstream::master())
        {
            Info<< "Creating residual.dat" << endl;
            residualFilePtr_.set
            (
                new OFstream(runTime.path()/"residual.dat")
            );
        }
    }

    // Create stabilisation object

    if (!solidModelDict().found("stabilisation"))
    {
        // If the stabilisation sub-dict is not found, we will add it with
        // default settings
        dictionary stabDict;
        stabDict.add("type", "RhieChow");
        stabDict.add("scaleFactor", 0.1);
        solidModelDict().add("stabilisation", stabDict);
    }

    stabilisationPtr_.set
    (
        new momentumStabilisation
        (
            solidModelDict().subDict("stabilisation")
        )
    );

    // If the case is axisymmetric, we will disable solving in the out-of-plane
    // direction
    // PC, 12-Nov-18: disabling the 3rd direction slows down convergence a lot
    // in some elastic cases: disabled for now
    //checkWedges();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{
    thermalPtr_.clear();
    mechanicalPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solidModel::rho() const
{
    if (rhoPtr_.empty())
    {
        makeRho();
    }

    return rhoPtr_();
}

const Foam::thermalModel& Foam::solidModel::thermal() const
{
    if (thermalPtr_.empty())
    {
        makeThermalModel();
    }

    return thermalPtr_();
}


const Foam::mechanicalModel& Foam::solidModel::mechanical() const
{
    if (mechanicalPtr_.empty())
    {
        makeMechanicalModel();
    }

    return mechanicalPtr_();
}


void Foam::solidModel::DisRequired()
{
#ifdef OPENFOAMESIORFOUNDATION
    if (!Dheader_.typeHeaderOk<volVectorField>(true))
#else
    if (!Dheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::DisRequired()")
            << "This solidModel requires the 'D' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::DDisRequired()
{
#ifdef OPENFOAMESIORFOUNDATION
    if (!DDheader_.typeHeaderOk<volVectorField>(true))
#else
    if (!DDheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::DDisRequired()")
            << "This solidModel requires the 'DD' field to be specified!"
            << abort(FatalError);
    }
}

void Foam::solidModel::UisRequired()
{
#ifdef OPENFOAMESIORFOUNDATION
    if (!Uheader_.typeHeaderOk<volVectorField>(true))
#else
    if (!Uheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::UisRequired()")
            << "This solidModel requires the 'U' field to be specified!"
            << abort(FatalError);
    }
}

void Foam::solidModel::pointDisRequired()
{
#ifdef OPENFOAMESIORFOUNDATION
    if (!pointDheader_.typeHeaderOk<pointVectorField>(true))
#else
    if (!pointDheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::pointDisRequired()")
            << "This solidModel requires the 'pointD' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::solidModel::makeGlobalPatches
(
    const wordList& patchNames,
    const bool currentConfiguration
) const
{
    globalPatchesPtrList_.setSize(patchNames.size());

    forAll(patchNames, i)
    {
        if (globalPatchesPtrList_.set(i))
        {
            FatalErrorIn
            (
                type() + "::makeGlobalPatches(const wordList&) const"
            )   << "Pointer already set for global patch: "
                << patchNames[i] << "!"
                << abort(FatalError);
        }

        if (currentConfiguration)
        {
            // The global patch will create a standAlone zone based on the
            // current point positions. So we will temporarily move the mesh to
            // the deformed position, then create the globalPatch, then move the
            // mesh back
            const pointField pointsBackup = mesh().points();

            // Lookup patch index
            const label patchID =
                mesh().boundaryMesh().findPatchID(patchNames[i]);
            if (patchID == -1)
            {
                FatalErrorIn("void Foam::solidModel::makeGlobalPatches(...)")
                    << "Patch not found!" << abort(FatalError);
            }

            // Patch point displacement
            const vectorField pointDisplacement
            (
                pointDorPointDD().internalField(),
                mesh().boundaryMesh()[patchID].meshPoints()
            );

            // Calculate deformation point positions
            const pointField newPoints
            (
                mesh().points() + pointDorPointDD().internalField()
            );

            // Move the mesh to deformed position
            // const_cast is justified as it is not our intention to permanently
            // move the mesh; however, it would be better if we did not need it
            mesh().V();
            const_cast<dynamicFvMesh&>(mesh()).movePoints(newPoints);
            const_cast<dynamicFvMesh&>(mesh()).moving(false);
#ifndef OPENFOAMESIORFOUNDATION
            const_cast<dynamicFvMesh&>(mesh()).changing(false);
#endif
            const_cast<dynamicFvMesh&>(mesh()).setPhi().writeOpt() =
                IOobject::NO_WRITE;

            // Create global patch based on deformed mesh
            globalPatchesPtrList_.set
            (
                i,
                new globalPolyPatch(patchNames[i], mesh())
            );

            // Force creation of standAlonePatch
            globalPatchesPtrList_[i].globalPatch();

            // Move the mesh back
            const_cast<dynamicFvMesh&>(mesh()).movePoints(pointsBackup);
            mesh().V();
            const_cast<dynamicFvMesh&>(mesh()).moving(false);
#ifndef OPENFOAMESIORFOUNDATION
            const_cast<dynamicFvMesh&>(mesh()).changing(false);
#endif
            const_cast<dynamicFvMesh&>(mesh()).setPhi().writeOpt() =
                IOobject::NO_WRITE;
        }
        else
        {
            globalPatchesPtrList_.set
            (
                i,
                new globalPolyPatch(patchNames[i], mesh())
            );
        }
    }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidModel::globalPatches() const
{
    if (globalPatchesPtrList_.empty())
    {
        FatalErrorIn(type() + "::globalPatches() const")
            << "makeGlobalPatches(const wordList&) must be called "
            << "before globalPatch can be called!"
            << abort(FatalError);
    }

    return globalPatchesPtrList_;
}


void Foam::solidModel::clearGlobalPatches() const
{
    globalPatchesPtrList_.clear();
}


Foam::vector Foam::solidModel::pointU(const label pointID) const
{
    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    mechanical().volToPoint().interpolate(U(), pointU);

    return pointU.internalField()[pointID];
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementIncrement
(
    const label interfaceI
) const
{
    // Create patch point field
    const vectorField patchPointDispIncr
    (
        pointDD().internalField(),
        globalPatches()[interfaceI].patch().meshPoints()
    );

    // Return the global patch field
    return globalPatches()[interfaceI].patchPointToGlobal(patchPointDispIncr);
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementOld
(
    const label interfaceI
) const
{
    // Create patch point field
    const vectorField patchPointDispOld
    (
        pointD().oldTime().internalField(),
        globalPatches()[interfaceI].patch().meshPoints()
    );

    // Return the global patch field
    return globalPatches()[interfaceI].patchPointToGlobal(patchPointDispOld);
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZoneAcceleration
(
    const label interfaceI
) const
{
    const volVectorField a(fvc::d2dt2(D()));

    return globalPatches()[interfaceI].patchFaceToGlobal
    (
        a.boundaryField()[globalPatches()[interfaceI].patch().index()]
    );
}


void Foam::solidModel::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void Foam::solidModel::end()
{
    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum momentum correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The momentum equation converged in all time-steps"
            << nl << endl;
    }
}


Foam::autoPtr<Foam::solidModel> Foam::solidModel::New
(
    Time& runTime,
    const word& region
)
{
    word solidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid model is created, otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                bool(region == dynamicFvMesh::defaultRegion)
              ? fileName(runTime.caseConstant())
              : fileName(runTime.caseConstant()/region),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidModel")
            >> solidModelTypeName;
    }

    Info<< "Selecting solidModel " << solidModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::New(Time&, const word&)"
        )   << "Unknown solidModel type " << solidModelTypeName
            << endl << endl
            << "Valid solidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidModel>(cstrIter()(runTime, region));
}


void Foam::solidModel::setTraction
(
    fvPatchVectorField& tractionPatch,
    const vectorField& traction
)
{
    if (tractionPatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(tractionPatch);

        patchD.traction() = traction;
    }
#ifndef OPENFOAMESIORFOUNDATION
    else if
    (
        tractionPatch.type() == blockSolidTractionFvPatchVectorField::typeName
    )
    {
        blockSolidTractionFvPatchVectorField& patchD =
            refCast<blockSolidTractionFvPatchVectorField>(tractionPatch);

        patchD.traction() = traction;
    }
#endif
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << tractionPatch.type()
            << " for patch " << tractionPatch.patch().name()
            << " should instead be type "
            << solidTractionFvPatchVectorField::typeName
#ifndef OPENFOAMESIORFOUNDATION
            << " or "
            << blockSolidTractionFvPatchVectorField::typeName
#endif
            << abort(FatalError);
    }
}


void Foam::solidModel::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    const vectorField patchTraction
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTraction)
    );

#ifdef OPENFOAMESIORFOUNDATION
    setTraction(solutionD().boundaryFieldRef()[patchID], patchTraction);
#else
    setTraction(solutionD().boundaryField()[patchID], patchTraction);
#endif
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const volScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const surfaceScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


void Foam::solidModel::writeFields(const Time& runTime)
{
    // Write strain fields
    // Currently only defined for linear geometry
    if (nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
    {
        // Total strain
        volSymmTensorField epsilon("epsilon", symm(gradD()));
        epsilon.write();

        // Equivalent strain
        volScalarField epsilonEq
        (
            "epsilonEq", sqrt((2.0/3.0)*magSqr(dev(epsilon)))
        );
        epsilonEq.write();

        Info<< "Max epsilonEq = " << gMax(epsilonEq) << endl;
    }

    // Calculate equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        "sigmaEq", sqrt((3.0/2.0)*magSqr(dev(sigma())))
    );
    sigmaEq.write();

    Info<< "Max sigmaEq (von Mises stress) = " << gMax(sigmaEq) << endl;

    // If asked, write the residual field
    if (writeResidualField_)
    {
        const volVectorField& D = solutionD();
#ifdef OPENFOAMESIORFOUNDATION
        scalar denom =
            gMax(mag(D.primitiveField() - D.oldTime().primitiveField()));
        if (denom < SMALL)
        {
            denom = max(gMax(mag(D.primitiveField())), SMALL);
        }
#else
        scalar denom =
            gMax(mag(D.internalField() - D.oldTime().internalField()));
        if (denom < SMALL)
        {
            denom = max(gMax(mag(D.internalField())), SMALL);
        }
#endif

        const volVectorField residualD
        (
            "residualD",
            (D - D.prevIter())/denom
        );

        Info<< "Writing residualD field" << endl;
        residualD.write();
    }

    physicsModel::writeFields(runTime);
}


Foam::scalar Foam::solidModel::newDeltaT()
{
    return min
    (
        runTime().deltaTValue(),
        mechanical().newDeltaT()
    );
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

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD, pointDD);

    // Fix, AW/PC, 22-Dec-20,
    // correctBoundaryConditions should not be called as it causes (global?)
    // points to become out of sync. This results in the error "face area does
    // not match neighbour..."
    //pointDD.correctBoundaryConditions();

#ifdef OPENFOAMESIORFOUNDATION
    vectorField& pointDDI = pointDD.primitiveFieldRef();
#else
    vectorField& pointDDI = pointDD.internalField();
#endif

    vectorField newPoints = oldPoints;

    // Correct symmetryPlane points

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            if
            (
                returnReduce(mesh().boundaryMesh()[patchI].size(), sumOp<int>())
            )
            {
                continue;
            }

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

            if
            (
                returnReduce(mesh().boundaryMesh()[patchI].size(), sumOp<int>())
            )
            {
                continue;
            }

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
    forAll(pointDDI, pointI)
    {
        newPoints[pointI] += pointDDI[pointI];
    }

    // Move unused globalFaceZone points
    // Not need anymore as globalFaceZones are not used
    //updateGlobalFaceZoneNewPoints(pointDDI, newPoints);

    twoDPointCorrector twoDCorrector(mesh());
    twoDCorrector.correctPoints(newPoints);
    twoDCorrector.correctPoints(pointDDI);
    mesh().movePoints(newPoints);
    mesh().V00();
    mesh().moving(false);
#ifndef OPENFOAMESIORFOUNDATION
    mesh().changing(false);
#endif
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

#ifndef OPENFOAMESIORFOUNDATION
    // Tell the mechanical model to move the subMeshes, if they exist
    mechanical().moveSubMeshes();
#endif
}


const Foam::dictionary& Foam::solidModel::solidModelDict() const
{
    return solidProperties_.subDict(type_ + "Coeffs");
}

// ************************************************************************* //
