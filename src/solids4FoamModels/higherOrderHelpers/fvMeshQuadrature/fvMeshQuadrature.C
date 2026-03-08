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

#include "fvMeshQuadrature.H"
#include "triangle.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "triFace.H"
#include "triQuadrature.H"
#include "tetQuadrature.H"
#include "lineQuadrature.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fvMeshQuadrature, 0);

const scalar fvMeshQuadrature::minTriAreaRatio_
(
    debug::floatOptimisationSwitch("fvMeshQuadratureMinTriAreaRatio", 0.05)
);

const scalar fvMeshQuadrature::minTriQuality_
(
    debug::floatOptimisationSwitch("fvMeshQuadratureMinTriQuality", 0.01)
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool fvMeshQuadrature::poorTriangulation
(
    const faceList& triFaces,
    const pointField& pts
) const
{
    scalar totalArea = 0.0;

    // Compute total area and check area of each triangle
    forAll(triFaces, triI)
    {
        const face& triF = triFaces[triI];

        const triPoints tri
        (
            pts[triF[0]],
            pts[triF[1]],
            pts[triF[2]]
        );
        const scalar triArea = tri.mag();
        totalArea += triArea;

        if (triArea < SMALL)
        {
            return true;
        }
    }

    if (totalArea < SMALL)
    {
        return true;
    }

    // Check individual triangles
    forAll(triFaces, triI)
    {
        const face& triF = triFaces[triI];

        const triPoints tri
        (
            pts[triF[0]],
            pts[triF[1]],
            pts[triF[2]]
        );

        const scalar triArea = tri.mag();
        const scalar areaRatio = triArea/(totalArea + VSMALL);

        if (areaRatio < minTriAreaRatio_)
        {
            return true;
        }

        if (triPointRef(tri).quality() < minTriQuality_)
        {
            return true;
        }
    }

    return false;
}


void fvMeshQuadrature::calcQuadPointsAndWeights() const
{
    if (mesh_.nGeometricD() == 2)
    {
        calcQuadPointsAndWeights2D();
    }
    else if (mesh_.nGeometricD() == 3)
    {
        calcQuadPointsAndWeights3D();
    }
    else
    {
        FatalErrorIn("calcQuadPointsAndWeights()")
            << "Only implemented for 2D and 3D models, current model is "
            << mesh_.nGeometricD() << "D model!"
            << abort(FatalError);
    }
}


void fvMeshQuadrature::calcQuadPointsAndWeights2D() const
{
    if
    (
        faceQuadPointsPtr_.valid() || faceQuadWeightsPtr_.valid()
     || cellQuadPointsPtr_.valid() || cellQuadWeightsPtr_.valid()
    )
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();

    // Empty direction vector
    vector emptyDir(vector::zero);

    label emptyCmpt = -1;
    const Vector<label>& solD = mesh_.solutionD();
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (solD[cmpt] == -1)
        {
            emptyDir[cmpt] = 1.0;
            break;
        }
    }

    if (emptyCmpt == -1)
    {
        FatalErrorInFunction
            << "Mesh is not 2D! No empty direction detected!"
            << abort(FatalError);
    }

    emptyDir[emptyCmpt] = 1.0;

    // Domain midpoint in empty direction
    const boundBox bb(pts);
    const scalar thickness = bb.span()[emptyCmpt];
    const scalar mid = bb.centre()[emptyCmpt];

    // Determine number of quadrature points per face, for storage allocation
    labelList quadPtsPerFace(mesh.nFaces(), 0);
    const label nLineQP = lineQuadrature::nPoints(faceOrder_);

    forAll(faces, faceI)
    {
        if (faceI >= mesh.nInternalFaces())
        {
            const label patchID = mesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh.boundaryMesh()[patchID];

            if (isA<emptyPolyPatch>(pp))
            {
                quadPtsPerFace[faceI] = 0;
                continue;
            }

            if (isA<wedgePolyPatch>(pp))
            {
                FatalErrorInFunction
                    << "Not implemented for axisymmetric case!"
                    << abort(FatalError);
            }
        }

        quadPtsPerFace[faceI] = nLineQP;
    }

    // Initialise quadrature points and weights for faces
    faceQuadPointsPtr_.set(new CompactListList<point>(quadPtsPerFace));
    faceQuadWeightsPtr_.set(new CompactListList<scalar>(quadPtsPerFace));

    CompactListList<point>&  faceQP  = *faceQuadPointsPtr_;
    CompactListList<scalar>& faceQW = *faceQuadWeightsPtr_;

    forAll(faces, faceI)
    {
        if (!faceQP[faceI].size())
        {
            // Skip empty faces
            continue;
        }

        // Set face quadrature points and corresponding weights. We will loop
        // over face edges and take the first edge on the empty patch.

        const face& curFace = faces[faceI];
        const edgeList curFaceEdges = curFace.edges();

        bool found = false;

        forAll(curFaceEdges, edgeI)
        {
            const edge& curEdge = curFaceEdges[edgeI];
            const vector e  = curEdge.vec(pts);
            const vector eNorm = e / mag(e);

            // In-plane edge, perpendicular to empty direction
            if (mag(eNorm & emptyDir) < SMALL)
            {
                // This edge is perpendicular to empty direction, we will use it
                point a = pts[curEdge.start()];
                point b = pts[curEdge.end()];

               // Translate edge to mid-plane along empty direction
                a[emptyCmpt] = mid;
                b[emptyCmpt] = mid;

                // Construct line and lineQuadrature
                const linePointRef l(a,b);
                const lineQuadrature lq(l, faceOrder_);

                const List<point>& points = lq.points();
                const List<scalar>& weights = lq.weights();

                forAll(points, pI)
                {
                    faceQP[faceI][pI] = points[pI];
                    faceQW[faceI][pI] = weights[pI] * l.mag() * thickness;
                }

                // Go to next face, this face is done
                found = true;
                break;
            }
        }

        if (!found)
        {
            FatalErrorInFunction
                << "Could not find an in-plane edge for face " << faceI
                << ". Check 2D mesh, face topology and used tolerance"
                << abort(FatalError);
        }
    }

    // Cell quadrature points and weights
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const label nTriQP = triQuadrature::nPoints(cellOrder_);

    labelList selectedFace(mesh.nCells(), -1);
    labelList quadPtsPerCell(mesh.nCells(), 0);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        const vectorField& patchSf = pp.faceAreas();
        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;
            const vector& Sf = patchSf[patchFaceI];

            if ((Sf & emptyDir) > SMALL)
            {
                const label cellI = mesh.faceOwner()[faceI];

                if (selectedFace[cellI] != -1)
                {
                    FatalErrorInFunction
                        << "More than one positively-oriented empty face found "
                        << "for cell " << cellI
                        << abort(FatalError);
                }

                selectedFace[cellI] = faceI;
                quadPtsPerCell[cellI] = faces[faceI].nTriangles()*nTriQP;
            }
        }
    }

    forAll(selectedFace, cellI)
    {
        if (selectedFace[cellI] == -1)
        {
            FatalErrorInFunction
                << "No positively-oriented empty face found for cell "
                << cellI << abort(FatalError);
        }
    }

    // Allocate memory for compactListList for cell points and weights
    cellQuadPointsPtr_.set(new CompactListList<point>(quadPtsPerCell));
    cellQuadWeightsPtr_.set(new CompactListList<scalar>(quadPtsPerCell));

    CompactListList<point>&  cellQP = *cellQuadPointsPtr_;
    CompactListList<scalar>& cellQW = *cellQuadWeightsPtr_;

    // Triangulate selected faces and fill cell quadrature
    forAll(selectedFace, cellI)
    {
        const label faceI = selectedFace[cellI];
        const face& f = faces[faceI];

        const label nTri = f.nTriangles();
        faceList triFaces(nTri);

        label t2 = 0;
        const label t1 = f.triangles(pts, t2, triFaces);

        if (nTri != t1 || nTri != t2)
        {
            FatalErrorInFunction
                << "The numbers of reported triangles in the face do not "
                << "match that generated by the triangulation for face "
                << faceI << abort(FatalError);
        }

        if (poorTriangulation(triFaces, pts))
        {
            WarningInFunction
                << "Cell " << cellI  << " has poor triangulation."
                << "Consider improving the mesh or implement central point"
                << " decomposition." << endl;
        }

        // Now fill quadrature points and weights
        forAll(triFaces, triI)
        {
            point a = pts[triFaces[triI][0]];
            point b = pts[triFaces[triI][1]];
            point c = pts[triFaces[triI][2]];

            a[emptyCmpt] = mid;
            b[emptyCmpt] = mid;
            c[emptyCmpt] = mid;

            const triPoints tp(a, b, c);
            const scalar triArea = tp.mag();


            const triQuadrature tq(tp, cellOrder_);
            const List<point>& points = tq.points();
            const List<scalar>& weights = tq.weights();

            const label base = triI*nTriQP;

            forAll(points, i)
            {
                cellQP[cellI][base + i] = points[i];
                cellQW[cellI][base + i] = weights[i] * triArea * thickness;
            }
        }
    }
}


void fvMeshQuadrature::calcQuadPointsAndWeights3D() const
{
    if
    (
     faceQuadPointsPtr_.valid() || faceQuadWeightsPtr_.valid()
     || cellQuadPointsPtr_.valid() || cellQuadWeightsPtr_.valid()
    )
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    // Preliminaries
    const fvMesh& mesh = mesh_;
    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();
    const List<cell>& cells = mesh.cells();
    const vectorField& C = mesh_.C().internalField();


    // Decompose faces into triangles.
    // We have two options (N-number of face points):
    //                      - central point triangulation (N triangles)
    //                      - fan triangulation (N-2 triangles)
    // We will use fan triangulation for face decomposition. In the case that
    // resulting faces are small or invalid we will switch to more robust
    // central point decomposition (for each face separately).
    // Instead of central decomposition, the alternative is interior constrained
    // Delaunay triangulation which is much more complex to implement.
    // Note: We are using fan triangulation from face.H class which is smart
    // and adaptive triangulation.

    // Points of each triangle sub-element
    List<List<triPoints>> faceTri(mesh.nFaces());

    labelList quadPtsPerFace(mesh.nFaces(), 0);
    bool fallbackActivated = false;
    const label nTriQP = triQuadrature::nPoints(faceOrder_);

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];
        const label nTri = f.nTriangles();
        const label nPoints = f.size();

        // Fan triangulation
        faceList triFaces(nTri);
        label t2 = 0;
        const label t1 = f.triangles(pts, t2, triFaces);

        if (nTri != t1 || nTri != t2)
        {
            FatalErrorInFunction
                << "The numbers of reported triangles in the face do not "
                << "match that generated by the triangulation for face "
                << faceI << abort(FatalError);
        }

        faceTri[faceI].setSize(nTri);

        forAll(triFaces, triI)
        {
            const face& triF = triFaces[triI];
            const triPoints tp(pts[triF[0]], pts[triF[1]], pts[triF[2]]);
            faceTri[faceI][triI] = tp;
        }

        // Check decomposition and perform central point triangulation if needed
        if (allowDegenerateTriFallback_ && poorTriangulation(triFaces, pts))
        {
            // Triangulation using central point
            fallbackActivated = true;

            faceTri[faceI].setSize(nPoints);
            const point fc = f.centre(pts);

            for (label pI = 0; pI<nPoints; ++pI)
            {
                const label nextpI = (pI + 1 < nPoints ? pI + 1 : 0);

                const triPoints tp(pts[f[pI]], pts[f[nextpI]], fc);
                faceTri[faceI][pI] = tp;
            }
        }

        // Final number of quadrature points is:
        // number of face triangles * quadrature points per triangle
        quadPtsPerFace[faceI] = faceTri[faceI].size() * nTriQP;
    }

    if (fallbackActivated)
    {
        WarningInFunction
            << "Switching to central point triangulation for one or more faces!"
            << endl;
    }

    // Allocate face quadrature storage and fill
    faceQuadPointsPtr_.set(new CompactListList<point>(quadPtsPerFace));
    faceQuadWeightsPtr_.set(new CompactListList<scalar>(quadPtsPerFace));

    CompactListList<point>&  faceQP = *faceQuadPointsPtr_;
    CompactListList<scalar>& faceQW = *faceQuadWeightsPtr_;

    forAll(faceTri, faceI)
    {
        const List<triPoints>& fT = faceTri[faceI];

        forAll(fT, tI)
        {
            const triPoints& tp = fT[tI];
            const scalar triArea = tp.mag();

            const triQuadrature tq(tp, faceOrder_);
            const List<point>& points = tq.points();
            const List<scalar>& weights = tq.weights();

            forAll(points, i)
            {
                const label pos = tI*nTriQP + i;
                faceQP[faceI][pos] = points[i];
                faceQW[faceI][pos] = weights[i] * triArea;
            }
        }
    }

    // Decompose cells into tetrahedrals.
    // We will use existing face triangulation and cell centre to construct
    // tetrahedrals.
    labelList quadPtsPerCell(mesh.nCells(), 0);

    const label nTetQP = tetQuadrature::nPoints(cellOrder_);

    forAll(quadPtsPerCell, cellI)
    {
        const cellShape& shape = mesh.cellShapes()[cellI];

        if (shape.model() == cellModel::ref(cellModel::TET))
        {
            quadPtsPerCell[cellI] = nTetQP;
        }
        else
        {
            const cell& c = cells[cellI];
            label cellTets = 0;

            forAll(c, fI)
            {
                const label faceI = c[fI];
                cellTets += faceTri[faceI].size();
            }
            quadPtsPerCell[cellI] = nTetQP * cellTets;
        }
    }

    // Allocate memory for compactListList for points and weights
    cellQuadPointsPtr_.set(new CompactListList<point>(quadPtsPerCell));
    cellQuadWeightsPtr_.set(new CompactListList<scalar>(quadPtsPerCell));

    CompactListList<point>&  cellQP = *cellQuadPointsPtr_;
    CompactListList<scalar>& cellQW = *cellQuadWeightsPtr_;

    // Loop over cells

    forAll (cells, cellI)
    {
        const cellShape& shape = mesh.cellShapes()[cellI];
        const labelList& cellPoints = mesh.cellPoints()[cellI];
        const point& cellC = C[cellI];
        const cell& c = mesh.cells()[cellI];

        // Skip decomposition for tetrahedral cells
        if (shape.model() == cellModel::ref(cellModel::TET))
        {
            const tetPoints tet =
                tetPoints
                (
                    pts[cellPoints[0]],
                    pts[cellPoints[1]],
                    pts[cellPoints[2]],
                    pts[cellPoints[3]]
                 );

            // Get tet quadrature points and their weight
            const tetQuadrature tq(tet, cellOrder_);
            const List<point>& points = tq.points();
            const List<scalar>& weights = tq.weights();

            tetPointRef tetRef(tet);
            const scalar tetV = mag(tetRef.mag());

            // Store quadrature points and weights
            forAll(points, i)
            {
                cellQP[cellI][i] = points[i];
                cellQW[cellI][i] = weights[i] * tetV;
            }
        }
        else
        {
            label tetI = 0;

            // Loop over faces of the cell
            forAll(c, fI)
            {
                const label faceI = c[fI];
                const List<triPoints>& tris = faceTri[faceI];

                forAll(tris, triI)
                {
                    const triPoints& tp = tris[triI];
                    const tetPoints subTet(cellC, tp.a(), tp.b(), tp.c());
                    const scalar tetV = mag(tetPointRef(subTet).mag());

                    const tetQuadrature tq(subTet, cellOrder_);
                    const List<point>& points = tq.points();
                    const List<scalar>& weights = tq.weights();

                    // Store quadrature points and weights
                    const label base = tetI*nTetQP;
                    forAll(points, i)
                    {
                        cellQP[cellI][base + i] = points[i];
                        cellQW[cellI][base + i] = weights[i]*tetV;
                    }
                    ++tetI;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvMeshQuadrature::fvMeshQuadrature
(
    const fvMesh& mesh,
    const label cellOrder,
    const label faceOrder,
    const Switch allowDegenerateTriFallback
)
:
    mesh_(mesh),
    cellOrder_(cellOrder),
    faceOrder_(faceOrder),
    allowDegenerateTriFallback_(allowDegenerateTriFallback),
    faceQuadPointsPtr_(nullptr),
    faceQuadWeightsPtr_(nullptr),
    cellQuadPointsPtr_(nullptr),
    cellQuadWeightsPtr_(nullptr)
{
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

fvMeshQuadrature::~fvMeshQuadrature()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const CompactListList<point>& fvMeshQuadrature::faceQuadPoints() const
{
    if (!faceQuadPointsPtr_.valid())
    {
        calcQuadPointsAndWeights();
    }

    return faceQuadPointsPtr_();
}


const CompactListList<scalar>& fvMeshQuadrature::faceQuadWeights() const
{
    if (!faceQuadWeightsPtr_.valid())
    {
        calcQuadPointsAndWeights();
    }

    return faceQuadWeightsPtr_();
}

const CompactListList<point>& fvMeshQuadrature::cellQuadPoints() const
{
    if (!cellQuadPointsPtr_.valid())
    {
        calcQuadPointsAndWeights();
    }

    return cellQuadPointsPtr_();
}


const CompactListList<scalar>& fvMeshQuadrature::cellQuadWeights() const
{
    if (!cellQuadWeightsPtr_.valid())
    {
        calcQuadPointsAndWeights();
    }

    return cellQuadWeightsPtr_();
}

void fvMeshQuadrature::clear()
{
    faceQuadPointsPtr_.clear();
    faceQuadWeightsPtr_.clear();
    cellQuadPointsPtr_.clear();
    cellQuadWeightsPtr_.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
