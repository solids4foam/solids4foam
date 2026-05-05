/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twoDPointCorrectorS4f.H"
#include "polyMesh.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "SubField.H"
#include "meshTools.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoDPointCorrectorS4f, 0);
}

const Foam::scalar Foam::twoDPointCorrectorS4f::edgeOrthogonalityTol = 1.0 - 1e-4;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::twoDPointCorrectorS4f::calcAddressing() const
{
    // Find geometry normal
    planeNormalPtr_ = std::make_unique<vector>(0, 0, 0);
    auto& pn = *planeNormalPtr_;

    // Algorithm:
    // Attempt to find wedge patch and work out the normal from it.
    // If not found, find an empty patch with faces in it and use the
    // normal of the first face.  If neither is found, declare an
    // error.

    // Try and find a wedge patch
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const polyPatch& p : patches)
    {
        if (isA<wedgePolyPatch>(p))
        {
            isWedge_ = true;

            const wedgePolyPatch& wp = refCast<const wedgePolyPatch>(p);

            pn = wp.centreNormal();

            wedgeAxis_ = wp.axis();
            wedgeAngle_ = mag(acos(wp.cosAngle()));

            if (polyMesh::debug)
            {
                Pout<< "Found normal from wedge patch " << p.index() << nl;
            }

            break;
        }
    }

    // Try to find an empty patch with faces
    if (!isWedge_)
    {
        for (const polyPatch& p : patches)
        {
            if (isA<emptyPolyPatch>(p) && p.size())
            {
                pn = p.faceAreas()[0];

                if (polyMesh::debug)
                {
                    Pout<< "Found normal from empty patch " << p.index() << nl;
                }

                break;
            }
        }
    }


    const bool missingLocalNormal = (mag(pn) < VSMALL);

    if (Pstream::parRun() && returnReduceOr(missingLocalNormal, mesh_.comm()))
    {
        if (!missingLocalNormal)
        {
            pn /= mag(pn);

            direction maxCmpt = vector::X;
            if (mag(pn.y()) > mag(pn[maxCmpt]))
            {
                maxCmpt = vector::Y;
            }
            if (mag(pn.z()) > mag(pn[maxCmpt]))
            {
                maxCmpt = vector::Z;
            }

            if (pn[maxCmpt] < 0)
            {
                pn = -pn;
            }
        }

        reduce(pn, sumOp<vector>(), UPstream::msgType(), mesh_.comm());

        if (mag(pn) < VSMALL)
        {
            FatalErrorInFunction
                << "Cannot determine normal vector from patches."
                << abort(FatalError);
        }

        pn /= mag(pn);
    }
    else if (missingLocalNormal)
    {
        FatalErrorInFunction
            << "Cannot determine normal vector from patches."
            << abort(FatalError);
    }
    else
    {
        pn /= mag(pn);
    }

    if (polyMesh::debug)
    {
        Pout<< " twoDPointCorrectorS4f normal: " << pn << nl;
    }

    // Select edges to be included in check.
    normalEdgeIndicesPtr_ = std::make_unique<labelList>(mesh_.nEdges());
    auto& neIndices = *normalEdgeIndicesPtr_;

    const edgeList& meshEdges = mesh_.edges();

    const pointField& meshPoints = mesh_.points();

    label nNormalEdges = 0;

    forAll(meshEdges, edgeI)
    {
        const edge& e = meshEdges[edgeI];

        vector edgeVector = e.unitVec(meshPoints);

        if (mag(edgeVector & pn) > edgeOrthogonalityTol)
        {
            // this edge is normal to the plane. Add it to the list
            neIndices[nNormalEdges++] = edgeI;
        }
    }

    neIndices.setSize(nNormalEdges);

    // Construction check: number of points in a read 2-D or wedge geometry
    // should be odd and the number of edges normal to the plane should be
    // exactly half the number of points
    if (!isWedge_)
    {
        if (meshPoints.size() % 2 != 0)
        {
            WarningInFunction
                << "the number of vertices in the geometry "
                << "is odd - this should not be the case for a 2-D case. "
                << "Please check the geometry."
                << endl;
        }

        if (2*nNormalEdges != meshPoints.size())
        {
            WarningInFunction
                << "The number of points in the mesh is "
                << "not equal to twice the number of edges normal to the plane "
                << "- this may be OK only for wedge geometries.\n"
                << "    Please check the geometry or adjust "
                << "the orthogonality tolerance.\n" << endl
                << "Number of normal edges: " << nNormalEdges
                << " number of points: " << meshPoints.size()
                << endl;
        }
    }
}


void Foam::twoDPointCorrectorS4f::clearAddressing() const
{
    planeNormalPtr_.reset(nullptr);
    normalEdgeIndicesPtr_.reset(nullptr);
}


void Foam::twoDPointCorrectorS4f::snapToWedge
(
    const vector& n,
    const point& A,
    point& p
) const
{
    scalar ADash = mag(A - wedgeAxis_*(wedgeAxis_ & A));
    vector pDash = ADash*tan(wedgeAngle_)*planeNormal();

    p = A + sign(n & p)*pDash;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoDPointCorrectorS4f::twoDPointCorrectorS4f(const polyMesh& mesh)
:
    MeshObject_type(mesh),
    required_(mesh_.nGeometricD() == 2),
    isWedge_(false),
    wedgeAxis_(Zero),
    wedgeAngle_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoDPointCorrectorS4f::~twoDPointCorrectorS4f()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::direction Foam::twoDPointCorrectorS4f::normalDir() const
{
    const vector& pn = planeNormal();

    if (mag(pn.x()) >= edgeOrthogonalityTol)
    {
        return vector::X;
    }
    else if (mag(pn.y()) >= edgeOrthogonalityTol)
    {
        return vector::Y;
    }
    else if (mag(pn.z()) >= edgeOrthogonalityTol)
    {
        return vector::Z;
    }

    FatalErrorInFunction
        << "Plane normal not aligned with the coordinate system" << nl
        << "    pn = " << pn
        << abort(FatalError);

    return vector::Z;
}


const Foam::vector& Foam::twoDPointCorrectorS4f::planeNormal() const
{
    if (!planeNormalPtr_)
    {
        calcAddressing();
    }

    return *planeNormalPtr_;
}


const Foam::labelList& Foam::twoDPointCorrectorS4f::normalEdgeIndices() const
{
    if (!normalEdgeIndicesPtr_)
    {
        calcAddressing();
    }

    return *normalEdgeIndicesPtr_;
}


void Foam::twoDPointCorrectorS4f::correctPoints(pointField& p) const
{
    if (!required_) return;

    // Algorithm:
    // Loop through all edges. Calculate the average point position A for
    // the front and the back. Correct the position of point P (in two planes)
    // such that vectors AP and planeNormal are parallel

    // Get reference to edges
    const edgeList&  meshEdges = mesh_.edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();

    for (const label edgei : neIndices)
    {
        point& pStart = p[meshEdges[edgei].start()];

        point& pEnd = p[meshEdges[edgei].end()];

        // calculate average point position
        point A = 0.5*(pStart + pEnd);
        meshTools::constrainToMeshCentre(mesh_, A);

        if (isWedge_)
        {
            snapToWedge(pn, A, pStart);
            snapToWedge(pn, A, pEnd);
        }
        else
        {
            // correct point locations
            pStart = A + pn*(pn & (pStart - A));
            pEnd = A + pn*(pn & (pEnd - A));
        }
    }
}


void Foam::twoDPointCorrectorS4f::correctDisplacement
(
    const pointField& p,
    vectorField& disp
) const
{
    if (!required_) return;

    // Algorithm:
    // Loop through all edges. Calculate the average point position A for
    // the front and the back. Correct the position of point P (in two planes)
    // such that vectors AP and planeNormal are parallel

    // Get reference to edges
    const edgeList&  meshEdges = mesh_.edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();

    for (const label edgei : neIndices)
    {
        const edge& e = meshEdges[edgei];

        label startPointi = e.start();
        point pStart = p[startPointi] + disp[startPointi];

        label endPointi = e.end();
        point pEnd = p[endPointi] + disp[endPointi];

        // calculate average point position
        point A = 0.5*(pStart + pEnd);
        meshTools::constrainToMeshCentre(mesh_, A);

        if (isWedge_)
        {
            snapToWedge(pn, A, pStart);
            snapToWedge(pn, A, pEnd);
            disp[startPointi] = pStart - p[startPointi];
            disp[endPointi] = pEnd - p[endPointi];
        }
        else
        {
            // correct point locations
            disp[startPointi] = (A + pn*(pn & (pStart - A))) - p[startPointi];
            disp[endPointi] = (A + pn*(pn & (pEnd - A))) - p[endPointi];
        }
    }
}


void Foam::twoDPointCorrectorS4f::updateMesh(const mapPolyMesh&)
{
    clearAddressing();
}


bool Foam::twoDPointCorrectorS4f::movePoints()
{
    return true;
}


// ************************************************************************* //
