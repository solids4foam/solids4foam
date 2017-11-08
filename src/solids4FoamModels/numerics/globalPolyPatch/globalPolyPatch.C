/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "globalPolyPatch.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPolyPatch, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalPolyPatch::calcGlobalPatch() const
{
    if (debug)
    {
        InfoIn("void globalPolyPatch::calcGlobalPatch() const")
            << "Calculating primitive patch"
            << endl;
    }

    if (globalPatchPtr_ || pointToGlobalAddrPtr_ || faceToGlobalAddrPtr_)
    {
        FatalErrorIn
        (
            "void globalPolyPatch::calcGlobalPatch() const"
        )   << "primitive face zone patch and addressing already calculated"
            << abort(FatalError);
    }

    // Get patch
    polyPatchID patchID
    (
        patchName_,
        mesh_.boundaryMesh()
    );

    if (!patchID.active())
    {
        FatalErrorIn("void globalPolyPatch::calcGlobalPatch() const")
            << "Cannot find patch " << patchName_
            << abort(FatalError);
    }

    // Collect points and faces from all processors
    typedef List<point> pointList;
    typedef List<pointList> pointListList;

    pointListList procPatchPoints(Pstream::nProcs());
    faceListList procPatchFaces(Pstream::nProcs());

    // Add points and faces if the patch is not empty
    if (!mesh_.boundaryMesh()[patchID.index()].empty())
    {
        // Insert my points
        procPatchPoints[Pstream::myProcNo()] =
            mesh_.boundaryMesh()[patchID.index()].localPoints();

        // Insert my faces
        procPatchFaces[Pstream::myProcNo()] =
            mesh_.boundaryMesh()[patchID.index()].localFaces();
    }

    // Communicate points
    Pstream::gatherList(procPatchPoints);
    Pstream::scatterList(procPatchPoints);

    // Communicate faces
    Pstream::gatherList(procPatchFaces);
    Pstream::scatterList(procPatchFaces);

    // At this point, all points and faces for the current patch
    // are available.

    // Count the number of faces in the global face zone
    label nZoneFaces = 0;
    label nZonePoints = 0;

    // Sum up points and faces to add
    forAll (procPatchFaces, procI)
    {
        nZonePoints += procPatchPoints[procI].size();
        nZoneFaces += procPatchFaces[procI].size();
    }

    Info<< "Global zone size for patch " << patchID.name()
        << ": " << nZoneFaces << endl;

    if (nZoneFaces == 0)
    {
        FatalErrorIn("void globalPolyPatch::calcGlobalPatch() const")
            << "Patch " << patchID.name()
            << " appears to be globally empty.  "
            << "Please check definition."
            << abort(FatalError);
    }

    // Record current points to add
    pointField zonePoints(nZonePoints);
    faceList zoneFaces(nZoneFaces);
    
    label nCurPoints = 0;
    label nCurFaces = 0;

    // Collect all points and faces
    forAll (procPatchPoints, procI)
    {
        // Add points from all processors except self
        const pointList& curProcPoints = procPatchPoints[procI];

        // Label point map for the current processor
        labelList pointMap(curProcPoints.size());

        // Add points from all processors
        forAll (curProcPoints, pointI)
        {
            // Note: possible removal of duplicate points here
            // HJ, 28/Dec/2016

            // Add the point
            zonePoints[nCurPoints] = curProcPoints[pointI];

            // Record point mapping
            pointMap[pointI] = nCurPoints;
            nCurPoints++;
        }

        // Add faces from all processors
        const faceList& curProcFaces = procPatchFaces[procI];

        // Label face map for the current processor
        labelList faceMap(curProcFaces.size());

        forAll (curProcFaces, faceI)
        {
            // Renumber face into new points
            face curFace = curProcFaces[faceI];

            forAll (curFace, fI)
            {
                curFace[fI] = pointMap[curFace[fI]];
            }

            // Record the face into zone
            zoneFaces[nCurFaces] = curFace;
            faceMap[faceI] = nCurFaces;
            nCurFaces++;
        }

        if (procI == Pstream::myProcNo())
        {
            // Store point addressing
            pointToGlobalAddrPtr_ = new labelList(pointMap);

            // Store face addressing
            faceToGlobalAddrPtr_ = new labelList(faceMap);
        }
    }

    // All points and faces are collected.  Make a patch
    globalPatchPtr_ = new standAlonePatch(zoneFaces, zonePoints);

    if (debug)
    {
        Info<< "void globalPolyPatch::calcGlobalPatch() const : "
            << "Finished calculating primitive patch"
            << endl;
    }
}


void Foam::globalPolyPatch::check() const
{
    label patchIndex = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchIndex < 0)
    {
        FatalErrorIn("void globalPolyPatch::check() const")
            << "Patch " << patchName_ << " not found."
            << abort(FatalError);
    }
}


void Foam::globalPolyPatch::clearOut() const
{
    deleteDemandDrivenData(globalPatchPtr_);

    deleteDemandDrivenData(pointToGlobalAddrPtr_);
    deleteDemandDrivenData(faceToGlobalAddrPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::globalPolyPatch::globalPolyPatch
(
    const word& patchName,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    patchName_(patchName),
    patch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    globalPatchPtr_(NULL),
    pointToGlobalAddrPtr_(NULL),
    faceToGlobalAddrPtr_(NULL)
{
    check();
}


// Construct from dictionary
Foam::globalPolyPatch::globalPolyPatch
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    patchName_(dict.lookup("patch")),
    patch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    globalPatchPtr_(NULL),
    pointToGlobalAddrPtr_(NULL),
    faceToGlobalAddrPtr_(NULL)
{
    check();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalPolyPatch::~globalPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::globalPolyPatch::mesh() const
{
    return mesh_;
}


const Foam::standAlonePatch& Foam::globalPolyPatch::globalPatch() const
{
    if (!globalPatchPtr_)
    {
        calcGlobalPatch();
    }

    return *globalPatchPtr_;
}


const Foam::labelList& Foam::globalPolyPatch::pointToGlobalAddr() const
{
    if (!pointToGlobalAddrPtr_)
    {
        calcGlobalPatch();
    }

    return *pointToGlobalAddrPtr_;
}


const Foam::labelList& Foam::globalPolyPatch::faceToGlobalAddr() const
{
    if (!faceToGlobalAddrPtr_)
    {
        calcGlobalPatch();
    }

    return *faceToGlobalAddrPtr_;
}


void Foam::globalPolyPatch::updateMesh()
{
    clearOut();
}


void Foam::globalPolyPatch::movePoints(const pointField& p)
{
    if (globalPatchPtr_)
    {
        globalPatchPtr_->movePoints(p);
    }
}


// ************************************************************************* //
