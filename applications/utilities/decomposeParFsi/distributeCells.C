/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "decompositionMethod.H"
#include "cpuTime.H"
#include "cyclicPolyPatch.H"
#include "cellSet.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void domainDecomposition::distributeCells()
{
    Info<< "\nCalculating distribution of cells" << endl;

    cpuTime decompositionTime;


    // See if any faces need to have owner and neighbour on same processor
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet sameProcFaces;

    if (decompositionDict_.found("preservePatches"))
    {
        wordList pNames(decompositionDict_.lookup("preservePatches"));

        Info<< "Keeping owner and neighbour of faces in patches " << pNames
            << " on same processor" << endl;

        const polyBoundaryMesh& patches = boundaryMesh();

        forAll(pNames, i)
        {
            label patchI = patches.findPatchID(pNames[i]);

            if (patchI == -1)
            {
                FatalErrorIn("domainDecomposition::distributeCells()")
                    << "Unknown preservePatch " << pNames[i]
                    << nl << "Valid patches are " << patches.names()
                    << exit(FatalError);
            }

            const polyPatch& pp = patches[patchI];

            forAll(pp, i)
            {
                sameProcFaces.insert(pp.start() + i);
            }
        }
    }
    if (decompositionDict_.found("preserveFaceZones"))
    {
        wordList zNames(decompositionDict_.lookup("preserveFaceZones"));

        Info<< "Keeping owner and neighbour of faces in zones " << zNames
            << " on same processor" << endl;

        const faceZoneMesh& fZones = faceZones();

        forAll(zNames, i)
        {
            label zoneI = fZones.findZoneID(zNames[i]);

            if (zoneI == -1)
            {
                FatalErrorIn("domainDecomposition::distributeCells()")
                    << "Unknown preserveFaceZone " << zNames[i]
                    << endl << "Valid faceZones are " << fZones.names()
                    << exit(FatalError);
            }

            const faceZone& fz = fZones[zoneI];

            forAll(fz, i)
            {
                sameProcFaces.insert(fz[i]);
            }
        }
    }


    // Construct decomposition method and either do decomposition on
    // cell centres or on agglomeration


    autoPtr<decompositionMethod> decomposePtr = decompositionMethod::New
    (
        decompositionDict_,
        *this
    );

    if (sameProcFaces.empty())
    {
        cellToProc_ = decomposePtr().decompose(cellCentres());
    }
    else
    {
        Info<< "Selected " << sameProcFaces.size()
            << " faces whose owner and neighbour cell should be kept on the"
            << " same processor" << endl;

        // Faces where owner and neighbour are not 'connected' (= all except
        // sameProcFaces)
        boolList blockedFace(nFaces(), true);

        forAllConstIter(labelHashSet, sameProcFaces, iter)
        {
            blockedFace[iter.key()] = false;
        }

        // Connect coupled boundary faces
        const polyBoundaryMesh& patches =  boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    blockedFace[pp.start() + i] = false;
                }
            }
        }

        // Determine global regions, separated by blockedFaces
        regionSplit globalRegion(*this, blockedFace);


        // Determine region cell centres
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // This just takes the first cell in the region. Otherwise the problem
        // is with cyclics - if we'd average the region centre might be
        // somewhere in the middle of the domain which might not be anywhere
        // near any of the cells.

        const point greatPoint(GREAT, GREAT, GREAT);

        pointField regionCentres(globalRegion.nRegions(), greatPoint);

        forAll(globalRegion, cellI)
        {
            label regionI = globalRegion[cellI];

            if (regionCentres[regionI] == greatPoint)
            {
                regionCentres[regionI] = cellCentres()[cellI];
            }
        }

        // Do decomposition on agglomeration
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cellToProc_ = decomposePtr().decompose(globalRegion, regionCentres);
    }

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}


// ************************************************************************* //
