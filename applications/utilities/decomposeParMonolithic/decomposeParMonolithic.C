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
    decomposeParMonolithic

Description
    Decomposes multiple volume regions as one coupled graph for monolithic
    FSI solvers.  By default it also runs the standard decomposePar utility
    internally (via method manual) so that a single command produces the
    full processor directory tree for every region.

    This avoids the inefficiency of independent per-region decomposition
    when one region is much smaller than the other. With a coupled
    decomposition some ranks may have zero cells in one region, which the
    monolithic solver already supports.

Usage
    \b decomposeParMonolithic [OPTIONS]

    Options:
      - \par -decomposeParDict \<file\>
        Alternative decomposeParDict file with monolithicCoeffs section.

      - \par -cellDist
        Write cell distribution as a volScalarField for visualisation.

      - \par -decompose-only
        Only compute and write cellDecomposition files; do not run
        decomposePar to create processor directories.

      - \par -force
        Remove existing processor subdirectories before decomposing.

      - \par -copy-zero
        Copy 0 directory to processor directories rather than decompose
        the fields.

      - \par -no-fields
        Suppress conversion of fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "labelIOList.H"
#include "decompositionMethod.H"
#include "CompactListList.H"
#include "volFields.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Build face-centre matching between two conformal patches.
// Returns list of (ownerCellA, ownerCellB) pairs.
List<labelPair> matchConformalInterface
(
    const fvMesh& meshA,
    const word& patchNameA,
    const fvMesh& meshB,
    const word& patchNameB
)
{
    const label patchIA = meshA.boundaryMesh().findPatchID(patchNameA);
    if (patchIA < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchNameA << " not found in region "
            << meshA.name()
            << exit(FatalError);
    }

    const label patchIB = meshB.boundaryMesh().findPatchID(patchNameB);
    if (patchIB < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchNameB << " not found in region "
            << meshB.name()
            << exit(FatalError);
    }

    const polyPatch& ppA = meshA.boundaryMesh()[patchIA];
    const polyPatch& ppB = meshB.boundaryMesh()[patchIB];

    const vectorField cfA(ppA.faceCentres());
    const vectorField cfB(ppB.faceCentres());

    if (cfA.size() != cfB.size())
    {
        FatalErrorInFunction
            << "Interface patches have different sizes: "
            << patchNameA << " (" << cfA.size() << ") vs "
            << patchNameB << " (" << cfB.size() << ")" << nl
            << "Non-conformal interfaces are not supported in v1."
            << exit(FatalError);
    }

    // Compute a length scale for tolerance
    scalar tolSq = 0;
    if (cfA.size() > 0)
    {
        // Use a fraction of the mesh bounding box diagonal
        const boundBox bbA(cfA);
        const boundBox bbB(cfB);
        const scalar diagA = bbA.mag();
        const scalar diagB = bbB.mag();
        const scalar charLen = max(diagA, diagB);
        // Tolerance: 1e-6 of characteristic length
        tolSq = sqr(1e-6*charLen);

        if (tolSq < SMALL)
        {
            tolSq = SMALL;
        }
    }

    List<labelPair> pairs(cfA.size());

    // For each face in patchA, find the matching face in patchB
    // by nearest face centre. For conformal meshes this should be exact.
    labelList usedB(cfB.size(), -1);

    forAll(cfA, faceI)
    {
        const point& ptA = cfA[faceI];
        label bestJ = -1;
        scalar bestDistSq = GREAT;

        forAll(cfB, faceJ)
        {
            if (usedB[faceJ] >= 0) continue;

            const scalar dSq = magSqr(ptA - cfB[faceJ]);
            if (dSq < bestDistSq)
            {
                bestDistSq = dSq;
                bestJ = faceJ;
            }
        }

        if (bestJ < 0 || bestDistSq > tolSq)
        {
            FatalErrorInFunction
                << "Could not find matching face for face " << faceI
                << " of patch " << patchNameA
                << " at " << ptA << nl
                << "Best distance: " << Foam::sqrt(bestDistSq)
                << ", tolerance: " << Foam::sqrt(tolSq)
                << exit(FatalError);
        }

        usedB[bestJ] = faceI;

        // Get owner cells
        const label cellA = ppA.faceCells()[faceI];
        const label cellB = ppB.faceCells()[bestJ];

        pairs[faceI] = labelPair(cellA, cellB);
    }

    return pairs;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Decompose multiple regions as one coupled graph for monolithic FSI"
    );

    argList::noParallel();

    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Alternative decomposeParDict file with monolithicCoeffs"
    );

    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as volScalarField for visualisation"
    );

    argList::addBoolOption
    (
        "decompose-only",
        "Only write cellDecomposition files; skip decomposePar"
    );

    argList::addBoolOption
    (
        "force",
        "Remove existing processor*/ subdirs before decomposing"
    );

    argList::addBoolOption
    (
        "copy-zero",
        "Copy 0/ directory to processor*/ rather than decompose fields"
    );

    argList::addBoolOption
    (
        "no-fields",
        "Suppress conversion of fields"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool writeCellDist = args.found("cellDist");
    const bool decomposeOnly = args.found("decompose-only");
    const bool forceOverwrite = args.found("force");
    const bool copyZero = args.found("copy-zero");
    const bool noFields = args.found("no-fields");

    // -----------------------------------------------------------------------
    // Read dictionary
    // -----------------------------------------------------------------------

    fileName decompDictFile;
    if
    (
        args.readIfPresent("decomposeParDict", decompDictFile)
     && !decompDictFile.empty() && !decompDictFile.isAbsolute()
    )
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    // Find the dictionary
    IOdictionary decompDict
    (
        IOobject::selectIO
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            decompDictFile
        )
    );

    const label nDomains = decompDict.get<label>("numberOfSubdomains");

    Info<< "Number of subdomains: " << nDomains << nl << endl;

    // Read monolithicCoeffs
    const dictionary& monoDict = decompDict.subDict("monolithicCoeffs");
    const wordList regionNames(monoDict.get<wordList>("regions"));

    if (regionNames.size() < 2)
    {
        FatalErrorInFunction
            << "monolithicCoeffs::regions must list at least 2 regions, got "
            << regionNames.size()
            << exit(FatalError);
    }

    Info<< "Regions: " << regionNames << nl << endl;

    // Optional region weights
    dictionary regionWeightsDict;
    if (monoDict.found("regionWeights"))
    {
        regionWeightsDict = monoDict.subDict("regionWeights");
    }

    // Read interface definitions
    const List<dictionary> interfaceDicts
    (
        monoDict.lookup("interfaces")
    );

    // -----------------------------------------------------------------------
    // Load region meshes
    // -----------------------------------------------------------------------

    PtrList<fvMesh> regionMeshes(regionNames.size());

    labelList regionOffset(regionNames.size(), 0);
    label nTotalCells = 0;

    forAll(regionNames, ri)
    {
        Info<< "Loading mesh for region: " << regionNames[ri] << endl;

        regionMeshes.set
        (
            ri,
            new fvMesh
            (
                IOobject
                (
                    regionNames[ri],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        regionOffset[ri] = nTotalCells;
        nTotalCells += regionMeshes[ri].nCells();

        Info<< "  Cells: " << regionMeshes[ri].nCells()
            << "  (global offset: " << regionOffset[ri] << ")" << endl;
    }

    Info<< nl << "Total cells across all regions: " << nTotalCells << nl
        << endl;

    // -----------------------------------------------------------------------
    // Build combined adjacency graph
    // -----------------------------------------------------------------------

    Info<< "Building combined adjacency graph" << endl;

    // Start with internal adjacency per region
    labelListList globalCellCells(nTotalCells);

    forAll(regionMeshes, ri)
    {
        const fvMesh& mesh = regionMeshes[ri];
        const label offset = regionOffset[ri];
        const labelListList& cc = mesh.cellCells();

        forAll(cc, cellI)
        {
            const labelList& nbrs = cc[cellI];
            labelList& globalNbrs = globalCellCells[offset + cellI];
            globalNbrs.setSize(nbrs.size());

            forAll(nbrs, ni)
            {
                globalNbrs[ni] = nbrs[ni] + offset;
            }
        }
    }

    // Add cross-region interface edges
    forAll(interfaceDicts, ii)
    {
        const dictionary& ifDict = interfaceDicts[ii];

        const word regionAName = ifDict.get<word>("regionA");
        const word patchAName = ifDict.get<word>("patchA");
        const word regionBName = ifDict.get<word>("regionB");
        const word patchBName = ifDict.get<word>("patchB");

        // Find region indices
        const label riA = regionNames.find(regionAName);
        const label riB = regionNames.find(regionBName);

        if (riA < 0)
        {
            FatalErrorInFunction
                << "Interface regionA '" << regionAName
                << "' not found in regions list"
                << exit(FatalError);
        }
        if (riB < 0)
        {
            FatalErrorInFunction
                << "Interface regionB '" << regionBName
                << "' not found in regions list"
                << exit(FatalError);
        }

        Info<< "Matching interface: " << regionAName << "/" << patchAName
            << " <-> " << regionBName << "/" << patchBName << endl;

        const List<labelPair> pairs = matchConformalInterface
        (
            regionMeshes[riA], patchAName,
            regionMeshes[riB], patchBName
        );

        Info<< "  Matched " << pairs.size() << " face pairs" << endl;

        const label offsetA = regionOffset[riA];
        const label offsetB = regionOffset[riB];

        // Add bidirectional edges
        forAll(pairs, pi)
        {
            const label globalA = pairs[pi].first() + offsetA;
            const label globalB = pairs[pi].second() + offsetB;

            // Add B to A's neighbours
            labelList& nbrsA = globalCellCells[globalA];
            nbrsA.append(globalB);

            // Add A to B's neighbours
            labelList& nbrsB = globalCellCells[globalB];
            nbrsB.append(globalA);
        }
    }

    // -----------------------------------------------------------------------
    // Build combined cell centres and weights
    // -----------------------------------------------------------------------

    pointField combinedCC(nTotalCells);
    scalarField combinedWeights(nTotalCells, 1.0);

    forAll(regionMeshes, ri)
    {
        const fvMesh& mesh = regionMeshes[ri];
        const label offset = regionOffset[ri];
        const vectorField& cc = mesh.cellCentres();

        // Optionally apply region weight
        scalar regionWeight = 1.0;
        if (regionWeightsDict.size() > 0)
        {
            regionWeightsDict.readIfPresent(regionNames[ri], regionWeight);
        }

        forAll(cc, cellI)
        {
            combinedCC[offset + cellI] = cc[cellI];
            combinedWeights[offset + cellI] = regionWeight;
        }
    }

    // -----------------------------------------------------------------------
    // Call decomposition method
    // -----------------------------------------------------------------------

    Info<< nl << "Performing coupled decomposition" << endl;

    autoPtr<decompositionMethod> methodPtr =
        decompositionMethod::New(decompDict);

    const labelList cellToProc = methodPtr->decompose
    (
        globalCellCells,
        combinedCC,
        combinedWeights
    );

    // -----------------------------------------------------------------------
    // Split result per region and write cellDecomposition
    // -----------------------------------------------------------------------

    Info<< nl << "Writing per-region cellDecomposition files" << nl << endl;

    forAll(regionMeshes, ri)
    {
        const fvMesh& mesh = regionMeshes[ri];
        const label offset = regionOffset[ri];
        const label nCells = mesh.nCells();

        labelList regionCellToProc
        (
            SubList<label>(cellToProc, nCells, offset)
        );

        labelIOList cellDecomposition
        (
            IOobject
            (
                "cellDecomposition",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            regionCellToProc
        );

        cellDecomposition.write();

        Info<< "  Region " << regionNames[ri] << ": wrote "
            << cellDecomposition.objectRelPath() << endl;

        // Optional: write as volScalarField
        if (writeCellDist)
        {
            volScalarField cellDist
            (
                IOobject
                (
                    "cellDist",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar("cellDist", dimless, 0)
            );

            forAll(regionCellToProc, cellI)
            {
                cellDist[cellI] = regionCellToProc[cellI];
            }

            cellDist.write();

            Info<< "  Region " << regionNames[ri] << ": wrote "
                << cellDist.objectRelPath() << endl;
        }
    }

    // -----------------------------------------------------------------------
    // Print statistics
    // -----------------------------------------------------------------------

    Info<< nl << "Decomposition statistics:" << nl << endl;

    // Per-processor per-region cell counts
    labelListList procRegionCells(nDomains);
    forAll(procRegionCells, pi)
    {
        procRegionCells[pi].setSize(regionNames.size(), 0);
    }

    forAll(regionMeshes, ri)
    {
        const label offset = regionOffset[ri];
        const label nCells = regionMeshes[ri].nCells();

        for (label ci = 0; ci < nCells; ++ci)
        {
            const label proc = cellToProc[offset + ci];
            procRegionCells[proc][ri]++;
        }
    }

    // Print table header
    Info<< "    Proc";
    forAll(regionNames, ri)
    {
        Info<< "\t" << regionNames[ri];
    }
    Info<< "\tTotal" << nl;

    Info<< "    ----";
    forAll(regionNames, ri)
    {
        Info<< "\t----";
    }
    Info<< "\t-----" << nl;

    labelList maxCellsPerRegion(regionNames.size(), 0);
    labelList totalCellsPerRegion(regionNames.size(), 0);

    forAll(procRegionCells, pi)
    {
        Info<< "    " << pi;

        label procTotal = 0;

        forAll(regionNames, ri)
        {
            const label n = procRegionCells[pi][ri];
            Info<< "\t" << n;
            procTotal += n;
            maxCellsPerRegion[ri] = max(maxCellsPerRegion[ri], n);
            totalCellsPerRegion[ri] += n;
        }

        Info<< "\t" << procTotal << nl;
    }

    // Zero-cell rank counts
    Info<< nl;
    forAll(regionNames, ri)
    {
        label nZeroRanks = 0;
        forAll(procRegionCells, pi)
        {
            if (procRegionCells[pi][ri] == 0)
            {
                nZeroRanks++;
            }
        }

        const scalar avgCells =
            scalar(totalCellsPerRegion[ri]) / nDomains;

        Info<< "  Region " << regionNames[ri] << ":"
            << "  total=" << totalCellsPerRegion[ri]
            << "  avg=" << avgCells
            << "  max=" << maxCellsPerRegion[ri];

        if (avgCells > SMALL)
        {
            Info<< " (" << 100.0*(maxCellsPerRegion[ri]-avgCells)/avgCells
                << "% above avg)";
        }

        Info<< "  zero-cell ranks=" << nZeroRanks << nl;
    }

    // -----------------------------------------------------------------------
    // Run decomposePar for each region (unless -decompose-only)
    // -----------------------------------------------------------------------

    if (decomposeOnly)
    {
        Info<< nl
            << "Decompose-only mode: skipping processor mesh generation."
            << nl
            << "To complete the decomposition, run decomposePar for each"
            << " region with method=manual, e.g.:" << nl << nl;

        forAll(regionNames, ri)
        {
            Info<< "  decomposePar -region " << regionNames[ri]
                << " -decomposeParDict system/decomposeParDict.manual"
                << nl;
        }

        Info<< nl
            << "where system/decomposeParDict.manual contains:" << nl
            << "  numberOfSubdomains " << nDomains << ";" << nl
            << "  method manual;" << nl
            << "  manualCoeffs { dataFile \"cellDecomposition\"; }" << nl
            << endl;
    }
    else
    {
        // Write a temporary manual decomposeParDict
        const fileName tmpDictPath =
            runTime.path()/"system"/"decomposeParDict.monolithic_tmp_";

        {
            OFstream os(tmpDictPath);

            os  << "FoamFile" << nl
                << "{" << nl
                << "    version     2.0;" << nl
                << "    format      ascii;" << nl
                << "    class       dictionary;" << nl
                << "    object      decomposeParDict;" << nl
                << "}" << nl << nl
                << "numberOfSubdomains " << nDomains << ";" << nl << nl
                << "method          manual;" << nl << nl
                << "manualCoeffs" << nl
                << "{" << nl
                << "    dataFile    \"cellDecomposition\";" << nl
                << "}" << nl;
        }

        Info<< nl
            << "Running decomposePar for each region" << nl << endl;

        forAll(regionNames, ri)
        {
            Info<< "--- Decomposing region: " << regionNames[ri]
                << " ---" << nl << endl;

            std::string cmd = "decomposePar"
                " -region " + regionNames[ri]
                + " -decomposeParDict " + tmpDictPath;

            if (forceOverwrite)
            {
                cmd += " -force";
            }
            if (copyZero)
            {
                cmd += " -copyZero";
            }
            if (noFields)
            {
                cmd += " -no-fields";
            }

            const int ret = Foam::system(cmd);

            if (ret != 0)
            {
                FatalErrorInFunction
                    << "decomposePar failed for region "
                    << regionNames[ri]
                    << " (exit code " << ret << ")"
                    << exit(FatalError);
            }
        }

        // Clean up temporary dict
        Foam::rm(tmpDictPath);

        // Clean up cellDecomposition files
        forAll(regionMeshes, ri)
        {
            const fvMesh& mesh = regionMeshes[ri];
            const fileName cellDecompPath =
                mesh.time().path()
              / mesh.facesInstance()
              / mesh.dbDir()
              / "cellDecomposition";

            Foam::rm(cellDecompPath);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
