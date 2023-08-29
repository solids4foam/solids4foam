/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "solidPolyMesh.H"
#include "demandDrivenData.H"
#include "SortableList.H"
#include "pointFields.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "processorFvPatch.H"
#include "BlockLduMatrix.H"
#include "CompactListList.H"
#include "PackedList.H"
#include "blockFixedDisplacementFvPatchVectorField.H"
#include "blockFixedDisplacementZeroShearFvPatchVectorField.H"
//#include "solidSymmetryFvPatchVectorField.H"
#include "emptyFvPatchFields.H"
//#include "blockGlobalPolyPatch.H"
//#include "blockGlobalFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::solidPolyMesh, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void solidPolyMesh::addGlobalPatch()
// {
//     polyMesh& mesh = fMesh_;

//     // Create new list of patches
//     const polyBoundaryMesh& patches = mesh.boundaryMesh();
//     DynamicList<polyPatch*> newPatches(patches.size() + 1);

//     // Copy old patches
//     label startFaceI = mesh.nInternalFaces();
//     forAll(patches, patchI)
//     {
//         const polyPatch& pp = patches[patchI];

//         newPatches.append
//         (
//             pp.clone
//             (
//                 patches,
//                 patchI,
//                 pp.size(),
//                 startFaceI
//             ).ptr()
//         );

//         startFaceI += pp.size();
//     }

//     // Add new patch
//     newPatches.append
//     (
//         new blockGlobalPolyPatch
//         (
//             "blockGlobalSolid",
//             0,
//             startFaceI,
//             patches.size(),
//             patches
//         )
//         // polyPatch::New
//         // (
//         //     coupledPolyPatch::typeName, // patch type
//         //     "blockGlobalSolid",  // name
//         //     0,              // size
//         //     startFaceI,     // start
//         //     patches.size(), // patch index
//         //     patches         // polyBoundaryMesh
//         // ).ptr()
//     );

//     // Remove the old patches and add the new patches
//     newPatches.shrink();
//     mesh.removeBoundary();
//     //mesh.addPatches(newPatches);

//     fMesh_.removeFvBoundary();
//     fMesh_.addFvPatches(newPatches, true);
// }


void solidPolyMesh::clearOut() const
{
    deleteDemandDrivenData(lduPtr_);
    deleteDemandDrivenData(implicitBondsPtr_);
    deleteDemandDrivenData(faceImplicitBondsPtr_);
    deleteDemandDrivenData(fvMeshAddressingMapPtr_);
    deleteDemandDrivenData(cellPointCellsPtr_);
    deleteDemandDrivenData(cellPointBFacesPtr_);
    deleteDemandDrivenData(bFacePointBFacesPtr_);
    deleteDemandDrivenData(cellImplicitBondsPtr_);
    deleteDemandDrivenData(cellBFaceImplicitBondsPtr_);
    deleteDemandDrivenData(bFaceImplicitBondsPtr_);
    deleteDemandDrivenData(pointProcFacesCoeffsPtr_);
    deleteDemandDrivenData(pointProcBndFacesCoeffsPtr_);
    deleteDemandDrivenData(pointProcCellsCoeffsPtr_);
    deleteDemandDrivenData(gPtNgbProcBndFaceCoeffsPtr_);
    deleteDemandDrivenData(gPtNgbProcCellCoeffsPtr_);
    deleteDemandDrivenData(pointFixedComponentPtr_);
    deleteDemandDrivenData(pointFixedDirectionPtr_);
}


void solidPolyMesh::calcImplicitBonds() const
{
    if (debug)
    {
        Info<< "void solidPolyMesh::calcImplicitBonds() const "
            << "creating implicit bonds" << endl;
    }

    // It is an error to attempt to recalculate if the pointer is already set
    if
    (
        implicitBondsPtr_
     || faceImplicitBondsPtr_
     || cellPointBFacesPtr_
     || bFacePointBFacesPtr_
    )
    {
        FatalErrorIn("solidPolyMesh::calcImplicitBonds() const")
            << "implicit bonds, near implicit bonds, cell point cells, "
            << "coupled implicit bonds or coupled shadow patch size "
            << "lists already exist" << abort(FatalError);
    }

    // Define implicit bonds connecting cell centres and boundary face centres
    // to the cell centres and boundary face centres of surrounding
    // point-cell/point-boundary-faces neighbours
    // We will follow a similar procedure to primitiveMesh::calcEdges() to
    // calculate the edges.

    // Count the maximun number of implicit bonds
    label maxEdges = 0;

    const polyMesh& mesh = (*this)();

    const cellList& c = mesh.cells();
    const labelListList& pc = mesh.pointCells();
    const labelListList& cp = mesh.cellPoints();
    const labelListList& cc = mesh.cellCells();
    const labelListList& pf = mesh.pointFaces();
    const faceList& f = mesh.faces();
    const label nCells = mesh.nCells();
    const label nVar = this->nVariables();

    // We will create the list of implicit neighbours
    // Implicit neighbours are the neighbours cells/boundary-faecs connected
    // through a point

    // Cell neighbours of the cells
    labelListList implicitNeiOfCell(c.size());

    // Boundary face neighbours of the cells
    labelListList implicitBFaceOfCell(c.size());

    // Cells neighbours of the boundary faces
    labelListList implicitNeiOfBFace(nVar - nCells);

    // Create cellPointCells array
    cellPointCellsPtr_ = new labelListList(nCells);
    labelListList& cellPointCells = *cellPointCellsPtr_;
    cellPointBFacesPtr_ = new labelListList(nCells);
    labelListList& cellPointBFaces = *cellPointBFacesPtr_;

    // Count and store the neighbours

    forAll(c, cellI)
    {
        // Count cell cell-to-cell implicitBonds and store the neighbour IDs
        label ie = 0;
        labelHashSet curCellPointCells(implicitBondsPerCell_);

        // Count cell-to-boundary-face implicitBonds
        label bie = 0;
        labelHashSet curCellPointBFaces(implicitBondsPerCell_);

        forAll(cp[cellI], cellPointI)
        {
            // Count point cells for every point in the cell.
            // We don't count the current cell i.e. if there are 4
            // pointCells then we have 3 implicit bonds from our original
            // cell
            forAll(pc[cp[cellI][cellPointI]], pointCellI)
            {
                // Add the cell if it is not the original cell and it has
                // not already been added

                const label otherCellID = pc[cp[cellI][cellPointI]][pointCellI];

                if (otherCellID != cellI)
                {
                    if (!curCellPointCells.found(otherCellID))
                    {
                        curCellPointCells.insert(otherCellID);
                        ie++;
                    }
                }
            }

            // Count cell-to-boundary-face implicitBonds for every cell
            forAll(pf[cp[cellI][cellPointI]], pointFaceI)
            {
                // Add the boundary face if it has not already been added

                const label bFaceID = pf[cp[cellI][cellPointI]][pointFaceI];

                if (!mesh.isInternalFace(bFaceID))
                {
                    word pType =
                        mesh.boundaryMesh()
                        [
                            mesh.boundaryMesh().whichPatch(bFaceID)
                        ].type();

                    if
                    (
                        pType != emptyPolyPatch::typeName
                        && pType != processorPolyPatch::typeName
                    )
                    {
                        label varID = findOldVariableID(bFaceID);
                        if (!curCellPointBFaces.found(varID))
                        {
                            curCellPointBFaces.insert(varID);
                            bie++;
                        }
                    }
                }
            }
        }

        maxEdges += ie;
        maxEdges += bie;

        // Record implicit neighbours

        labelList cellImplicitNeiCells(ie);

        labelList& cpc = cellPointCells[cellI];
        cpc.setSize(curCellPointCells.size(), -1);
        cpc = curCellPointCells.toc();

        forAll(cellImplicitNeiCells, imBondI)
        {
            label otherCellID = cpc[imBondI];
            if (otherCellID != cellI)
            {
                cellImplicitNeiCells[imBondI] = otherCellID;
            }
        }

        implicitNeiOfCell[cellI] = cellImplicitNeiCells;

        labelList cellBFaceImplicitBFaces(bie);

        labelList& cpbf = cellPointBFaces[cellI];
        cpbf.setSize(curCellPointBFaces.size(), -1);
        cpbf = curCellPointBFaces.toc();

        forAll(cellBFaceImplicitBFaces, imBondI)
        {
            label faceID = cpbf[imBondI];
            cellBFaceImplicitBFaces[imBondI] = faceID;
        }

        implicitBFaceOfCell[cellI] = cellBFaceImplicitBFaces;
    }

    // Count boundary-face-to-boundary-face implicitBonds
    bFacePointBFacesPtr_ = new labelListList(nVar - nCells);
    labelListList& bFacePointBFaces = *bFacePointBFacesPtr_;

    forAll(mesh.boundaryMesh(), patchI)
    {
        const word pType = mesh.boundaryMesh()[patchI].type();

        if
        (
            pType != emptyPolyPatch::typeName
         && pType != processorPolyPatch::typeName
        )
        {
            const label start = mesh.boundaryMesh()[patchI].start();

            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
                label bbie = 0;
                labelHashSet curFaceBFaces(implicitBondsPerCell_);

                // Count boundary-face-to-boundary-face implicitBonds for every
                // boundary face
                const label faceID = start + faceI;

                forAll(f[faceID], facePointI)
                {
                    forAll(pf[f[faceID][facePointI]], pfI)
                    {
                        const label bFaceID = pf[f[faceID][facePointI]][pfI];

                        if
                        (
                            bFaceID != faceID
                         && !mesh.isInternalFace(bFaceID)
                        )
                        {
                            word otherPType =
                                mesh.boundaryMesh()
                                [
                                    mesh.boundaryMesh().whichPatch(bFaceID)
                                ].type();

                            if
                            (
                                otherPType != emptyPolyPatch::typeName
                                && otherPType != processorPolyPatch::typeName
                            )
                            {
                                label otherVarID = findOldVariableID(bFaceID);

                                // Add the boundary face if it has not already
                                // been added
                                if (!curFaceBFaces.found(otherVarID))
                                {
                                    curFaceBFaces.insert(otherVarID);
                                    bbie++;
                                }
                            }
                        }
                    }
                }

                maxEdges += bbie;

                // Record boundary face implicit neighbours for the current
                // boundary face
                labelList bFaceImplicitNeiFaces(bbie, -1);

                const label bVarID = findOldVariableID(faceID) - nCells;

                labelList& fpf = bFacePointBFaces[bVarID];
                fpf.setSize(curFaceBFaces.size(), -1);
                fpf = curFaceBFaces.toc();

                forAll(bFaceImplicitNeiFaces, bondI)
                {
                    bFaceImplicitNeiFaces[bondI] = fpf[bondI];
                }

                if (min(bFaceImplicitNeiFaces) == -1)
                {
                    FatalErrorIn("solidPolyMesh::calcImplicitBonds()")
                        << "problem calculating implicit bonds" << nl
                        << "bFaceImplicitNeiFaces = " << bFaceImplicitNeiFaces
                        << abort(FatalError);
                }

                implicitNeiOfBFace[bVarID] = bFaceImplicitNeiFaces;
            }
        }
    }


    // Create implicit bonds


    // The implicitBonds will be ordered in an analagous way to the faces in a
    // standard fvMesh i.e. we start at cell 0 and add all implicitBonds, then
    // go to cell 1 and add all implicitBonds not already added, then cell 2
    // etc.
    // Within a cell we add the implicitBonds in order of increasing neighbour
    // cell ID
    // One key difference with fvMesh is that boundary faces are included
    // implicitly and are interdispersed with cells on the diagonal, so we
    // must use a map when access the diagonal.

    implicitBondsPtr_ = new edgeList(maxEdges);
    edgeList& implicitBonds = *implicitBondsPtr_;
    faceImplicitBondsPtr_ = new boolList(maxEdges, false);
    boolList& faceImplicitBonds = *faceImplicitBondsPtr_;

    // coupledImplicitBondsPtr_ = new List<edgeList>(bMesh.size());
    // List<edgeList>& coupledImplicitBonds = *coupledImplicitBondsPtr_;
    // coupledShadImplicitBondsPtr_ = new List<edgeList>(bMesh.size());
    // List<edgeList>& coupledShadImplicitBonds = *coupledShadImplicitBondsPtr_;
    // coupledImplicitBondsShadowSizePtr_ = new List<int>(bMesh.size(), 0);
    // List<int>& coupledImplicitBondsShadowSize =
    //     *coupledImplicitBondsShadowSizePtr_;

    label nEdges = 0;

    forAll(c, cellI)
    {
        const labelList& cellCells = cc[cellI];

        // Sort the implicitNei of this cell
        SortableList<label> sortedImplicitNei(implicitNeiOfCell[cellI]);
        forAll(sortedImplicitNei, implicitNeiCellI)
        {
            const label otherCellID = sortedImplicitNei[implicitNeiCellI];
            if (otherCellID > cellI)
            {
                implicitBonds[nEdges] = edge(cellI, otherCellID);

                // Check if the other cell is a near neighbour
                forAll(cellCells, cellCellI)
                {
                    if (otherCellID == cellCells[cellCellI])
                    {
                        faceImplicitBonds[nEdges] = true;
                        break;
                    }
                }
                nEdges++;
            }
        }

        // Add cell-to-boundary-face implicitBonds

        // CHECK: we may consider adding boundary imBonds here to improve
        // matrix bandwidth

        // The implicit bonds designiate the conectivity of the matrix; the
        // boundary faces are additional degrees of freedom after the cell
        // centre degrees of freedom. So we offset the boundary faceID to get
        // its matrix ID. It is like the boundary faces are extra cells at the
        // end of the matrix
        // const label offsetIndex = mesh.nCells() - mesh.nInternalFaces();

        // // Sort the boundary faces implicitNei of this cell
        // SortableList<label> sortedBFaceImplicitNei
        // (implicitBFaceOfCell[cellI]);
        // forAll(sortedBFaceImplicitNei, implicitNeiCellI)
        // {
        //     const label bFaceID = sortedBFaceImplicitNei[implicitNeiCellI];

        //     implicitBonds[nEdges++] = edge(cellI, bFaceID + offsetIndex);
        // }
    }

    forAll(c, cellI)
    {
        // Add cell-to-boundary-face implicitBonds

        // The implicit bonds designiate the conectivity of the matrix; the
        // boundary faces are additional degrees of freedom after the cell
        // centre degrees of freedom. So we offset the boundary faceID to get
        // its matrix ID. It is like the boundary faces are extra cells at the
        // end of the matrix
        //const label offsetIndex = mesh.nCells() - nIntFaces;

        // Sort the boundary faces implicitNei of this cell
        SortableList<label> sortedBFaceImplicitNei(implicitBFaceOfCell[cellI]);
        forAll(sortedBFaceImplicitNei, implicitNeiCellI)
        {
            const label bFaceVarID = sortedBFaceImplicitNei[implicitNeiCellI];

            implicitBonds[nEdges++] = edge(cellI, bFaceVarID);
        }
    }

    // Add boundary-face-to-boundary-face implicitBonds
    forAll(mesh.boundaryMesh(), patchI)
    {
        const word pType = mesh.boundaryMesh()[patchI].type();

        if
        (
            pType != emptyPolyPatch::typeName
         && !mesh.boundaryMesh()[patchI].coupled()
        )
        {
            const label start = mesh.boundaryMesh()[patchI].start();

            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
                const label faceID = start + faceI;
                const label varID = findOldVariableID(faceID);
                const label bVarID = varID - nCells;

                // Sort the neighbour boundary faces of this face
                SortableList<label> sortedImplicitNei
                (
                    implicitNeiOfBFace[bVarID]
                );

                forAll(sortedImplicitNei, implicitNeiFaceI)
                {
                    const label otherFaceVarID =
                        sortedImplicitNei[implicitNeiFaceI];

                    if (otherFaceVarID > varID)
                    {
                        implicitBonds[nEdges++] = edge(varID, otherFaceVarID);
                    }
                }
            }
        }
    }

    // Resize to actual number of implicitBonds
    implicitBonds.resize(nEdges);
    faceImplicitBonds.resize(nEdges);

    // Check the implicit bonds were set correctly

    label minID = implicitBonds[0].start();
    label maxID = implicitBonds[0].start();

    forAll(implicitBonds, imBondI)
    {
        minID = min(minID, implicitBonds[imBondI].start());
        minID = min(minID, implicitBonds[imBondI].end());

        maxID = max(maxID, implicitBonds[imBondI].start());
        maxID = max(maxID, implicitBonds[imBondI].end());
    }

    if (minID < 0 || maxID >= nVariables())
    {
        FatalErrorIn("solidPolyMesh::calcImplicitBonds()")
            << "implicit bonds not set correctly" << nl
            << "minID < 0 or maxID >= nVariables" << nl
            << "minID: " << minID << nl
            << "maxID: " << maxID << nl
            << "nVariables: " << nVariables() << nl
            << "implicitBonds: " << implicitBonds << nl
            << "If you have empty patches, make sure they are at the end!"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "void solidPolyMesh::calcImplicitBonds() const "
            << "finished creating implicit bonds" << endl;
    }
}


void solidPolyMesh::calcFvMeshAddressingMap() const
{
    if (debug)
    {
        Info<< "void solidPolyMesh::calcFvMeshAddressingMap() const "
            << "calculating mapping to fvMesh" << endl;
    }

    // It is an error to attempt to recalculate if the pointer is already set
    if (fvMeshAddressingMapPtr_)
    {
        FatalErrorIn("solidPolyMesh::fvMeshAddressingMap() const")
            << "map already exists" << abort(FatalError);
    }

    const polyMesh& mesh = (*this)();

    fvMeshAddressingMapPtr_ = new labelList(mesh.nInternalFaces(), -1);
    labelList& fvMap = *fvMeshAddressingMapPtr_;

    const boolList& faceImBonds = faceImplicitBonds();

    label fvFaceI = 0;
    forAll(faceImBonds, imBondI)
    {
        if (faceImBonds[imBondI])
        {
            fvMap[fvFaceI] = imBondI;
            fvFaceI++;
        }
    }

    if (min(fvMap) < 0)
    {
        FatalErrorIn
        (
            "void solidPolyMesh::calcFvMeshAddressingMap() const"
        )   << "All faces not set" << abort(FatalError);
    }
}


void solidPolyMesh::calcCellImplicitBonds() const
{
    if (debug)
    {
        Info<< "solidPolyMesh::calcCellImplicitBonds : "
            << "constructing cell implicit bonds arrays"
            << endl;
    }

    // It is an error to attempt to recalculate if the pointer is already set
    if
    (
        cellImplicitBondsPtr_
     || cellBFaceImplicitBondsPtr_
     || bFaceImplicitBondsPtr_
    )
    {
        FatalErrorIn("solidPolyMesh::calcCellImplicitBonds() const")
            << "arrays already exist" << abort(FatalError);
    }

    const polyMesh& mesh = (*this)();

    const label nCells = mesh.nCells();

    cellImplicitBondsPtr_ = new labelListList(nCells);
    labelListList& cellImplicitBonds = *cellImplicitBondsPtr_;

    cellBFaceImplicitBondsPtr_ = new labelListList(nCells);
    labelListList& cellBFaceImplicitBonds = *cellBFaceImplicitBondsPtr_;

    bFaceImplicitBondsPtr_ = new labelListList(this->nVariables() - nCells);
    labelListList& bFaceImplicitBonds = *bFaceImplicitBondsPtr_;

    const labelListList& cpc = cellPointCells();
    const labelListList& cpbf = cellPointBFaces();

    forAll(cellImplicitBonds, cellI)
    {
        cellImplicitBonds[cellI].setSize(cpc[cellI].size(), -1);

        cellBFaceImplicitBonds[cellI].setSize(cpbf[cellI].size(), -1);
    }

    const labelListList& bfpbf = bFacePointBFaces();
    forAll(bFaceImplicitBonds, faceI)
    {
        bFaceImplicitBonds[faceI].setSize(bfpbf[faceI].size(), -1);
    }

    const edgeList& imBonds = implicitBonds();

    // We will count how many bonds have been found for each cell/face as we
    // go
    labelList imBondsFound(mesh.nCells(), 0);
    labelList bFaceImBondsFound(mesh.nCells(), 0);
    labelList bFbFImBondsFound(bFaceImplicitBonds.size(), 0);

    forAll(imBonds, imBondI)
    {
        const edge& e = imBonds[imBondI];
        const label startID = e.start();
        const label endID = e.end();

        if (endID < nCells)
        {
            // Implicit bond connecting two cells

            cellImplicitBonds[startID][imBondsFound[startID]++] = imBondI;
            cellImplicitBonds[endID][imBondsFound[endID]++] = imBondI;

            // Info<< "    " << (imBondsFound[startID] - 1) << nl
            //     << "    " << (imBondsFound[endID] - 1) << nl
            //     << "        to " << imBondI << nl << endl;
        }
        else if (startID < nCells)
        {
            // Implicit bond connecting a cell with a boundary face

            cellBFaceImplicitBonds[startID][bFaceImBondsFound[startID]++] =
                imBondI;

            // label eID = endID - nCells;
            // bFaceImplicitBonds[eID][bFbFImBondsFound[eID]++] = imBondI;

            // Info<< "    " << (bFaceImBondsFound[startID] - 1) << nl
            //     << "    " << (bFaceImBondsFound[endID] - 1) << nl
            //     << "        to " << imBondI << nl << endl;
        }
        else
        {
            // Implicit bond connecting a boundary face with a boundary face

            label sID = startID - nCells;
            label eID = endID - nCells;

            if
            (
                sID >= bFaceImplicitBonds.size()
             || eID >= bFaceImplicitBonds.size()
            )
            {
                FatalErrorIn("solidPolyMesh::calcCellImplicitBonds()")
                    << "Something is wrong: this may happen for "
                    << "very coarse meshes" << endl
                    << "sID: " << sID << nl
                    << "eID: " << eID << nl
                    << "nCells: " << nCells << nl
                    << "startID: " << startID << nl
                    << "endID: " << endID << nl
                    << "bFaceImplicitBonds.size(): "
                    << bFaceImplicitBonds.size() << nl
                    << abort(FatalError);
            }

            bFaceImplicitBonds[sID][bFbFImBondsFound[sID]++] = imBondI;
            bFaceImplicitBonds[eID][bFbFImBondsFound[eID]++] = imBondI;

            // Info<< "bFaceImplicitBonds[" << sID << "]["
            //     << bFbFImBondsFound[sID] << "]: "
            //     << bFaceImplicitBonds[sID][bFbFImBondsFound[sID]] << nl
            //     << "bFaceImplicitBonds[" << eID << "]["
            //     << bFbFImBondsFound[eID] << "]: "
            //     << bFaceImplicitBonds[eID][bFbFImBondsFound[eID]] << nl
            //     << "    set to " << imBondI << nl
            //     << endl;
        }
    }

    // Check all edges were set

    label minImBondID = 0;

    forAll(cellImplicitBonds, cellI)
    {
        minImBondID = min(min(cellImplicitBonds[cellI]), minImBondID);
    }

    if (minImBondID < 0)
    {
        FatalErrorIn
        (
            "void solidPolyMesh::calcCellImplicitBonds() const"
        )   << "Not all cell implicit bonds set" << abort(FatalError);
    }

    minImBondID = 0;

    forAll(cellBFaceImplicitBonds, cellI)
    {
        if (cellBFaceImplicitBonds[cellI].size())
        {
            minImBondID = min(min(cellBFaceImplicitBonds[cellI]), minImBondID);
        }
    }

    if (minImBondID < 0)
    {
        FatalErrorIn
        (
            "void solidPolyMesh::calcCellImplicitBonds() const"
        )   << "Not all cell boundary face implicit bonds set"
            //<< "cellBFaceImplicitBonds: " << cellBFaceImplicitBonds
            << abort(FatalError);
    }

    minImBondID = 0;

    forAll(bFaceImplicitBonds, faceI)
    {
        if (bFaceImplicitBonds[faceI].size())
        {
            minImBondID = min(min(bFaceImplicitBonds[faceI]), minImBondID);
        }
    }

    if (minImBondID < 0)
    {
        FatalErrorIn
        (
            "void solidPolyMesh::calcCellImplicitBonds() const"
        )   << "Not all boundary-face-to-boundary-face implicit bonds set"
            << nl << "bFaceImplicitBonds: " << bFaceImplicitBonds
            //<< nl << "bfpbf: " << bfpbf
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "solidPolyMesh::calcCellImplicitBonds : "
            << "finished constructing cell implicit bonds arrays"
            << endl;
    }
}


void solidPolyMesh::makeGlobalCoeffs
(
    const newLeastSquaresVolPointInterpolation& volToPointInterp
) const
{
    if (debug)
    {
        Info<< "void solidPolyMesh::makeGlobalCoeffs() const "
            << "creating global coefficient fields" << endl;
    }

    // Return if the pointer has already been set
    if
    (
        pointProcFacesCoeffsPtr_
     || pointProcBndFacesCoeffsPtr_
     || pointProcCellsCoeffsPtr_
     || gPtNgbProcBndFaceCoeffsPtr_
     || gPtNgbProcCellCoeffsPtr_
    )
    {
        // WarningIn("solidPolyMesh::makeGlobalCoeffs() const")
        //     << "globalCoeffs already made"
        //     << endl;

        return;
    }

    const polyMesh& mesh = (*this)();
    const label nVar = this->nVariables();

    pointProcFacesCoeffsPtr_ = new List< Map<tensorField> >(nVar);
    List< Map<tensorField> >& pointProcFacesCoeffs = *pointProcFacesCoeffsPtr_;

    pointProcBndFacesCoeffsPtr_ = new List< Map<tensorField> >(nVar);
    List< Map<tensorField> >& pointProcBndFacesCoeffs =
        *pointProcBndFacesCoeffsPtr_;

    pointProcCellsCoeffsPtr_ = new List< Map<tensorField> >(nVar);
    List< Map<tensorField> >& pointProcCellsCoeffs = *pointProcCellsCoeffsPtr_;

    gPtNgbProcBndFaceCoeffsPtr_ = new List< Map<tensorField> >(nVar);
    List< Map<tensorField> >& gPtNgbProcBndFaceCoeffs =
        *gPtNgbProcBndFaceCoeffsPtr_;

    gPtNgbProcCellCoeffsPtr_ = new List< Map<tensorField> >(nVar);
    List< Map<tensorField> >& gPtNgbProcCellCoeffs = *gPtNgbProcCellCoeffsPtr_;

    const label nCells = mesh.nCells();
    const label nInternalFaces = mesh.nInternalFaces();
    //const cellList& cells = mesh.cells();
    const labelListList& cellPoints = mesh.cellPoints();
    const faceList& faces = mesh.faces();
    const labelListList& pointProcFaces = volToPointInterp.pointProcFaces();

    // We need the working field to create the global point fields
    const volVectorField* vfPtr = NULL;
    if (mesh.foundObject<volVectorField>("DD"))
    {
        vfPtr = &mesh.lookupObject<volVectorField>("DD");
    }
    else if (mesh.foundObject<volVectorField>("D"))
    {
        vfPtr = &mesh.lookupObject<volVectorField>("D");
    }
    else
    {
        FatalErrorIn("solidPolyMesh::makeGlobalCoeffs()")
            << "Currently the working field must be called D or DD"
            << abort(FatalError);
    }
    const volVectorField& vf = *vfPtr;

    Map<Field<vector> > gPtNgbProcBndFaceFieldData;
    volToPointInterp.globalPointNgbProcBndFaceFieldData
    (
        vf, gPtNgbProcBndFaceFieldData
    );
    Map<Field<vector> > gPtNgbProcCellFieldData;
    volToPointInterp.globalPointNgbProcCellFieldData
    (
        vf, gPtNgbProcCellFieldData
    );
    const Map< List<labelPair> >& pointProcCells =
        volToPointInterp.pointProcCells();
    const List< List<labelPair> >& pointProcBndFaces =
        volToPointInterp.pointProcBndFaces();

    for (label varI = 0; varI < nVar; varI++)
    {
        if (varI < nCells)
        {
            const labelList& cp = cellPoints[varI];

            forAll(cp, cellPointI)
            {
                const label pointID = cp[cellPointI];

                if (pointProcFaces[pointID].size())
                {
                    pointProcFacesCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            pointProcFaces[pointID].size(), tensor::zero
                        )
                    );
                }

                if (gPtNgbProcBndFaceFieldData.found(pointID))
                {
                    gPtNgbProcBndFaceCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            gPtNgbProcBndFaceFieldData[pointID].size(),
                            tensor::zero
                        )
                    );
                }

                if (gPtNgbProcCellFieldData.found(pointID))
                {
                    gPtNgbProcCellCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            gPtNgbProcCellFieldData[pointID].size(),
                            tensor::zero
                        )
                    );
                }

                if (pointProcCells.found(pointID))
                {
                    pointProcCellsCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            pointProcCells[pointID].size(), tensor::zero
                        )
                    );
                }

                if (pointProcBndFaces[pointID].size())
                {
                    pointProcBndFacesCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            pointProcBndFaces[pointID].size(), tensor::zero
                        )
                    );
                }
            }
        }
        else
        {
            const label faceI = varI - nCells + nInternalFaces;

            // Check faceI is correct
            // if (varI != findOldVariableID(faceI))
            // {
            //     FatalErrorIn("makeGlobalCoeffsNew")
            //         << "faceI is incorrectly calculated from varI"
            //         << abort(FatalError);
            // }

            const labelList& fp = faces[faceI];

            forAll(fp, facePointI)
            {
                const label pointID = fp[facePointI];

                if (pointProcFaces[pointID].size())
                {
                    pointProcFacesCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            pointProcFaces[pointID].size(), tensor::zero
                        )
                    );
                }

                if (gPtNgbProcBndFaceFieldData.found(pointID))
                {
                    gPtNgbProcBndFaceCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            gPtNgbProcBndFaceFieldData[pointID].size(),
                            tensor::zero
                        )
                    );
                }

                if (gPtNgbProcCellFieldData.found(pointID))
                {
                    gPtNgbProcCellCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            gPtNgbProcCellFieldData[pointID].size(),
                            tensor::zero
                        )
                    );
                }

                if (pointProcCells.found(pointID))
                {
                    pointProcCellsCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            pointProcCells[pointID].size(), tensor::zero
                        )
                    );
                }

                if (pointProcBndFaces[pointID].size())
                {
                    pointProcBndFacesCoeffs[varI].insert
                    (
                        pointID,
                        tensorField
                        (
                            pointProcBndFaces[pointID].size(), tensor::zero
                        )
                    );
                }
            }
        }
    }
}


void solidPolyMesh::calcPointFixed(const volVectorField& D) const
{
    if (pointFixedComponentPtr_ || pointFixedDirectionPtr_)
    {
        FatalErrorIn("void solidPolyMesh::calcPointFixed() const")
            << "pointers already set" << abort(FatalError);
    }

    const fvMesh& mesh = D.mesh();

    pointFixedComponentPtr_ = new Map<vector>(0.1*mesh.nPoints());
    pointFixedDirectionPtr_ = new Map<symmTensor>(0.1*mesh.nPoints());

    Map<vector>& pointFixedComp = *pointFixedComponentPtr_;
    Map<symmTensor>& pointFixedDir = *pointFixedDirectionPtr_;

    // We need the point field of the working vol field
    const pointVectorField* pfPtr = NULL;
    if (mesh.foundObject<pointVectorField>("pointDD"))
    {
        pfPtr = &mesh.lookupObject<pointVectorField>("pointDD");
    }
    else if (mesh.foundObject<pointVectorField>("pointD"))
    {
        pfPtr = &mesh.lookupObject<pointVectorField>("pointD");
    }
    else
    {
        FatalErrorIn("solidPolyMesh::calcPointFixed()")
            << "Currently the working field must be called D or DD"
            << " and there should be a corresponding pointD or pointDD field"
            << abort(FatalError);
    }
    const pointVectorField& pf = *pfPtr;

    // Force pointD to correct boundary condition e.g. so fixedValue boundaries
    // will set the point internalField values
    const_cast<pointVectorField&>(pf).correctBoundaryConditions();

    forAll(D.boundaryField(), patchI)
    {
        const word pType = D.boundaryField()[patchI].type();

        if (pType == blockFixedDisplacementFvPatchVectorField::typeName)
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // The points are fixed in every direction
            const symmTensor fixedDir = symmTensor(I);

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                // Add point and overwrite if is is already there
                pointFixedComp.set(pointID, pf[pointID]);
                pointFixedDir.set(pointID, fixedDir);
            }
        }
        else if
        (
            pType == blockFixedDisplacementZeroShearFvPatchVectorField::typeName
         || pType == "solidSymmetry"
         || pType == "symmetry"
        )
        {
            // WarningIn("solidPolyMesh.C")
            //     << "Applying point fixed components on patch "
            //     << mesh.boundaryMesh()[patchI].name() << endl;

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            const vectorField& pointNormals =
                mesh.boundaryMesh()[patchI].pointNormals();

            const symmTensor symmTensorI = symmTensor(I);

            // The points are fixed in the normal direction

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                if (!pointFixedComp.found(pointID))
                {
                    pointFixedComp.insert(pointID, pf[pointID]);
                    pointFixedDir.insert(pointID, sqr(pointNormals[pI]));
                }
                else
                {
                    // Previous direction
                    const symmTensor prevDir = pointFixedDir[pointID];

                    // Do nothing if the point is already fully fixed
                    if (prevDir != symmTensorI)
                    {
                        // Overwrite the new direction and component but keep
                        // any previously fixed components
                        const vector prevComp = pointFixedComp[pointID];
                        const symmTensor newDir = sqr(pointNormals[pI]);

                        pointFixedComp.set
                        (
                            pointID,
                            (newDir & pf[pointID]) + ((I - newDir) & prevComp)
                        );

                        pointFixedDir.set
                        (
                            pointID,
                            symm(newDir + ((I - newDir) & prevDir))
                        );

                        if (debug)
                        {
                            Info<< pointID << " new dir "
                                << newDir + ((I - newDir) & prevDir)
                                << " new comp "
                                << (newDir & pf[pointID])
                                 + ((I - newDir) & prevComp)
                                << endl;
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidPolyMesh::solidPolyMesh(fvMesh& fMesh)
:
    GeoMesh<polyMesh>(fMesh),
    IOdictionary
    (
        IOobject
        (
            "solidPolyMesh",
            fMesh.time().timeName(),
            fMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    fMesh_(fMesh),
    boundary_(*this, fMesh.boundaryMesh()),
    lduPtr_(NULL),
    implicitBondsPtr_(NULL),
    faceImplicitBondsPtr_(NULL),
    fvMeshAddressingMapPtr_(NULL),
    cellPointCellsPtr_(NULL),
    cellPointBFacesPtr_(NULL),
    bFacePointBFacesPtr_(NULL),
    cellImplicitBondsPtr_(NULL),
    cellBFaceImplicitBondsPtr_(NULL),
    bFaceImplicitBondsPtr_(NULL),
    pointProcFacesCoeffsPtr_(NULL),
    pointProcBndFacesCoeffsPtr_(NULL),
    pointProcCellsCoeffsPtr_(NULL),
    gPtNgbProcBndFaceCoeffsPtr_(NULL),
    gPtNgbProcCellCoeffsPtr_(NULL),
    pointFixedComponentPtr_(NULL),
    pointFixedDirectionPtr_(NULL),
    comm_(Pstream::worldComm)
{
    if (debug)
    {
        Info<< "solidPolyMesh::solidPolyMesh(const fvMesh&) : "
            << "Creating solidPolyMesh" << endl;
    }

    // WIP - forcing updateGlobalFields to be called by
    // blockGlobalSolidFvPatchVectorField coupled patch but I get a MPI receive
    // error in newLeastSquaresInterpolate... I'm not sure why because it works
    // when I called updateGlobalFields directly from within
    // BlockLduVectorMatrixATmul.C...
    //addGlobalPatch();
    // Check if a globalPatch exsits
    // if (Pstream::parRun())
    // {
    //     bool globalPatchFound = false;

    //     forAll(fMesh_.boundary(), patchI)
    //     {
    //         if
    //         (
    //         fMesh_.boundary()[patchI].type() == blockGlobalFvPatch::typeName
    //         )
    //         {
    //             globalPatchFound = true;
    //             break;
    //         }
    //     }

    //     if (!globalPatchFound)
    //     {
    //         FatalErrorIn("solidPolyMesh::solidPolyMesh")
    //             << "blockGlobalSolid patch not found:" << nl
    //         << "Please add a patch called blockGlobalSolid to the end of "
    //             << "the boundary list"
    //             << " with zero faces, before decomposing the case," << nl
    //             << " with 'type blockGlobalSolid;'" << nl
    //             << "Also add a blockGlobalSolid patch to the D field"
    //             << abort(FatalError);
    //     }
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

solidPolyMesh::~solidPolyMesh()
{
    if (debug)
    {
        Info<< "solidPolyMesh::~solidPolyMesh() : "
            << "Deleting solidPolyMesh" << endl;
    }

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const lduAddressing& solidPolyMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new solidPolyMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


const edgeList& solidPolyMesh::implicitBonds() const
{
    if (!implicitBondsPtr_)
    {
        calcImplicitBonds();
    }

    return *implicitBondsPtr_;
}


const boolList& solidPolyMesh::faceImplicitBonds() const
{
    if (!faceImplicitBondsPtr_)
    {
        calcImplicitBonds();
    }

    return *faceImplicitBondsPtr_;
}


const labelList& solidPolyMesh::fvMeshAddressingMap() const
{
    if (!fvMeshAddressingMapPtr_)
    {
        calcFvMeshAddressingMap();
    }

    return *fvMeshAddressingMapPtr_;
}


const labelListList& solidPolyMesh::cellPointCells() const
{
    if (!cellPointCellsPtr_)
    {
        calcImplicitBonds();
    }

    return *cellPointCellsPtr_;
}


const labelListList& solidPolyMesh::cellPointBFaces() const
{
    if (!cellPointBFacesPtr_)
    {
        calcImplicitBonds();
    }

    return *cellPointBFacesPtr_;
}


const labelListList& solidPolyMesh::bFacePointBFaces() const
{
    if (!bFacePointBFacesPtr_)
    {
        calcImplicitBonds();
    }

    return *bFacePointBFacesPtr_;
}

const labelListList& solidPolyMesh::cellImplicitBonds() const
{
    if (!cellImplicitBondsPtr_)
    {
        calcCellImplicitBonds();
    }

    return *cellImplicitBondsPtr_;
}


const labelListList& solidPolyMesh::cellBFaceImplicitBonds() const
{
    if (!cellBFaceImplicitBondsPtr_)
    {
        calcCellImplicitBonds();
    }

    return *cellBFaceImplicitBondsPtr_;
}


const labelListList& solidPolyMesh::bFaceImplicitBonds() const
{
    if (!bFaceImplicitBondsPtr_)
    {
        calcCellImplicitBonds();
    }

    return *bFaceImplicitBondsPtr_;
}


label solidPolyMesh::nVariables() const
{
    // The number of solution variables =
    // nCells + nBoundaryFaces - emptyBoundaryFaces - procBoundaryFaces

    label emptyFaces = 0;
    label procFaces = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const word& patchType = mesh_.boundaryMesh()[patchI].type();

        if (patchType == emptyPolyPatch::typeName)
        {
            emptyFaces += mesh_.boundaryMesh()[patchI].size();
        }
        else if (patchType == processorPolyPatch::typeName)
        {
            procFaces += mesh_.boundaryMesh()[patchI].size();
        }
    }

    return
    (
        mesh_.nCells() + mesh_.nFaces() - mesh_.nInternalFaces() - emptyFaces
        - procFaces
    );
}


label solidPolyMesh::findCellCellImplicitBond
(
    const label firstCellID,
    const label secondCellID
) const
{
    // Check which cell has fewer implicit bonds connected to it
    label curCellID = firstCellID;
    label otherCellID = secondCellID;
    if
    (
        cellImplicitBonds()[firstCellID].size()
      > cellImplicitBonds()[secondCellID].size()
    )
    {
        curCellID = secondCellID;
        otherCellID = firstCellID;
    }
    const labelList& curCellImplicitBonds = cellImplicitBonds()[curCellID];

    // We will loop through all implicit bonds of curCell until we find the
    // edge which connects the curCell and otherCell

    label foundEdgeID = -1;
    const edgeList& imBonds = implicitBonds();

    forAll(curCellImplicitBonds, imBondI)
    {
        const edge curImBond = imBonds[curCellImplicitBonds[imBondI]];
        if
        (
            curImBond.start() == otherCellID || curImBond.end() == otherCellID
        )
        {
            foundEdgeID = imBondI;
            break;
        }
    }

    if (foundEdgeID == -1)
    {
        return -1;
    }

    return curCellImplicitBonds[foundEdgeID];
}


label solidPolyMesh::findCellFaceImplicitBond
(
    const label cellID, const label bFaceID
) const
{
    const labelList& curCellBFaceImplicitBonds =
        cellBFaceImplicitBonds()[cellID];

    const polyMesh& mesh = (*this)();

    // We will loop through all boundary face implicit bonds of current cell
    // until we find the edge which connects the current cell and boundary face

    label foundEdgeID = -1;
    const edgeList& imBonds = implicitBonds();

    // Calculate matrix index of boundary face
    // CHECK: we may reconsider how we calculate the index for boundary faces
    // due to emptyFaces; maybe store the pre empty faces
    label preEmptyFaces = 0;
    label patchID = mesh.boundaryMesh().whichPatch(bFaceID);
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == emptyPolyPatch::typeName)
        {
            if (patchI < patchID)
            {
                preEmptyFaces += mesh.boundaryMesh()[patchI].size();
            }
        }
    }
    const label offset = mesh.nCells() - mesh.nInternalFaces() - preEmptyFaces;
    const label bFaceImBondID = bFaceID + offset;

    forAll(curCellBFaceImplicitBonds, imBondI)
    {
        const edge curImBond = imBonds[curCellBFaceImplicitBonds[imBondI]];

        if (curImBond.end() == bFaceImBondID)
        {
            foundEdgeID = imBondI;
            break;
        }
    }

    if (foundEdgeID == -1)
    {
        Warning
            << "cell-to-bFace not found" << nl
            << "curCellBFaceImplicitBonds[" << cellID << "] " << nl
            << "    C " << mesh.cellCentres()[cellID] << nl
            << "    Cf " << mesh.faceCentres()[bFaceID] << nl
            << "    bFaceID " << bFaceID << nl
            << "    bFaceImBondID " << bFaceImBondID << nl
            << "    cellImBonds " << curCellBFaceImplicitBonds
            << endl;

        return -1;
    }

    return curCellBFaceImplicitBonds[foundEdgeID];
}


label solidPolyMesh::findFaceFaceImplicitBond
(
    const label firstFaceID,
    const label secondFaceID
) const
{
    const polyMesh& mesh = (*this)();
    const label nIntFaces = mesh.nInternalFaces();

    // Calculate matrix index of boundary face
    // CHECK: we may reconsider how we calculate the index for boundary faces
    // due to emptyFaces; maybe store the pre empty faces
    label firstPreEmptyFaces = 0;
    label secondPreEmptyFaces = 0;
    label firstPatchID = mesh.boundaryMesh().whichPatch(firstFaceID);
    label secondPatchID = mesh.boundaryMesh().whichPatch(secondFaceID);
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == emptyPolyPatch::typeName)
        {
            if (patchI < firstPatchID)
            {
                firstPreEmptyFaces += mesh.boundaryMesh()[patchI].size();
            }

            if (patchI < secondPatchID)
            {
                secondPreEmptyFaces += mesh.boundaryMesh()[patchI].size();
            }
        }
    }
    const label firstOffset = mesh.nCells() - nIntFaces - firstPreEmptyFaces;
    const label secondOffset = mesh.nCells() - nIntFaces - secondPreEmptyFaces;
    const label firstID = firstFaceID + firstOffset;
    const label secondID = secondFaceID + secondOffset;
    const label firstBFaceID = firstFaceID - nIntFaces - firstPreEmptyFaces;
    const label secondBFaceID = secondFaceID - nIntFaces - secondPreEmptyFaces;

    // We will loop through all boundary implicit bonds of curFace until we
    // find the edge which connects the curFace and otherFace

    const labelList& curFaceImplicitBonds = bFaceImplicitBonds()[firstBFaceID];

    label foundEdgeID = -1;
    const edgeList& imBonds = implicitBonds();

    forAll(curFaceImplicitBonds, imBondI)
    {
        const edge curImBond = imBonds[curFaceImplicitBonds[imBondI]];
        if
        (
            curImBond.start() == secondID || curImBond.end() == secondID
        )
        {
            foundEdgeID = imBondI;
            break;
        }
    }

    if (foundEdgeID == -1)
    {
        Warning
            << "bFace-to-bFace imBond not found" << nl
            << "curFaceID " << firstID << nl
            << "    Cf " << mesh.faceCentres()[firstFaceID] << nl;

        forAll(bFaceImplicitBonds()[firstID], ie)
        {
            Info<< "    imBonds "
                << imBonds[bFaceImplicitBonds()[firstBFaceID][ie]]
                << endl;
        }

        Info<< "otherFaceID " << secondID << nl
            << "    Cf " << mesh.faceCentres()[secondFaceID] << nl;

        forAll(bFaceImplicitBonds()[secondID], ie)
        {
            Info<< "    imBonds "
                << imBonds[bFaceImplicitBonds()[secondBFaceID][ie]]
                << endl;
        }

        return -1;
    }

    return curFaceImplicitBonds[foundEdgeID];
}


label solidPolyMesh::findOldVariableID(const label bFaceID) const
{
    const polyMesh& mesh = (*this)();
    const label nIntFaces = mesh.nInternalFaces();
    const label nCells = mesh.nCells();

    label preEmptyFaces = 0;

    forAll(mesh.boundaryMesh(), patchI)
    {
        const word pType = mesh.boundaryMesh()[patchI].type();

        if (pType == emptyPolyPatch::typeName)
        {
            preEmptyFaces += mesh.boundaryMesh()[patchI].size();
        }
        else if (pType != processorPolyPatch::typeName)
        {
            return nCells + bFaceID - nIntFaces - preEmptyFaces;
        }
    }

    WarningIn("label findOldVariableID(const label bFaceID) const")
        << "bFaceID is not on a standard boundary!" << endl;

    return -1;
}


void solidPolyMesh::addFvMatrix
(
    BlockLduMatrix<vector>& blockM,
    vectorField& blockB,
    const fvVectorMatrix& fvM,
    const bool flipSign
) const
{
    if (!fvM.hasDiag() && !fvM.hasUpper())
    {
        // Copy source constributions
        const vectorField& fvSource = fvM.source();
        if (flipSign)
        {
            forAll(fvSource, cellI)
            {
                // Source contributions
                blockB[cellI] -= fvSource[cellI];
            }
        }
        else
        {
            forAll(fvSource, cellI)
            {
                // Source contributions
                blockB[cellI] += fvSource[cellI];
            }
        }

        return;
    }

    if (!fvM.diagonal())
    {
        FatalErrorIn
        (
            "void solidPolyMesh::addFvMatrix"
            "("
            "    BlockLduMatrix<vector>& blockM,"
            "    vectorField& blockB,"
            "    const fvVectorMatrix& fvM,"
            "    const bool flipSign"
            ") const"
        )   << "only implemented for a diagonal fvMatrix!" << abort(FatalError);
    }

    const scalarField& fvD = fvM.diag();
    const vectorField& fvSource = fvM.source();
    Field<tensor>& blockD = blockM.diag().asSquare();

    if (flipSign)
    {
        forAll(fvD, cellI)
        {
            // Diagonal contributions
            blockD[cellI] -= tensor(fvD[cellI]*I);

            // Source contributions
            blockB[cellI] -= fvSource[cellI];
        }
    }
    else
    {
        forAll(fvD, cellI)
        {
            // Diagonal contributions
            blockD[cellI] += tensor(fvD[cellI]*I);

            // Source contributions
            blockB[cellI] += fvSource[cellI];
        }
    }
}


void solidPolyMesh::addFvSource
(
    vectorField& blockB,
    const volVectorField& fvB
) const
{
    const vectorField& fvBI = fvB.internalField();

    forAll(fvBI, cellI)
    {
        blockB[cellI] += fvBI[cellI];
    }
}


void solidPolyMesh::insertBoundaryConditions
(
    BlockLduMatrix<vector>& blockM,
    vectorField& blockB,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const volVectorField& D
) const
{
    if (debug)
    {
        Info<< "Adding boundary condition terms to matrix" << endl;
    }

    // Const reference to fvMesh
    const fvMesh& mesh = D.mesh();

    // Insert discretised boundary condition equations

    forAll(mesh.boundary(), patchI)
    {
        const word pType = D.boundaryField()[patchI].type();

        if (pType == "solidSymmetry" || pType == "symmetry")
        {
            FatalErrorIn("void solidPolyMesh::insertBoundaryConditions")
                << "symmetry disabled: use instead "
                << blockFixedDisplacementZeroShearFvPatchVectorField::typeName
                << abort(FatalError);
        }
        else if
        (
            pType != emptyFvPatchVectorField::typeName
         && !mesh.boundaryMesh()[patchI].coupled()
        )
        {
            const blockFvPatchVectorField& blockPatch =
                refCast<const blockFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            blockPatch.insertBlockCoeffs
            (
                *this,
                muf,
                lambdaf,
                D,
                blockB,
                blockM
            );
        }
    }
}


List< Map<tensorField> >& solidPolyMesh::pointProcFacesCoeffs() const
{
    if (!pointProcFacesCoeffsPtr_)
    {
        FatalErrorIn
        (
            "const List< Map<tensorField> >& solidPolyMesh::"
            "pointProcFacesCoeffs() const"
        )   << "pointer not set" << abort(FatalError);
    }

    return *pointProcFacesCoeffsPtr_;
}


List< Map<tensorField> >& solidPolyMesh::pointProcBndFacesCoeffs() const
{
    if (!pointProcBndFacesCoeffsPtr_)
    {
        FatalErrorIn
        (
            "const List< Map<tensorField> >& solidPolyMesh::"
            "pointProcBndFacesCoeffs() const"
        )   << "pointer not set" << abort(FatalError);
    }

    return *pointProcBndFacesCoeffsPtr_;
}


List< Map<tensorField> >& solidPolyMesh::pointProcCellsCoeffs() const
{
    if (!pointProcCellsCoeffsPtr_)
    {
        FatalErrorIn
        (
            "const List< Map<tensorField> >& solidPolyMesh::"
            "pointProcCellsCoeffs() const"
        )   << "pointer not set" << abort(FatalError);
    }

    return *pointProcCellsCoeffsPtr_;
}


List< Map<tensorField> >& solidPolyMesh::gPtNgbProcBndFaceCoeffs() const
{
    if (!gPtNgbProcBndFaceCoeffsPtr_)
    {
        FatalErrorIn
        (
            "const List< Map<tensorField> >& solidPolyMesh::"
            "gPtNgbProcBndFaceCoeffs() const"
        )   << "pointer not set" << abort(FatalError);
    }

    return *gPtNgbProcBndFaceCoeffsPtr_;
}


List< Map<tensorField> >& solidPolyMesh::gPtNgbProcCellCoeffs() const
{
    if (!gPtNgbProcCellCoeffsPtr_)
    {
        FatalErrorIn
        (
            "const List< Map<tensorField> >& solidPolyMesh::"
            "gPtNgbProcCellCoeffs() const"
        )   << "pointer not set" << abort(FatalError);
    }

    return *gPtNgbProcCellCoeffsPtr_;
}


void solidPolyMesh::updateGlobalFields
(
    vectorField& Ax,
    const vectorField& x
) const
{
    if (Pstream::parRun())
    {
        const polyMesh& mesh = (*this)();

        const newLeastSquaresVolPointInterpolation& volToPointInterp =
            this->volToPointInterp();

        const labelListList& pointProcFaces = volToPointInterp.pointProcFaces();
        const Map< List<labelPair> >& pointProcCells =
            volToPointInterp.pointProcCells();
        const List< List<labelPair> >& pointProcBndFaces =
            volToPointInterp.pointProcBndFaces();

        // We need the working field to create the global point fields
        // We should figure out a better way where we don't need to const_cast!
        const volVectorField* vfPtr = NULL;
        if (mesh.foundObject<volVectorField>("DD"))
        {
            vfPtr = &mesh.lookupObject<volVectorField>("DD");
        }
        else if (mesh.foundObject<volVectorField>("D"))
        {
            vfPtr = &mesh.lookupObject<volVectorField>("D");
        }
        else
        {
            FatalErrorIn("solidPolyMesh::makeGlobalCoeffs()")
                << "Currently the working field must be called D or DD"
                << abort(FatalError);
        }
        volVectorField& vf = const_cast<volVectorField&>(*vfPtr);


        // Update vf
        // The volToPoint parallel functions need vf to be up-to-date, so we
        // will copy the solution vector into the D field
        copySolutionVector(x, vf);

        Map<Field<vector> > gPtNgbProcBndFaceFieldData;
        volToPointInterp.globalPointNgbProcBndFaceFieldData
        (
            vf, gPtNgbProcBndFaceFieldData
        );

        Map<Field<vector> > gPtNgbProcCellFieldData;
        volToPointInterp.globalPointNgbProcCellFieldData
        (
            vf, gPtNgbProcCellFieldData
        );

        const FieldField<Field, vector> procCellVfI =
            volToPointInterp.procCellsFieldData(vf.internalField());
        const FieldField<Field, vector> procBndFaceVf =
            volToPointInterp.procBndFacesFieldData(vf);

        const List< Map<tensorField> >& pointProcFacesCoeffs =
            this->pointProcFacesCoeffs();
        const List< Map<tensorField> >& pointProcBndFacesCoeffs =
            this->pointProcBndFacesCoeffs();
        const List< Map<tensorField> >& pointProcCellsCoeffs =
            this->pointProcCellsCoeffs();
        const List< Map<tensorField> >& gPtNgbProcBndFaceCoeffs =
            this->gPtNgbProcBndFaceCoeffs();
        const List< Map<tensorField> >& gPtNgbProcCellCoeffs =
            this->gPtNgbProcCellCoeffs();

        const label nVar = x.size();

        for (label varI = 0; varI < nVar; varI++)
        {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~              proc boundary faces                 ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            // Cells on other procs which are direclty accessible through proc
            // faces

            if (pointProcFacesCoeffs[varI].size())
            {
                // List of point IDs
                const labelList pointIDs = pointProcFacesCoeffs[varI].toc();

                // List of coeffs
                const Map<tensorField>& coeffsMap = pointProcFacesCoeffs[varI];

                forAll(pointIDs, pI)
                {
                    const label pointID = pointIDs[pI];
                    const tensorField& coeffs = coeffsMap[pointID];
                    const labelList& interpProcFaces =
                        pointProcFaces[pointID];

                    forAll(coeffs, cI)
                    {
                        const label faceID = interpProcFaces[cI];
                        const label patchID =
                            mesh.boundaryMesh().whichPatch(faceID);
                        const label localFaceID =
                            faceID - mesh.boundaryMesh()[patchID].start();

                        // Add contribution
                        Ax[varI] +=
                            coeffs[cI]
                            & vf.boundaryField()[patchID][localFaceID];
                    }
                }
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~             global boundary faces                ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            if (gPtNgbProcBndFaceCoeffs[varI].size())
            {
                // List of point IDs
                const labelList pointIDs = gPtNgbProcBndFaceCoeffs[varI].toc();

                // List of coeffs
                const Map<tensorField>& coeffsMap =
                    gPtNgbProcBndFaceCoeffs[varI];

                forAll(pointIDs, pI)
                {
                    const label pointID = pointIDs[pI];
                    const tensorField& coeffs = coeffsMap[pointID];

                    // Assemble values
                    Field<vector> glNgbProcBndFaceData =
                        gPtNgbProcBndFaceFieldData[pointID];

                    forAll(coeffs, cI)
                    {
                        // Add contribution
                        Ax[varI] += coeffs[cI] & glNgbProcBndFaceData[cI];
                    }
                }
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~                 global cells                     ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            if (gPtNgbProcCellCoeffs[varI].size())
            {
                // List of point IDs
                const labelList pointIDs = gPtNgbProcCellCoeffs[varI].toc();

                // List of coeffs
                const Map<tensorField>& coeffsMap = gPtNgbProcCellCoeffs[varI];

                forAll(pointIDs, pI)
                {
                    const label pointID = pointIDs[pI];
                    const tensorField& coeffs = coeffsMap[pointID];

                    // Assemble values
                    Field<vector> glNgbProcCellData =
                        gPtNgbProcCellFieldData[pointID];

                    forAll(coeffs, cI)
                    {
                        // Add contribution
                        Ax[varI] += coeffs[cI] & glNgbProcCellData[cI];
                    }
                }
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~                  proc cells                      ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            if (pointProcCellsCoeffs[varI].size())
            {
                // List of point IDs
                const labelList pointIDs = pointProcCellsCoeffs[varI].toc();

                // List of coeffs
                const Map<tensorField>& coeffsMap = pointProcCellsCoeffs[varI];

                forAll(pointIDs, pI)
                {
                    const label pointID = pointIDs[pI];
                    const tensorField& coeffs = coeffsMap[pointID];

                    // Assemble values
                    Field<vector> interpNgbProcCellData(0);
                    const List<labelPair>& pc = pointProcCells[pointID];
                    interpNgbProcCellData.setSize(pc.size());

                    forAll(pc, cI)
                    {
                        interpNgbProcCellData[cI] =
                            procCellVfI
                            [
                                pc[cI].first()
                            ]
                            [
                                pc[cI].second()
                            ];
                    }

                    forAll(coeffs, cI)
                    {
                        // Add contribution
                        Ax[varI] += coeffs[cI] & interpNgbProcCellData[cI];
                    }
                }
            }

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~         boundary faces on other procs            ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            if (pointProcBndFacesCoeffs[varI].size())
            {
                // List of point IDs
                const labelList pointIDs = pointProcBndFacesCoeffs[varI].toc();

                // List of coeffs
                const Map<tensorField>& coeffsMap =
                    pointProcBndFacesCoeffs[varI];

                forAll(pointIDs, pI)
                {
                    const label pointID = pointIDs[pI];
                    const tensorField& coeffs = coeffsMap[pointID];

                    // Assemble values
                    Field<vector> interpNgbProcBndFaceData(0);
                    const List<labelPair>& pf = pointProcBndFaces[pointID];
                    interpNgbProcBndFaceData.setSize(pf.size());

                    forAll(pf, fI)
                    {
                        interpNgbProcBndFaceData[fI] =
                            procBndFaceVf
                            [
                                pf[fI].first()
                            ]
                            [
                                pf[fI].second()
                            ];
                    }

                    forAll(coeffs, cI)
                    {
                        // Add contribution
                        Ax[varI] += coeffs[cI] & interpNgbProcBndFaceData[cI];
                    }
                }
            }
        }
    } // if parRun
}


void solidPolyMesh::copySolutionVector
(
    const vectorField& x,
    volVectorField& vf
) const
{
    // Internal field
    vectorField& vfI = vf.internalField();
    forAll(vfI, cellI)
    {
        vfI[cellI] = x[cellI];
    }

    // Boundary field
    label preEmptyFaces = 0;
    const polyMesh& mesh = (*this)();
    const label nCells = mesh.nCells();
    const label nIntFaces = mesh.nInternalFaces();
    forAll(vf.boundaryField(), patchI)
    {
        const word patchType = mesh.boundaryMesh()[patchI].type();

        if (patchType == emptyPolyPatch::typeName)
        {
            preEmptyFaces += mesh.boundaryMesh()[patchI].size();
        }
        else if (patchType == processorPolyPatch::typeName)
        {
            // nothing
        }
        else
        {
            const label start = mesh.boundaryMesh()[patchI].start();

            forAll(vf.boundaryField()[patchI], faceI)
            {
                vf.boundaryField()[patchI][faceI] =
                    x[nCells + start + faceI - preEmptyFaces - nIntFaces];
            }
        }
    }

    // Update coupled boundraies
    vf.correctBoundaryConditions();
}


const Map<vector>& solidPolyMesh::pointFixedComponent
(
    const volVectorField& vf
) const
{
    if (!pointFixedComponentPtr_)
    {
        calcPointFixed(vf);
    }

    return *pointFixedComponentPtr_;
}


const Map<symmTensor>& solidPolyMesh::pointFixedDirection
(
    const volVectorField& vf
) const
{
    if (!pointFixedDirectionPtr_)
    {
        calcPointFixed(vf);
    }

    return *pointFixedDirectionPtr_;
}


void solidPolyMesh::clearOutGlobalCoeffs() const
{
    deleteDemandDrivenData(pointProcFacesCoeffsPtr_);
    deleteDemandDrivenData(pointProcBndFacesCoeffsPtr_);
    deleteDemandDrivenData(pointProcCellsCoeffsPtr_);
    deleteDemandDrivenData(gPtNgbProcBndFaceCoeffsPtr_);
    deleteDemandDrivenData(gPtNgbProcCellCoeffsPtr_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool solidPolyMesh::operator!=(const solidPolyMesh& bm) const
{
    return &bm != this;
}


bool solidPolyMesh::operator==(const solidPolyMesh& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
