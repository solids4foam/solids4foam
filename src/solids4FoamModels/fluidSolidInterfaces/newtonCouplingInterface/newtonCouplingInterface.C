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

#include "newtonCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "directMapInterfaceToInterfaceMapping.H"
#include "fixedValueFvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "newtonIcoFluid.H"
#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "meshMotionSolidModelFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newtonCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, newtonCouplingInterface, dictionary
);


// * * * * * * * * * * * Private Member Functions* * * * * * * * * * * * * * //


foamPetscSnesHelper& newtonCouplingInterface::motion()
{
    // meshMotionSolidModelFvMotionSolver is the only currently supported
    // mesh motion solver

    if (!isA<dynamicMotionSolverFvMesh>(fluidMesh()))
    {
        FatalErrorInFunction
            << "meshMotionSolidModelFvMotionSolver is the only currently "
            << "supported mesh motion solver" << exit(FatalError);
    }

    if
    (
        !isA<meshMotionSolidModelFvMotionSolver>
        (
            refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
        )
    )
    {
        FatalErrorInFunction
            << "meshMotionSolidModelFvMotionSolver is the only currently "
            << "supported mesh motion solver" << exit(FatalError);
    }

    // Cast the fluid mesh to a dynamicMotionSolverFvMesh then access its motion
    // solver and cast it to a foamPetscSnesHelper
    // This will only suceed if using a motion solver based on a
    // foamPetscSnesHelper solver
    const foamPetscSnesHelper& motion =
        refCast<const foamPetscSnesHelper>
        (
            refCast<const meshMotionSolidModelFvMotionSolver>
            (
                refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
            ).model()
        );

    // Cast away the const-ness (sorry)
    return const_cast<foamPetscSnesHelper&>(motion);
}


solidModel& newtonCouplingInterface::motionSolid()
{
    return refCast<solidModel>(motion());
}


void newtonCouplingInterface::createSubMatsAndMat
(
    Mat& jac,
    Mat*& subMatsPtr,
    const labelPairList& nBlocksAndBlockSize,
    const labelPairHashSet& nullSubMats
) const
{
    // Set the number of regions
    const label nRegions = nBlocksAndBlockSize.size();

    // Create arrays (vectors) of IS objects for rows and columns.
    // These will indicate where in the matrix the different regions are
    // located
    std::vector<IS> isRow(nRegions), isCol(nRegions);
    label globalOffset = 0;
    for (label r = 0; r < nRegions; ++r)
    {
        const label nBlocks = nBlocksAndBlockSize[r].first();
        const label blockSize = nBlocksAndBlockSize[r].second();
        const label regionSize = nBlocks*blockSize;

        // Create an IS that covers the indices for region r
        ISCreateStride
        (
            PETSC_COMM_WORLD, regionSize, globalOffset, 1, &isRow[r]
        );
        ISCreateStride
        (
            PETSC_COMM_WORLD, regionSize, globalOffset, 1, &isCol[r]
        );
        globalOffset += regionSize;
    }

    // Create an array of submatrices.
    // Each submatrix represents the coupling between region i and
    // region j.
    // For example, subMat[i*nRegions + i] might be the matrix for
    // region i, while subMat[i*nRegions + j] (i != j) are the coupling
    // matrices.
    subMatsPtr = new Mat[nRegions*nRegions];
    for (label i = 0; i < nRegions; ++i)
    {
        for (label j = 0; j < nRegions; ++j)
        {
            // Create the submatrix for regions i and j.
            Mat subA = nullptr;

            // Skip this submatrix as it is specified as null
            labelPair subMatID(i, j);
            if (nullSubMats.found(subMatID))
            {
                continue;
            }

            MatCreate(PETSC_COMM_WORLD, &subA);

            // Number of rows
            const label nRowBlocks = nBlocksAndBlockSize[i].first();
            const label rowBlockSize = nBlocksAndBlockSize[i].second();
            const label nRowsLocal = nRowBlocks*rowBlockSize;
            const label nRowsGlobal = returnReduce(nRowsLocal, sumOp<label>());

            // Number of columns
            const label nColBlocks = nBlocksAndBlockSize[j].first();
            const label colBlockSize = nBlocksAndBlockSize[j].second();
            const label nColsLocal = nColBlocks*colBlockSize;
            const label nColsGlobal = returnReduce(nColsLocal, sumOp<label>());

            if (debug)
            {
                Info<< "subMat(" << i << "," << j << ") " << nl
                    << "nRowsLocal = " << nRowsLocal << nl
                    << "nColsLocal = " << nColsLocal << endl;
            }

            MatSetSizes
            (
                subA,
                nRowsLocal,
                nColsLocal,
                nRowsGlobal,
                nColsGlobal
            );
            MatSetFromOptions(subA);
            MatSetType(subA, MATMPIAIJ);

            // Store subA in row-major order at index [i*nRegions + j]
            subMatsPtr[i*nRegions + j] = subA;
        }
    }

    // Create the nest matrix using the index sets and submatrices
    jac = Mat();
    MatCreateNest
    (
        PETSC_COMM_WORLD,
        nRegions,
        isRow.data(),
        nRegions,
        isCol.data(),
        (Mat*)subMatsPtr,
        &jac
    );

    // Cleanup IS objects.
    for (label r = 0; r < nRegions_; ++r)
    {
        ISDestroy(&isRow[r]);
        ISDestroy(&isCol[r]);
    }
}


label newtonCouplingInterface::initialiseAfm
(
    Mat Afm,
    const fvMesh& fluidMesh,
    const label fluidBlockSize,
    const label motionBlockSize,
    const bool twoD
) const
{
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // The motion appears as the mesh flux in the momentum convection and
    // pressure equation flux divergence terms.
    // This means each fluid equation is coupled to the mesh motion in all
    // cells in the compact stencil

    // CAREFUL: we are setting non-zeros here based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Afm, MATMPIAIJ));

    // Total number of scalar rows in Amf (same as in the fluid region)
    const label scalarRowN = fluidMesh.nCells()*fluidBlockSize;

    // Allocate per-scalar-row nonzeros
    // d is initialised to 1 to count the diagonal
    // while o is initialised to 0
    std::vector<int> d_nnz(scalarRowN, fluidBlockSize);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Count neighbours sharing an internal face
    const Foam::labelUList& own = fluidMesh.owner();
    const Foam::labelUList& nei = fluidMesh.neighbour();
    forAll(own, faceI)
    {
        const Foam::label ownCellID = own[faceI];
        const Foam::label neiCellID = nei[faceI];

        // Increase count for all scalar rows
        label ownRowID = ownCellID*fluidBlockSize;
        d_nnz[ownRowID++] += motionBlockSize;
        d_nnz[ownRowID++] += motionBlockSize;
        d_nnz[ownRowID++] += motionBlockSize;

        label neiRowID = neiCellID*fluidBlockSize;
        d_nnz[neiRowID++] += motionBlockSize;
        d_nnz[neiRowID++] += motionBlockSize;
        d_nnz[neiRowID++] += motionBlockSize;

        if (!twoD)
        {
            d_nnz[ownRowID++] += motionBlockSize;
            d_nnz[neiRowID] += motionBlockSize;
        }
    }

    // Count off-processor neighbour cells
    forAll(fluidMesh.boundary(), patchI)
    {
        if (fluidMesh.boundary()[patchI].type() == "processor")
        {
            const Foam::labelUList& faceCells =
                fluidMesh.boundary()[patchI].faceCells();

            forAll(faceCells, fcI)
            {
                const Foam::label cellID = faceCells[fcI];
                label rowID = cellID*fluidBlockSize;
                o_nnz[rowID] += motionBlockSize;
                o_nnz[rowID] += motionBlockSize;
                o_nnz[rowID] += motionBlockSize;

                if (!twoD)
                {
                    o_nnz[rowID] += motionBlockSize;
                }
            }
        }
        else if (fluidMesh.boundary()[patchI].coupled())
        {
            // Other coupled boundaries are not implemented
            FatalErrorInFunction
                << "Coupled boundary are not implemented, except for"
                << " processor boundaries" << exit(FatalError);
        }
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Afm, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Afm, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    //CHKERRQ(MatSetOption(Afm, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    return 0;
}


label newtonCouplingInterface::initialiseAfs
(
    Mat Afs,
    const fvMesh& fluidMesh,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    // Initialise Afs solid-in-fluid-coupling matrix
    if (debug)
    {
        InfoInFunction
            << "Initialising Afs" << endl;
    }

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Afs (and Asf) matrix is equal to the number of cells at the
    // interface

    // Initially, we only allow one FSI interface
    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Currently, only one interface is allowed when using "
            << typeName << abort(FatalError);
    }

    // Allow only a direct map (conformal interface)
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // CAREFUL: we are setting non-zeros here based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Afs, MATMPIAIJ));

    // Total number of scalar rows in the fluid region
    const label scalarRowN = fluidMesh.nCells()*fluidBlockSize;

    // Allocate per-scalar-row nonzeros, initialised to 0
    std::vector<int> d_nnz(scalarRowN, 0);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface fluid cells
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const labelList& fluidFaceCells =
        fluidMesh.boundary()[fluidPatchID].faceCells();
    forAll(fluidFaceMap, fluidFaceI)
    {
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Calculate the row index for this cell's first scalar equation
        const label rowI = fluidCellID*fluidBlockSize;

        // Set the number of non-zeros to be the number of solid equations
        // e.g., The x-momentum equation could have a coefficient for the
        // solid x/y/z displacement

        // Momentum equation
        d_nnz[rowI] = 1;
        d_nnz[rowI + 1] = 1;
        if (!twoD)
        {
            d_nnz[rowI + 2] = 1;
        }

        // Pressure equation
        d_nnz[rowI + fluidBlockSize - 1] = solidBlockSize;
    }

    // Parallel: we need to decide how to deal with this, e.g. do we allow
    // the same general decompositions as in the partitioned approach
    // For now, only allow serial
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Currently, serial runs are allowed in "
            << typeName << abort(FatalError);
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Afs, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Afs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    //CHKERRQ(MatSetOption(Afs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    return 0;
}


label newtonCouplingInterface::initialiseAms
(
    Mat Ams,
    const fvMesh& fluidMesh,
    const label motionBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Ams matrix is equal to the number of cells at the
    // interface

    // Initially, we only allow one FSI interface
    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Currently, only one interface is allowed when using "
            << typeName << abort(FatalError);
    }

    // Allow only a direct map (conformal interface)
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // CAREFUL: we are setting non-zeros here based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Ams, MATMPIAIJ));

    // Total number of scalar rows in the mesh motion region
    const label scalarRowN = fluidMesh.nCells()*motionBlockSize;

    // Allocate per-scalar-row nonzeros, initialised to 0
    std::vector<int> d_nnz(scalarRowN, 0);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface fluid cells
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const labelList& fluidFaceCells =
        fluidMesh.boundary()[fluidPatchID].faceCells();
    forAll(fluidFaceMap, fluidFaceI)
    {
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Calculate the row index for this cells first scalar equation
        label rowI = fluidCellID*motionBlockSize;

        // Set the number of non-zeros to be the number of solid equations
        d_nnz[rowI++] = solidBlockSize;
        d_nnz[rowI++] = solidBlockSize;

        if (!twoD)
        {
            d_nnz[rowI++] = solidBlockSize;
        }
    }

    // Parallel: we need to decide how to deal with this, e.g. do we allow
    // the same general decompositions as in the partitioned approach
    // For now, only allow serial
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Currently, serial runs are allowed in "
            << typeName << abort(FatalError);
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Ams, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Ams, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    //CHKERRQ(MatSetOption(Ams, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    return 0;
}


label newtonCouplingInterface::initialiseAsf
(
    Mat Asf,
    const fvMesh& solidMesh,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    if (debug)
    {
        InfoInFunction
            << "Initialising" << endl;
    }

    // Allow only a direct map (conformal interface)
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // CAREFUL: we are setting non-zeros her based on the scalar rows, not
    // the block rows

    // Set matrix type to AIJ (since BAIJ does not support non-square
    // blocks)
    CHKERRQ(MatSetType(Asf, MATMPIAIJ));

    // Total number of scalar rows in the solid region
    const label scalarRowN = solidMesh.nCells()*solidBlockSize;

    // Allocate per-scalar-row nonzeros, initialised to 0
    std::vector<int> d_nnz(scalarRowN, 0);
    std::vector<int> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface solid cells
    const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const labelList& solidFaceCells =
        solidMesh.boundary()[solidPatchID].faceCells();
    forAll(solidFaceMap, solidFaceI)
    {
        const label solidCellID = solidFaceCells[solidFaceI];

        // Calculate the row index for this cells first scalar equation
        label rowI = solidCellID*solidBlockSize;

        // Set the number of non-zeros to be the number of fluid equations
        // e.g., The x-momentum equation could have a coefficient for the
        // fluid x/y/z velocity and fluid pressure
        d_nnz[rowI++] = fluidBlockSize;
        d_nnz[rowI++] = fluidBlockSize;

        if (!twoD)
        {
            d_nnz[rowI++] = fluidBlockSize;
        }
    }

    // Parallel: we need to decide how to deal with this, e.g. do we allow
    // the same general decompositions as in the partitioned approach
    // For now, only allow serial
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Currently, serial runs are allowed in "
            << typeName << abort(FatalError);
    }

    // Allocate parallel matrix using AIJ
    CHKERRQ
    (
        MatMPIAIJSetPreallocation(Asf, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Asf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    //CHKERRQ(MatSetOption(Asf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    return 0;
}


label newtonCouplingInterface::formAfm
(
    Mat Afm,
    const label fluidBlockSize,
    const label motionBlockSize,
    const bool twoD
)
{
    // Next, we will manually insert to mesh-in-fluid coupling
    // This represents how the fluid equations depend on the mesh motion
    // The coupling in the momentum equation comes from the convection term
    // where the mesh flux appears. The coupling in the continuity equation
    // appears in the divergence of mesh motion

    // Take references
    const scalarField& phiI = fluid().phi();
    const vectorField& UI = fluid().U();
    const vectorField& SfI = fluidMesh().Sf();
    const scalarField& wI = fluidMesh().weights();
    const labelList& own = fluidMesh().owner();
    const labelList& nei = fluidMesh().neighbour();

    // The fluid and mesh motion refer to the same mesh so they can use the same
    // global cells object
    const globalIndex& globalCells =
        refCast<foamPetscSnesHelper>(fluid()).globalCells();

    // We insert the momentum coupling first
    // For the momentum, we will assume the momentum is 1st order upwind
    // discretised
    // The term is
    //      -div(phi, U) == -div(Sf & (U - meshU)*U)
    // Discretising the meshU term (minuses cancel):
    //     div(Sf & meshU*U)
    //       = sum_faces Sf*meshUf*Uf
    // where meshUf = w*meshUown + (1 - w)*meshUnei
    //
    // Assuming an upwind discretisation,
    // if (Sf & Uf) > 0:
    //     Uf = Uown
    // So, the term for the face becomes
    //     Sf*meshUf*Uown
    //       = Sf*[w*meshUown + (1 - w)*meshUnei]*Uown
    // The coefficient contributions become
    //     for the mesh own cell (differentiating wrt meshUown):
    //         coeff += Sf*w*Uown
    //     for the mesh nei cell (differentiating wrt meshUnei):
    //         coeff += Sf*(1 - w)*Uown
    // We must also add contributions considering face f from the
    // perspective of the fluid nei cell (i.e. Sf in opposite direction)
    //
    // Similarly
    // if (Sf & Uf) <= 0:
    //     Uf = Unei
    // The term for the face becomes
    //     Sf*meshUf*Unei
    //       = Sf*[w*meshUown + (1 - w)*meshUnei]*Unei
    // The coefficient contributions become
    //     for the mesh own cell (differentiating wrt meshUown):
    //         coeff += Sf*w*Unei
    //     for the mesh nei cell (differentiating wrt meshUnei):
    //         coeff += Sf*(1 - w)*Unei
    // We must also add contributions considering face f from the
    // perspective of the fluid nei cell (i.e. Sf in opposite direction)

    tensor coeff = tensor::zero;
    label rowI = -1;
    label colI = -1;
    forAll(own, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI = globalCells.toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI = globalCells.toGlobal(neiCellID);

        // Note: we cannot insert the coefficients as blocks since they are
        // not square (a requirement of MatSetValuesBlocked). Instead we
        // insert scalar coefficients (using MatSetValues)
        if (phiI[faceI] > 0)
        {
            // Add Sf*w*Uown to (fluid own, motion own)
            coeff = SfI[faceI]*wI[faceI]*UI[ownCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add Sf*(1 - w)*Uown to (fluid own, motion nei)
            coeff = SfI[faceI]*(1.0 - wI[faceI])*UI[ownCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*w*Uown to (fluid nei, motion nei)
            coeff = -SfI[faceI]*wI[faceI]*UI[ownCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*(1 - w)*Uown to (fluid nei, motion own)
            coeff = -SfI[faceI]*(1.0 - wI[faceI])*UI[ownCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }
        }
        else
        {
            // Add Sf*w*Unei to (fluid own, motion own)
            coeff = SfI[faceI]*wI[faceI]*UI[neiCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add Sf*(1 - w)*Unei to (fluid own, motion nei)
            coeff = SfI[faceI]*(1.0 - wI[faceI])*UI[neiCellID];
            rowI = globalBlockRowI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*w*Unei to (fluid nei, motion nei)
            coeff = -SfI[faceI]*wI[faceI]*UI[neiCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockColI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }

            // Add -Sf*(1 - w)*Unei to (fluid nei, motion own)
            coeff = -SfI[faceI]*(1.0 - wI[faceI])*UI[neiCellID];
            rowI = globalBlockColI*fluidBlockSize;
            colI = globalBlockRowI*motionBlockSize;
            for (label cmptI = 0; cmptI < motionBlockSize; ++cmptI)
            {
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                    )
                );

                rowI++;
                colI++;
            }
        }
    }

    // Next we insert the continuity coupling
    // -div(phi - meshPhi) + pressure terms
    // Hence, the coefficients are (d/dD)(-div(meshPhi)), where D are the
    // cell centred mesh motion (displacements)
    // div(meshPhi) = sum_faces Sf & (w*Down + (1 - w)*Dnei)
    // So coefficients contributions per face are:
    //     own cell: Sf*w
    //     nei cell: Sf*(1 - w)
    vector coeff2 = vector::zero;
    forAll(own, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI = globalCells.toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI = globalCells.toGlobal(neiCellID);

        // The scalar row index for the pressure equation
        const label rowI =
            globalBlockRowI*fluidBlockSize + (fluidBlockSize - 1);

        // The scalar column index for the first component of the mesh
        // motion equation
        label colI = globalBlockRowI*motionBlockSize;

        // Add w*Sf to owner eqn
        coeff2 = wI[faceI]*SfI[faceI];
        for (int cmptI = 0; cmptI < motionBlockSize; ++cmptI)
        {
            CHKERRQ
            (
                MatSetValues
                (
                    Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                )
            );

            colI++;
        }

        // Scalar column index for neighbour
        colI = globalBlockColI*motionBlockSize;

        // Add (1 - w)*Sf to owner eqn
        coeff2 = (1.0 - wI[faceI])*SfI[faceI];
        for (int cmptI = 0; cmptI < motionBlockSize; ++cmptI)
        {
            CHKERRQ
            (
                MatSetValues
                (
                    Afm, 1, &rowI, 1, &colI, &coeff[cmptI], ADD_VALUES
                )
            );

            colI++;
        }
    }

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "div(meshPhi) coeffs for parallel not yet implemented"
            << exit(FatalError);
    }

    return 0;
}


label newtonCouplingInterface::formAfs
(
    Mat Afs,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
)
{
    if (debug)
    {
        InfoInFunction
            << "Start" << endl;
    }

    // The fluid interface is a prescribed velocity (fixedValue) condition
    // where we assume the fluid wall velocity is equal to the velocity of
    // the adjacent solid cell centre. This approximatation is sufficiently
    // accurate as a preconditioner for the matrix and will not affected the
    // converged solution (which is entirely governed by formResidual)

    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalError
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // Lookup the interface map from the fluid faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each fluid face
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatch = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatch.faceCells();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // Lookup the fluid patch
    const fvPatchVectorField& fluidPatchU =
        fluid().U().boundaryField()[fluidPatchID];
    if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
    {
        FatalErrorInFunction
            << "The fluid interface patch must be of type 'fixedValue'"
            << abort(FatalError);
    }

    if (!isA<fluidModels::newtonIcoFluid>(fluid()))
    {
        FatalErrorInFunction
            << "Currently, the fluid model must be of type 'newtonIcoFluid'"
            << abort(FatalError);
    }

    // Patch viscosity
    const scalarField fluidPatchNuEff
    (
        refCast<fluidModels::newtonIcoFluid>
        (
            fluid()
        ).turbulence().nuEff(fluidPatchID)
    );

    // Fluid interface area vectors
    const vectorField& fluidPatchSf = fluidPatch.Sf();

    // First we will insert the contribution to the fluid momentum equation
    // coming from the diffusion term
    // The known fluid boundary face value is now replaced by the adjacent
    // solid cell velocity
    // Diffusion coefficient is nu*|Sf|/(|n & d|*dt), where dt comes
    // from converting the solid displacement to velocity
    const scalar deltaT = solid().time().deltaTValue();
    const scalarField fluidPatchCoeffs
    (
        fluidPatch.magSf()*fluidPatch.deltaCoeffs()*fluidPatchNuEff/deltaT
    );

    // Second we will insert the contribution to the fluid continuity
    // (pressure) equation, where the div(U) term should use the adjacent
    // solid cell velocity instead of the known boundary face velocity

    forAll(fluidPatch, fluidFaceI)
    {
        const label solidFaceID = fluidFaceMap[fluidFaceI];

        // Fluid and solid cells, which are coupled
        const label fluidCellID = fluidFaceCells[fluidFaceI];
        const label solidCellID = solidFaceCells[solidFaceID];

        // We will add coefficients at block row "fluidCellID" and block
        // column "solidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of fluid matrix
        const label globalBlockRowI =
            refCast<foamPetscSnesHelper>(fluid()).globalCells().toGlobal
            (
                fluidCellID
            );

        // Global block column ID of solid matrix
        const label globalBlockColI =
            refCast<foamPetscSnesHelper>(solid()).globalCells().toGlobal
            (
                solidCellID
            );

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*fluidBlockSize;
        label globalColI = globalBlockColI*solidBlockSize;

        // Momentum coefficient for this face
        PetscScalar value = fluidPatchCoeffs[fluidFaceI];

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        globalColI++;
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            globalRowI++;
            globalColI++;
            CHKERRQ
            (
                MatSetValues
                (
                    Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }

        // Secondly we will insert the contributions for the pressure
        // equation

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        value = -fluidPatchSf[fluidFaceI][vector::X]/deltaT;
        globalRowI++; // pressure equation
        globalColI = globalBlockColI*solidBlockSize; // x solid displacement
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        value = -fluidPatchSf[fluidFaceI][vector::Y]/deltaT;
        globalColI++; // y solid displacement
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            value = -fluidPatchSf[fluidFaceI][vector::Z]/deltaT;
            globalColI++; // y solid displacement
            CHKERRQ
            (
                MatSetValues
                (
                    Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }

    return 0;
}


label newtonCouplingInterface::formAms
(
    Mat Ams,
    const label solidBlockSize,
    const label motionBlockSize,
    const bool twoD
)
{
    // The mesh motion interface is a prescribed displacement (fixedValue)
    // condition where we assume the mesh motion interface displacement is
    // equal to the displacement of the adjacent solid cell centre. This
    // approximatation is sufficiently accurate as a preconditioner for the
    // matrix and will not affected the converged solution (which is entirely
    // governed by formResidual)

    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // Note: the fluid mesh is the same as the motion mesh, but we will use the
    // notation "motion" below rather than fluid to avoid confusion and
    // emphasise that we are referring to the mesh motion equations

    // Lookup the interface map from the motion faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each motion face
    const labelList& motionFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the motion interface patch
    const label motionPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& motionPatch = fluidMesh().boundary()[motionPatchID];
    const labelList& motionFaceCells = motionPatch.faceCells();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // Lookup the motion patch
    const fvPatchVectorField& motionPatchD =
        motionSolid().D().boundaryField()[motionPatchID];
    if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
    {
        FatalErrorInFunction
            << "The motion interface patch must be of type 'fixedValue'"
            << abort(FatalError);
    }

    // Patch impK
    const scalarField motionPatchImpK
    (
        motionSolid().mechanical().impK()().boundaryField()[motionPatchID]
    );

    // For the motion momentum equation, the known interface displacement is now
    // replaced by the adjacent solid cell displacement
    // Diffusion coefficient is impK*|Sf|/(|n & d|)
    const scalarField motionPatchCoeffs
    (
        motionPatch.magSf()*motionPatch.deltaCoeffs()*motionPatchImpK
    );

    forAll(motionPatch, motionFaceI)
    {
        const label solidFaceID = motionFaceMap[motionFaceI];

        // Motion and solid cells, which are coupled
        const label motionCellID = motionFaceCells[motionFaceI];
        const label solidCellID = solidFaceCells[solidFaceID];

        // We will add coefficients at block row "motionCellID" and block
        // column "solidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of motion matrix
        const label globalBlockRowI =
            refCast<foamPetscSnesHelper>(motionSolid()).globalCells().toGlobal
            (
                motionCellID
            );

        // Global block column ID of solid matrix
        const label globalBlockColI =
            refCast<foamPetscSnesHelper>(solid()).globalCells().toGlobal
            (
                solidCellID
            );

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*motionBlockSize;
        label globalColI = globalBlockColI*solidBlockSize;

        // Momentum coefficient for this face
        PetscScalar value = motionPatchCoeffs[motionFaceI];

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        CHKERRQ
        (
            MatSetValues
            (
                Ams, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        globalColI++;
        CHKERRQ
        (
            MatSetValues
            (
                Ams, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            globalRowI++;
            globalColI++;
            CHKERRQ
            (
                MatSetValues
                (
                    Ams, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }

    return 0;
}


label newtonCouplingInterface::formAsf
(
    Mat Asf,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
)
{
    if (debug)
    {
        Info<< "Forming Asf" << endl;
    }

    // The solid interface is a prescribed traction condition, where we
    // approximate the traction on the fluid side of the interface using a
    // compact stencil. For this approximate Jacobian, we assume the traction
    // on a fluid interface face is equal to the pressure at the centre of the
    // adjacent fluid cell. This approximatation is sufficiently
    // accurate as a preconditioner for the matrix and will not affected the
    // converged solution (which is entirely governed by formResidual)

    if (fluidSolidInterface::solidPatchIndices().size() != 1)
    {
        FatalError
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // Lookup the interface map from the solid faces to the fluid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the fluid face ID for each solid face
    const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();

    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatch = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatch.faceCells();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // The approximate force on the solid interface is the solid face area
    // vector multiplied by the adjacent fluid cell centre pressure. We also
    // need to multiply by the fluid density, if the kinematic pressure is
    // used. So the coefficient is the solid face area vector times the
    // fluid density
    // To-do: determine on the fly whether kinematic or dynamic pressure is
    // used
    const fvPatchVectorField& solidPatchD =
        solid().D().boundaryField()[solidPatchID];
    if (!isA<solidTractionFvPatchVectorField>(solidPatchD))
    {
        FatalErrorInFunction
            << "The solid interface patch must be of type 'solidTraction'"
            << abort(FatalError);
    }
    // Todo: add rho() to the fluidModel base class
    const vectorField patchCoeffs
    (
        solidPatch.Sf()
       *refCast<fluidModels::newtonIcoFluid>(fluid()).rho().value()
    );

    forAll(solidPatch, solidFaceI)
    {
        const label fluidFaceID = solidFaceMap[solidFaceI];

        // Fluid and solid cells, which are coupled
        const label solidCellID = solidFaceCells[solidFaceI];
        const label fluidCellID = fluidFaceCells[fluidFaceID];

        // We will add a coefficient at block row "solidCellID" and block
        // column "fluidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of solid matrix
        const label globalBlockRowI =
            refCast<foamPetscSnesHelper>(solid()).globalCells().toGlobal
            (
                solidCellID
            );

        // Global block column ID of fluid matrix
        const label globalBlockColI =
            refCast<foamPetscSnesHelper>(fluid()).globalCells().toGlobal
            (
                fluidCellID
            );

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        // The column corresponds to the pressure equation in the fluid cell
        label globalRowI = globalBlockRowI*solidBlockSize;
        const label globalColI =
            globalBlockColI*fluidBlockSize + (fluidBlockSize - 1);

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        PetscScalar value = -patchCoeffs[solidFaceI][vector::X];
        CHKERRQ
        (
            MatSetValues
            (
                Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        //globalColI++;
        value = -patchCoeffs[solidFaceI][vector::Y];
        CHKERRQ
        (
            MatSetValues
            (
                Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            globalRowI++;
            //globalColI++;
            value = -patchCoeffs[solidFaceI][vector::Z];
            CHKERRQ
            (
                MatSetValues
                (
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }

    return 0;
}


void newtonCouplingInterface::mapInterfaceSolidUToFluidU()
{
    // Lookup the interface map from the fluid faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each fluid face
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    fvPatchVectorField& fluidPatchU =
        fluid().U().boundaryFieldRef()[fluidPatchID];
    if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
    {
        FatalErrorInFunction
            << "The fluid interface patch must be of type 'fixedValue'"
            << abort(FatalError);
    }

    // Lookup the solid interface velocity and displacement
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];

    // Map the solid interface velocity to the fluid interface and the solid
    // displacement to the motion interface
    if (extrapolateSolidInterfaceDisplacement_)
    {
        const fvPatchVectorField& solidPatchU =
            solid().U().boundaryField()[solidPatchID];

        forAll(fluidPatchU, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Extrapolated patch value (larger stencil)
            fluidPatchU[fluidFaceI] = solidPatchU[solidFaceID];
        }
    }
    else
    {
        const labelList& solidFaceCells =
            solidMesh().boundary()[solidPatchID].faceCells();
        const vectorField& solidUI = solid().U();

        forAll(fluidPatchU, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Adjacent cell value
            const label solidCellID = solidFaceCells[solidFaceID];
            fluidPatchU[fluidFaceI] = solidUI[solidCellID];
        }
    }

    fluid().phi() = fvc::interpolate(fluid().U()) & fluidMesh().Sf();
}


void newtonCouplingInterface::mapInterfaceSolidUToMeshMotionD()
{
    // Map solid interface motion to the mesh motion interface

    // Lookup the interface map from the fluid faces to the solid faces
    const interfaceToInterfaceMappings::
        directMapInterfaceToInterfaceMapping& interfaceMap =
        refCast
        <
            const interfaceToInterfaceMappings::
            directMapInterfaceToInterfaceMapping
        >
        (
            interfaceToInterfaceList()[0]
        );

    // The face map gives the solid face ID for each fluid face
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the fluid mesh interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];

    // Lookup the mesh motion displacement field
    fvPatchVectorField& motionPatchD =
        motionSolid().D().boundaryFieldRef()[fluidPatchID];
    if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
    {
        FatalErrorInFunction
            << "The meshMotionFluid interface patch must be of type "
            << "'fixedValue'" << abort(FatalError);
    }

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];

    // Map the solid interface displacement to the motion interface
    if (extrapolateSolidInterfaceDisplacement_)
    {
        const fvPatchVectorField& solidPatchD =
            solid().D().boundaryField()[solidPatchID];

        forAll(motionPatchD, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Extrapolated patch value (larger stencil)
            motionPatchD[fluidFaceI] = solidPatchD[solidFaceID];
        }
    }
    else
    {
        const labelList& solidFaceCells =
            solidMesh().boundary()[solidPatchID].faceCells();
        const vectorField& solidDI = solid().D();

        forAll(motionPatchD, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Adjacent cell value
            const label solidCellID = solidFaceCells[solidFaceID];
            motionPatchD[fluidFaceI] = solidDI[solidCellID];
        }
    }

    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface allowed" << abort(FatalError);
    }

    // Take references to zones
    const standAlonePatch& fluidZone =
        fluid().globalPatches()[0].globalPatch();
    const standAlonePatch& solidZone =
        solid().globalPatches()[0].globalPatch();

    // Get the solid interface patch pointD field
    const vectorField solidPatchPointD
    (
        solid().pointD().boundaryField()[solidPatchID].patchInternalField()
    );

    // Solid point zone pointD
    const vectorField solidZonePointD
    (
        solid().globalPatches()[0].patchPointToGlobal(solidPatchPointD)
    );

    // Prepare the fluid mesh interface pointD field
    vectorField fluidZonePointD(fluidZone.nPoints(), vector::zero);

    // Transfer the point displacements
    interfaceToInterfaceList()[0].transferPointsZoneToZone
    (
        solidZone,
        fluidZone,
        solidZonePointD,
        fluidZonePointD
    );

    // Map the zone field to the patch
    const vectorField meshPatchPointD
    (
        fluid().globalPatches()[0].globalPointToPatch(fluidZonePointD)
    );

    // Set the mesh interface pointD
    // Use "==" to reassign fixedValue
    motionSolid().pointD().boundaryFieldRef()[fluidPatchID] ==
        meshPatchPointD;

    // Correct boundary conditions to enforce the new patch values on the
    // internal field
    motionSolid().pointD().correctBoundaryConditions();
}


void newtonCouplingInterface::predict()
{
    if (nRegions_ != 2)
    {
        notImplemented("Only implemented for nRegions = 2");
    }

    Info<< "Linear predictor" << endl;

    // Predict solution using previous time steps

    // Velocity
    fluid().U() =
        fluid().U().oldTime() + fluid().A()*runTime().deltaT();

    // Pressure
    // Predicting pressure seems to cause instabilities, where the pressure
    // jumps from positive to negative
    // fluid().p() =
    //     fluid().p().oldTime() + fluid().dpdt()*runTime().deltaT();

    // Displacement
    solid().D() =
        solid().D().oldTime() + solid().U()*runTime().deltaT();
    // solid().D() =
    //     solid().D().oldTime() + solid().U()*runTime().deltaT()
    //   + 0.5*sqr(runTime().deltaT())*solid().A();

    // Insert the OpenFOAM fields into the PETSc solution vector

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // The scalar row at which the solid equations start
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

    // Access the raw solution data
    PetscScalar *xx;
    VecGetArray(foamPetscSnesHelper::solution(), &xx);

    // Insert the fluid velocity
    foamPetscSnesHelper::InsertFieldComponents
    (
        fluid().U().primitiveFieldRef(),
        xx,
        0, // Location of U
        fluidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );

    // Insert the fluid pressure
    foamPetscSnesHelper::InsertFieldComponents
    (
        fluid().p().primitiveFieldRef(),
        xx,
        fluidBlockSize - 1, // Location of p component
        fluidBlockSize
    );

    // Insert the displacement
    foamPetscSnesHelper::InsertFieldComponents
    (
        solid().D().primitiveFieldRef(),
        &xx[solidFirstEqnID],
        0, // Location of first component
        solidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );

    // Restore the solution vector
    VecRestoreArray(foamPetscSnesHelper::solution(), &xx);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonCouplingInterface::newtonCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    foamPetscSnesHelper
    (
        fileName
        (
            fsiProperties().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        2*fluid().mesh().nCells() + solid().mesh().nCells(),
        fsiProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true // Will PETSc be used
    ),
    motionSystemScaleFactor_
    (
        fsiProperties().lookupOrDefault<scalar>
        (
            "motionSystemScaleFactor",
            gAverage
            (
                1.0/motionSolid().mechanical().shearModulus()().primitiveField()
            )
        )
    ),
    solidSystemScaleFactor_
    (
        fsiProperties().lookupOrDefault<scalar>
        (
            "solidSystemScaleFactor",
            gAverage(1.0/solid().mechanical().shearModulus()().primitiveField())
        )
    ),
    fluidToSolidCoupling_(fsiProperties().lookup("fluidToSolidCoupling")),
    solidToFluidCoupling_(fsiProperties().lookup("solidToFluidCoupling")),
    meshToFluidCoupling_(fsiProperties().lookup("meshToFluidCoupling")),
    solidToMeshCoupling_(fsiProperties().lookup("solidToMeshCoupling")),
    extrapolateSolidInterfaceDisplacement_
    (
        fsiProperties().lookupOrDefault<Switch>
        (
            "extrapolateSolidInterfaceDisplacement",
            true
        )
    ),
    passViscousStress_(fsiProperties().lookup("passViscousStress")),
    monolithicMeshMotion_(fsiProperties().lookup("monolithicMeshMotion")),
    nRegions_(monolithicMeshMotion_ ? 3 : 2),
    subMatsPtr_(nullptr)
{
    if (monolithicMeshMotion_)
    {
        FatalErrorInFunction
            << "Monolithic mesh motion is still work in progress."
            << " For now, monolithicMeshMotion must be 'off'"
            << exit(FatalError);
    }

    if (solid().twoD() != fluid().twoD())
    {
        FatalErrorInFunction
            << "Either the solid and fluid are both 2-D or both not 2-D"
            << exit(FatalError);
    }

    if
    (
        !isA<foamPetscSnesHelper>(fluid()) || !isA<foamPetscSnesHelper>(solid())
    )
    {
        FatalErrorInFunction
            << "You must use solid and fluid models derived from the "
            << "foamPetscSnesHelper class" << exit(FatalError);
    }

    if (!isA<dynamicMotionSolverFvMesh>(fluidMesh()))
    {
        FatalErrorInFunction
            << "The fluid dynamic mesh must be dynamicMotionSolverFvMesh"
            << exit(FatalError);
    }

    if
    (
        !isA<meshMotionSolidModelFvMotionSolver>
        (
            refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
        )
    )
    {
        FatalErrorInFunction
            << "The fluid mesh motion solver must be "
            << meshMotionSolidModelFvMotionSolver::typeName
            << exit(FatalError);
    }

    if (predictor())
    {
        // Check ddt schemes are not steadyState
        wordList solDdtScheme(3);
        solDdtScheme[0] =
            word(fluidMesh().ddtScheme("ddt(" + fluid().U().name() +')'));
        solDdtScheme[1] =
            word(fluidMesh().ddtScheme("ddt(" + fluid().p().name() +')'));
        solDdtScheme[2] =
            word(solidMesh().ddtScheme("ddt(" + solid().D().name() +')'));

        forAll(solDdtScheme, i)
        {
            if (solDdtScheme[i] == "steadyState")
            {
                FatalErrorIn(type() + "::" + type())
                    << "If predictor is turned on, then the ddtScheme should "
                    << "not be 'steadyState'!" << nl
                    << "Set the predictor to 'off' in fsiProperties or change"
                    << " the ddtSchemes to something other than 'steadyState'"
                    << exit(FatalError);
            }
        }
    }

    Info<< "fluidToSolidCoupling: " << fluidToSolidCoupling_ << nl
        << "solidToFluidCoupling: " << solidToFluidCoupling_ << nl
        << "meshToFluidCoupling: " << meshToFluidCoupling_ << nl
        << "solidToMeshCoupling: " << solidToMeshCoupling_ << nl
        << "extrapolateSolidInterfaceDisplacement: "
        << extrapolateSolidInterfaceDisplacement_ << nl
        << "passViscousStress: " << passViscousStress_ << nl
        << "monolithicMeshMotion: " << monolithicMeshMotion_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void newtonCouplingInterface::setDeltaT(Time& runTime)
{
    if
    (
        runTime.controlDict().getOrDefault("adjustTimeStep", false)
     && runTime.timeIndex() > 0
    )
    {
        // Adjust the time step based on the number of outer (Newton) iterations
        // required by the PETSc SNES solver
        // We will aim to keep the number of iterations within the range of
        // minTargetNIter to maxTargetNIter

        const scalar maxDeltaT =
            runTime.controlDict().get<scalar>("maxDeltaT");
        const scalar minDeltaT =
            runTime.controlDict().get<scalar>("minDeltaT");

        const int minTargetNIter =
            runTime.controlDict().getOrDefault<int>("minTargetNIter", 3);
        const int maxTargetNIter =
            runTime.controlDict().getOrDefault<int>("maxTargetNIter", 6);

        PetscInt numIter;
        SNESGetIterationNumber(foamPetscSnesHelper::snes(), &numIter);

        scalar newDeltaT = runTime.deltaTValue();

        if (numIter > maxTargetNIter)
        {
            newDeltaT = max(0.5*newDeltaT, minDeltaT);
        }
        else if (numIter < minTargetNIter)
        {
            newDeltaT = min(1.5*newDeltaT, maxDeltaT);
        }

        Info<< "Old time step = " << runTime.deltaTValue() << nl
            << "New time step = " << newDeltaT << nl << endl;

        runTime.setDeltaT(newDeltaT);
    }
}


bool newtonCouplingInterface::evolve()
{
    // Steps
    // 1. Optional: predict solution field using old-time fields
    // 2. Map foam solution fields to PETSc
    // 3. Solve PETSc system
    // 4. Map PETSc solution to foam fields
    // 5. Update an secondary foam fields

    // Check if coupling switch needs to be updated
    if (!coupled())
    {
        updateCoupled();
    }

    // Solution predictor
    if (predictor() && newTimeStep())
    {
        predict();
    }

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    if (nRegions_ == 2)
    {
        Info<< "Solving the fluid-solid equations for U, p and D using PETSc "
            << "SNES" << endl;

        // The scalar row at which the solid equations start
        const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

        // fluid-solid + mesh motion loop
        PetscReal initResidualNorm = 0;
        scalar residualNorm = 0;
        scalar solutionNorm = 0;
        Vec solutionDelta;
        scalar solutionDeltaNorm = 0;
        label iMeshCorr = 1;
        const label nMeshCorr(readInt(fsiProperties().lookup("nMeshCorr")));
        const scalar meshCorrRelTol
        (
            readScalar(fsiProperties().lookup("meshCorrRelTol"))
        );
        const scalar meshCorrSolTol
        (
            readScalar(fsiProperties().lookup("meshCorrSolTol"))
        );
        const scalar meshCorrAbsTol
        (
            readScalar(fsiProperties().lookup("meshCorrAbsTol"))
        );
        do
        {
            Info<< "Time = " << runTime().value()
                << ", Mesh loop corr = " << iMeshCorr << endl;

            // Solve the monolithic fluid-solid system
            foamPetscSnesHelper::solve();

            // Access the raw solution data
            const PetscScalar *xx;
            VecGetArrayRead(foamPetscSnesHelper::solution(), &xx);

            // Retrieve the fluid velocity
            foamPetscSnesHelper::ExtractFieldComponents
            (
                xx,
                fluid().U().primitiveFieldRef(),
                0, // Location of U
                fluidBlockSize,
                twoD ? labelList({0,1}) : labelList({0,1,2})
            );
            fluid().U().correctBoundaryConditions();

            // Retreive the fluid pressure
            foamPetscSnesHelper::ExtractFieldComponents
            (
                xx,
                fluid().p().primitiveFieldRef(),
                fluidBlockSize - 1, // Location of p component
                fluidBlockSize
            );
            fluid().p().correctBoundaryConditions();

            // Extract the displacement, which starts at row fluid.nCells*fluidBlockSize
            foamPetscSnesHelper::ExtractFieldComponents
            (
                &xx[solidFirstEqnID],
                solid().D().primitiveFieldRef(),
                0, // Location of first component
                solidBlockSize,
                twoD ? labelList({0,1}) : labelList({0,1,2})
            );
            solid().D().correctBoundaryConditions();

            // Restore the solution vector
            VecRestoreArrayRead(foamPetscSnesHelper::solution(), &xx);

            // Update gradient of displacement
            solid().mechanical().grad(solid().D(), solid().gradD());

            // Interpolate cell displacements to vertices
            solid().mechanical().interpolate(solid().D(), solid().gradD(), solid().pointD());

            // Increment of displacement
            solid().DD() = solid().D() - solid().D().oldTime();

            // Increment of point displacement
            solid().pointDD() = solid().pointD() - solid().pointD().oldTime();

            // Velocity
            solid().U() = fvc::ddt(solid().D());

            // Update the fluid interface velocity
            if (solidToFluidCoupling_)
            {
                mapInterfaceSolidUToFluidU();
            }

            // Update the mesh motion interface
            if (solidToMeshCoupling_)
            {
                mapInterfaceSolidUToMeshMotionD();
            }

            // Solve the mesh motion equations
            fluidMesh().update();

            // Calculate the solution and correction (delta) norms
            VecNorm(foamPetscSnesHelper::solution(), NORM_2, &solutionNorm);
            SNESGetSolutionUpdate(foamPetscSnesHelper::snes(), &solutionDelta);
            VecNorm(solutionDelta, NORM_2, &solutionDeltaNorm);

            // Extract the residual
            SNESGetFunctionNorm(foamPetscSnesHelper::snes(), &residualNorm);

            // Record the initial residual and solution norms
            if (iMeshCorr == 1)
            {
                initResidualNorm = residualNorm;
            }

            Info<< "Mesh corrector loop residual norm = " << residualNorm
                << ", solution delta norm = " << solutionDeltaNorm
                << ", solution norm = " << solutionNorm << nl
                << endl;

            if (iMeshCorr == (nMeshCorr - 1))
            {
                FatalErrorInFunction
                    << "Max mesh correctors reached!" << exit(FatalError);
            }
        }
        while
        (
            iMeshCorr++ < nMeshCorr
         && residualNorm > meshCorrRelTol*initResidualNorm
         && solutionDeltaNorm > meshCorrSolTol*solutionNorm
         && residualNorm > meshCorrAbsTol
        );
    }
    else if (nRegions_ == 3)
    {
        // Mesh motion block size
        const label motionBlockSize = solidBlockSize;

        // The scalar row at which the motion equations start
        const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

        // The scalar row at which the solid equations start
        const label solidFirstEqnID =
            motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;

        // Solve the nonlinear system and check the convergence
        foamPetscSnesHelper::solve();

        // Access the raw solution data
        const PetscScalar *xx;
        VecGetArrayRead(foamPetscSnesHelper::solution(), &xx);

        // Retrieve the fluid velocity
        foamPetscSnesHelper::ExtractFieldComponents
        (
            xx,
            fluid().U().primitiveFieldRef(),
            0, // Location of U
            fluidBlockSize,
            twoD ? labelList({0,1}) : labelList({0,1,2})
        );
        fluid().U().correctBoundaryConditions();

        // Retrieve the fluid pressure
        foamPetscSnesHelper::ExtractFieldComponents
        (
            xx,
            fluid().p().primitiveFieldRef(),
            fluidBlockSize - 1, // Location of p component
            fluidBlockSize
        );
        fluid().p().correctBoundaryConditions();

        // Retrieve the fluid mesh motion
        foamPetscSnesHelper::ExtractFieldComponents
        (
            &xx[motionFirstEqnID],
            motionSolid().D().primitiveFieldRef(),
            0, // Location of first component
            motionBlockSize,
            twoD ? labelList({0,1}) : labelList({0,1,2})
        );
        motionSolid().D().correctBoundaryConditions();

        // Retrieve the solid displacement
        foamPetscSnesHelper::ExtractFieldComponents
        (
            &xx[solidFirstEqnID],
            solid().D().primitiveFieldRef(),
            0, // Location of first component
            solidBlockSize,
            twoD ? labelList({0,1}) : labelList({0,1,2})
        );
        solid().D().correctBoundaryConditions();

        // Restore the x vector
        VecRestoreArrayRead(foamPetscSnesHelper::solution(), &xx);

        // Update gradient of displacement
        solid().mechanical().grad(solid().D(), solid().gradD());
        motionSolid().mechanical().grad(motionSolid().D(), motionSolid().gradD());

        // Interpolate cell displacements to vertices
        solid().mechanical().interpolate(solid().D(), solid().gradD(), solid().pointD());
        motionSolid().mechanical().interpolate
        (
            motionSolid().D(), motionSolid().gradD(), motionSolid().pointD()
        );

        // Increment of displacement
        solid().DD() = solid().D() - solid().D().oldTime();

        // Increment of point displacement
        solid().pointDD() = solid().pointD() - solid().pointD().oldTime();

        // Velocity
        solid().U() = fvc::ddt(solid().D());

        // Update the fluid interface velocity
        if (solidToFluidCoupling_)
        {
            mapInterfaceSolidUToFluidU();
        }

        // Update the mesh motion interface
        if (solidToMeshCoupling_)
        {
            mapInterfaceSolidUToMeshMotionD();
        }

        // Move the fluid mesh
        if (meshToFluidCoupling_)
        {
            const vectorField& points0 =
                refCast<const meshMotionSolidModelFvMotionSolver>
                (
                    refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
                ).points0();
            fluidMesh().movePoints(points0 + motionSolid().pointD());
            fluid().U().correctBoundaryConditions();
        }
    }
    else
    {
        FatalErrorInFunction
            << "nRegions must be 2 or 3!" << abort(FatalError);
    }

    // Update fluid velocity acceleration
    fluid().A() = fvc::ddt(fluid().U());

    // Update fluid rate of change of pressure
    fluid().dpdt() = fvc::ddt(fluid().p());

    return 0;
}


label newtonCouplingInterface::initialiseJacobian(Mat& jac)
{
    // A fluid-solid interaction problem with a moving mesh (arbitrary
    // Lagrangian Eulerian) fluid can be expressed in Ax=b form, where the
    // matrix A is
    //
    //     *-----------------*
    //     | Aff | Afm | Afs |
    //     |-----------------|
    //     |  0  | Amm | Ams |
    //     |-----------------|
    //     | Asf |  0  | Ass |
    //     *-----------------*
    //
    // The diagonal submatrices are
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations
    //
    // The off diagonal submatrices, which represent coupling between
    // equations, are
    //  - Afm: mesh motion terms in the fluid equations
    //  - Afs: solid terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations
    //
    // Note
    //  - Amf is empty as the the mesh motion does not directly
    //    depend on the fluid solution
    //  - Asm is empty as the solid equations do not depend on the fluid
    //    mesh motion

    // The full matrix can be formed, although the mesh-to-fluid coupling (Afm)
    // is not trivial to derive. As the coupling between the solid and fluid
    // (given by Afs and Asf) is typically stronger than the mesh-to-fluid
    // coupling (Afm), we can also solve a smaller fluid-solid monolithic system
    // followed sequentially by the fluid mesh motion.
    // We will consider both approaches here, starting with the sequential
    // fluid-solid + mesh motion first (monolithicMeshMotion = "off") and later
    // add the fluid-solid-mesh option (monolithicMeshMotion = "on"). The
    // fluid-solid monolithic matrix takes the form:
    //
    //     *-----------*
    //     | Aff | Afs |
    //     |-----------|
    //     | Asf | Ass |
    //     *-----------*


    // In both cases, we will create a nested matrix, where the overall
    // monolithic matrix is formed by smaller sub-matrices, e.g., Aff, Ass, etc.

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    if (nRegions_ == 2)
    {
        // Our monolithic system matrix will take the form:
        //
        //     *-----------*
        //     | Aff | Afs |
        //     |-----------|
        //     | Asf | Ass |
        //     *-----------*

        // For brevity and convenience, we will store the size and blockSize of
        // the regions in a labelPairList
        const labelPairList nBlocksAndBlockSize
        (
            {
                {fluidMesh().nCells(), fluidBlockSize},
                {solidMesh().nCells(), solidBlockSize}
            }
        );

        // Create the empty submatrices
        createSubMatsAndMat(jac, subMatsPtr_, nBlocksAndBlockSize);

        // Get access to the sub-matrices
        PetscInt nr, nc;
        Mat **subMats;
        CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

        // Initialise the diagonal submatrices:
        //  - Aff: fluid equations (momentum and pressure)
        //  - Ass: solid equations

        // Aff
        Foam::initialiseJacobian
        (
            subMats[0][0], fluid().mesh(), fluidBlockSize, false
        );
        PetscObjectSetName((PetscObject)subMats[0][0], "Aff");

        // Ass
        Foam::initialiseJacobian
        (
            subMats[1][1], solid().mesh(), solidBlockSize, false
        );
        PetscObjectSetName((PetscObject)subMats[1][1], "Ass");

        // Initialise the off diagonal submatrices:
        //  - Afs: solid terms in the fluid equations
        //  - Asf: fluid terms (pressure on interface) in solid equations

        // Afs
        initialiseAfs
        (
            subMats[0][1], fluidMesh(), fluidBlockSize, solidBlockSize, twoD
        );
        PetscObjectSetName((PetscObject)subMats[0][1], "Afs");

        // Asf
        initialiseAsf
        (
            subMats[1][0], solidMesh(), fluidBlockSize, solidBlockSize, twoD
        );
        PetscObjectSetName((PetscObject)subMats[1][0], "Asf");
    }
    else if (nRegions_ == 3)
    {
        // Our monolithic system matrix will take the form:
        //
        //     *-----------------*
        //     | Aff | Afm | Afs |
        //     |-----------------|
        //     |  0  | Amm | Ams |
        //     |-----------------|
        //     | Asf |  0  | Ass |
        //     *-----------------*
        //
        // TODO: move the mesh motion equation to the end so the upper left
        // fluid-solid system is the same as the nRegions=2 case

        // Mesh motion block size
        const label motionBlockSize = solidBlockSize;

        // For brevity and convenience, we will store the size and blockSize of the
        // regions in a labelPairList
        const labelPairList nBlocksAndBlockSize
        (
            {
                {fluidMesh().nCells(), fluidBlockSize},
                {fluidMesh().nCells(), motionBlockSize},
                {solidMesh().nCells(), solidBlockSize}
            }
        );

        // Create the empty submatrices and skip Amf and Asm, which will be set
        // to null ptrs
        labelPairHashSet nullSubMats;
        nullSubMats.insert(labelPair(1, 0)); // Amf
        nullSubMats.insert(labelPair(2, 1)); // Asm
        createSubMatsAndMat(jac, subMatsPtr_, nBlocksAndBlockSize, nullSubMats);

        // Get access to the sub-matrices
        PetscInt nr, nc;
        Mat **subMats;
        CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

        // Initialise the diagonal submatrices:
        //  - Aff: fluid equations (momentum and pressure)
        //  - Amm: mesh motion equations in the fluid domain
        //  - Ass: solid equations

        // Aff
        Foam::initialiseJacobian
        (
            subMats[0][0], fluid().mesh(), fluidBlockSize, false
        );
        PetscObjectSetName((PetscObject)subMats[0][0], "Aff");

        // Amm
        Foam::initialiseJacobian
        (
            subMats[1][1], fluid().mesh(), motionBlockSize, false
        );
        PetscObjectSetName((PetscObject)subMats[1][1], "Amm");

        // Ass
        Foam::initialiseJacobian
        (
            subMats[2][2], solid().mesh(), solidBlockSize, false
        );
        PetscObjectSetName((PetscObject)subMats[2][2], "Ass");

        // Initialise the four off diagonal submatrices:
        //  - Afm: mesh motion terms in the fluid equations
        //  - Afs: solid terms in the fluid equations
        //  - Ams: solid terms (interface motion) in the mesh motion equation
        //  - Asf: fluid terms (pressure on interface) in solid equations

        // Afm
        initialiseAfm
        (
            subMats[0][1], fluidMesh(), fluidBlockSize, motionBlockSize, twoD
        );
        PetscObjectSetName((PetscObject)subMats[0][1], "Afm");

        // Afs
        initialiseAfs
        (
            subMats[0][2], fluidMesh(), fluidBlockSize, solidBlockSize, twoD
        );
        PetscObjectSetName((PetscObject)subMats[0][2], "Afs");

        // Ams
        initialiseAms
        (
            subMats[1][2], fluidMesh(), motionBlockSize, solidBlockSize, twoD
        );
        PetscObjectSetName((PetscObject)subMats[1][2], "Ams");

        // Asf
        initialiseAsf
        (
            subMats[2][0], solidMesh(), fluidBlockSize, solidBlockSize, twoD
        );
        PetscObjectSetName((PetscObject)subMats[2][0], "Asf");
    }
    else
    {
        FatalErrorInFunction
            << "nRegions must be 2 or 3!" << abort(FatalError);
    }

    return 0;
}


label newtonCouplingInterface::initialiseSolution(Vec& x)
{
    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // Count the number of unknowns
    label n = 0;
    label N = 0;
    if (nRegions_ == 2)
    {
        // For brevity and convenience, we will store the size and blockSize of the
        // regions in a labelPairList
        const labelPairList nBlocksAndBlockSize
        (
            {
                {fluidMesh().nCells(), fluidBlockSize},
                {solidMesh().nCells(), solidBlockSize}
            }
        );

        // Count number of local blocks and local scalar equations
        forAll(nBlocksAndBlockSize, regionI)
        {
            const label nBlocksRegionI = nBlocksAndBlockSize[regionI].first();
            const label blockSizeRegionI = nBlocksAndBlockSize[regionI].second();
            n += nBlocksRegionI*blockSizeRegionI;
        }

        // Global system size: total number of scalar equation across all
        // processors
        N = returnReduce(n, sumOp<label>());
    }
    else if (nRegions_ == 3)
    {
        // Mesh motion block size
        const label motionBlockSize = solidBlockSize;

        // For brevity and convenience, we will store the size and blockSize of the
        // regions in a labelPairList
        const labelPairList nBlocksAndBlockSize
        (
            {
                {fluidMesh().nCells(), fluidBlockSize},
                {fluidMesh().nCells(), motionBlockSize},
                {solidMesh().nCells(), solidBlockSize}
            }
        );

        // Count number of local blocks and local scalar equations
        forAll(nBlocksAndBlockSize, regionI)
        {
            const label nBlocksRegionI = nBlocksAndBlockSize[regionI].first();
            const label blockSizeRegionI = nBlocksAndBlockSize[regionI].second();
            n += nBlocksRegionI*blockSizeRegionI;
        }

        // Global system size: total number of scalar equation across all
        // processors
        N = returnReduce(n, sumOp<label>());
    }
    else
    {
        FatalErrorInFunction
            << "nRegions must be 2 or 3!" << abort(FatalError);
    }

    // Create solution vector
    x = Vec();
    CHKERRQ(VecCreate(PETSC_COMM_WORLD, &x));
    CHKERRQ(VecSetSizes(x, n, N));
    CHKERRQ(VecSetType(x, VECMPI));
    CHKERRQ(PetscObjectSetName((PetscObject) x, "Solution"));
    CHKERRQ(VecSetFromOptions(x));
    CHKERRQ(VecZeroEntries(x));

    return 0;
}


label newtonCouplingInterface::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
)
{
    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    if (nRegions_ == 2)
    {
        // Considerations on the order of updating the fluid and solid residuals
        //  - solid depends on the fluid interface traction
        //  - fluid depends on the solid interface velocity
        //
        // Approach
        // 1. Update the fluid velocity and pressure fields and calculate the
        //    traction at the fluid interface
        // 2. Map the fluid interface traction to the solid interface
        // 3. Update the solid residual, which now has the correct interface
        //    traction
        // 4. Map the solid interface velocity to the fluid interface
        // 5. Update the fluid residual, which now has the correct interface
        //    velocity
        // 6. Apply scaling factor to solid equations to preserve the
        //    condition number of the monolithic system

        // The scalar row at which the solid equations start
        const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

        // Currently limited to one interface: it should be straight-forward to add
        // a loop over multiple interface => todo
        if (fluidSolidInterface::fluidPatchIndices().size() != 1)
        {
            FatalErrorInFunction
                << "Only one interface patch is currently allowed"
                << abort(FatalError);
        }

        // 1. Update the fluid velocity and pressure fields and calculate the
        //    traction at the fluid interface
        {
            // Take references
            volVectorField& U = fluid().U();
            volScalarField& p = fluid().p();

            // Retrieve the solution
            // Map the PETSc solution to the U field
            foamPetscSnesHelper::ExtractFieldComponents<vector>
            (
                x,
                U.primitiveFieldRef(),
                0, // Location of U
                fluidBlockSize,
                fluid().twoD() ? labelList({0,1}) : labelList({0,1,2})
            );

            U.correctBoundaryConditions();

            const_cast<volTensorField&>
            (
                fluidMesh().lookupObject<volTensorField>("grad(U)")
            ) = fvc::grad(U);

            U.correctBoundaryConditions();

            // Map the PETSc solution to the p field
            // p is located in the final component
            foamPetscSnesHelper::ExtractFieldComponents<scalar>
            (
                x,
                p.primitiveFieldRef(),
                fluidBlockSize - 1, // Location of p component
                fluidBlockSize
            );

            p.correctBoundaryConditions();
        }

        // 2. Map the fluid interface traction to the solid interface
        if (fluidToSolidCoupling_ && coupled())
        {
            // Fluid interface traction
            const label fluidPatchID =
                fluidSolidInterface::fluidPatchIndices()[0];
            const vectorField fluidNf
            (
                fluidMesh().boundary()[fluidPatchID].nf()
            );
            vectorField fluidTraction
            (
              - fluid().patchPressureForce(fluidPatchID)*fluidNf
            );
            if (passViscousStress_)
            {
                fluidTraction += fluid().patchViscousForce(fluidPatchID);
            }

            // Lookup the interface map from the fluid faces to the solid faces
            const interfaceToInterfaceMappings::
                directMapInterfaceToInterfaceMapping& interfaceMap =
                refCast
                <
                    const interfaceToInterfaceMappings::
                    directMapInterfaceToInterfaceMapping
                >
                (
                    interfaceToInterfaceList()[0]
                );

            // The face map gives the solid face ID for each fluid face
            const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

            // Calculate the solid interface traction
            // Flip the sign as the solid normals point in the opposite direction to
            // the fluid normals
            const label solidPatchID =
                fluidSolidInterface::solidPatchIndices()[0];
            vectorField solidTraction
            (
                solidMesh().boundary()[solidPatchID].size()
            );
            forAll(fluidTraction, fluidFaceI)
            {
                solidTraction[fluidFaceMap[fluidFaceI]] =
                    -fluidTraction[fluidFaceI];
            }

            // Lookup the displacement interface traction patch and set the traction
            fvPatchVectorField& solidPatchD =
                solid().D().boundaryFieldRef()[solidPatchID];
            if (!isA<solidTractionFvPatchVectorField>(solidPatchD))
            {
                FatalErrorInFunction
                    << "The solidinterface patch must be of type "
                    << "'solidTraction'"
                    << abort(FatalError);
            }

            solidTractionFvPatchVectorField& solidTractionPatch =
                refCast<solidTractionFvPatchVectorField>(solidPatchD);

            solidTractionPatch.traction() = solidTraction;
        }

        // 3. Update the solid residual, which now has the correct interface
        //    traction
        refCast<foamPetscSnesHelper>(solid()).formResidual
        (
            &f[solidFirstEqnID], &x[solidFirstEqnID]
        );

        // 4. Map the solid interface velocity to the fluid interface
        // Note: it is assumed the solid.U() is dD/dt even for steady-state cases
        if (solidToFluidCoupling_ && coupled())
        {
            mapInterfaceSolidUToFluidU();
        }

        // 5. Update the fluid residual, which now has the correct interface
        // velocity
        // Note that the fluid equations are first in the f (residual) and x
        // (solution) lists
        refCast<foamPetscSnesHelper>(fluid()).formResidual(f, x);

        // 6. Apply scaling factor to solid equations to preserve the
        // condition number of the monolithic system
        const label solidSystemEnd =
            solidFirstEqnID + solidMesh().nCells()*solidBlockSize;
        for (int i = solidFirstEqnID; i < solidSystemEnd; ++i)
        {
            f[i] *= solidSystemScaleFactor_;
        }
    }
    else if (nRegions_ == 3)
    {
        // Considerations on the order of updating the fluid, solid, and mesh motion
        // residuals
        //  - solid depends on the fluid interface traction
        //  - mesh motion depends on the solid interface displacement
        //  - fluid depends on the solid interface velocity
        //  - fluid depends on the entire mesh motion flux field
        //
        // Approach
        // 1. Update the fluid velocity and pressure fields and calculate the
        //    traction at the fluid interface
        // 2. Map the fluid interface traction to the solid interface
        // 3. Update the solid residual, which now has the correct interface
        //    traction
        // 4. Map the solid interface velocity to the fluid interface
        // 5. Map the solid interface displacement to the mesh motion interface
        // 6. Update the mesh motion residual, which now has the correct interface
        //    displacement
        // 7. Move the fluid mesh using the mesh motion field
        // 8. Update the fluid residual, which now has the correct interface
        //    velocity and mesh motion
        // 9. Apply scaling factor to solid and motion equations to preserve the
        //    condition number of the monolithic system

        const label motionBlockSize = solidBlockSize;

        // The scalar row at which the motion equations start
        const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

        // The scalar row at which the solid equations start
        const label solidFirstEqnID =
            motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;

        // Currently limited to one interface: it should be straight-forward to add
        // a loop over multiple interface => todo
        if (fluidSolidInterface::fluidPatchIndices().size() != 1)
        {
            FatalErrorInFunction
                << "Only one interface patch is currently allowed"
                << abort(FatalError);
        }

        // 1. Update the fluid velocity and pressure fields and calculate the
        //    traction at the fluid interface
        {
            // Take references
            volVectorField& U = fluid().U();
            volScalarField& p = fluid().p();

            // Retrieve the solution
            // Map the PETSc solution to the U field
            foamPetscSnesHelper::ExtractFieldComponents<vector>
            (
                x,
                U.primitiveFieldRef(),
                0, // Location of U
                fluidBlockSize,
                fluid().twoD() ? labelList({0,1}) : labelList({0,1,2})
            );

            U.correctBoundaryConditions();

            const_cast<volTensorField&>
            (
                fluidMesh().lookupObject<volTensorField>("grad(U)")
            ) = fvc::grad(U);

            U.correctBoundaryConditions();

            // Map the PETSc solution to the p field
            // p is located in the final component
            foamPetscSnesHelper::ExtractFieldComponents<scalar>
            (
                x,
                p.primitiveFieldRef(),
                fluidBlockSize - 1, // Location of p component
                fluidBlockSize
            );

            p.correctBoundaryConditions();
        }

        // 2. Map the fluid interface traction to the solid interface
        if (fluidToSolidCoupling_ && coupled())
        {
            // Fluid interface traction
            const label fluidPatchID =
                fluidSolidInterface::fluidPatchIndices()[0];
            const vectorField fluidNf
            (
                fluidMesh().boundary()[fluidPatchID].nf()
            );
            vectorField fluidTraction
            (
              - fluid().patchPressureForce(fluidPatchID)*fluidNf
            );
            if (passViscousStress_)
            {
                fluidTraction += fluid().patchViscousForce(fluidPatchID);
            }

            // Lookup the interface map from the fluid faces to the solid faces
            const interfaceToInterfaceMappings::
                directMapInterfaceToInterfaceMapping& interfaceMap =
                refCast
                <
                    const interfaceToInterfaceMappings::
                    directMapInterfaceToInterfaceMapping
                >
                (
                    interfaceToInterfaceList()[0]
                );

            // The face map gives the solid face ID for each fluid face
            const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

            // Calculate the solid interface traction
            // Flip the sign as the solid normals point in the opposite direction to
            // the fluid normals
            const label solidPatchID =
                fluidSolidInterface::solidPatchIndices()[0];
            vectorField solidTraction
            (
                solidMesh().boundary()[solidPatchID].size()
            );
            forAll(fluidTraction, fluidFaceI)
            {
                solidTraction[fluidFaceMap[fluidFaceI]] =
                    -fluidTraction[fluidFaceI];
            }

            // Lookup the displacement interface traction patch and set the traction
            fvPatchVectorField& solidPatchD =
                solid().D().boundaryFieldRef()[solidPatchID];
            if (!isA<solidTractionFvPatchVectorField>(solidPatchD))
            {
                FatalErrorInFunction
                    << "The solidinterface patch must be of type "
                    << "'solidTraction'"
                    << abort(FatalError);
            }

            solidTractionFvPatchVectorField& solidTractionPatch =
                refCast<solidTractionFvPatchVectorField>(solidPatchD);

            solidTractionPatch.traction() = solidTraction;
        }

        // 3. Update the solid residual, which now has the correct interface
        //    traction
        refCast<foamPetscSnesHelper>(solid()).formResidual
        (
            &f[solidFirstEqnID], &x[solidFirstEqnID]
        );

        // 4. Map the solid interface velocity to the fluid interface
        // Note: it is assumed the solid.U() is dD/dt even for steady-state cases
        if (solidToFluidCoupling_ && coupled())
        {
            mapInterfaceSolidUToFluidU();
        }

        // 5. Map the solid interface displacement to the mesh motion interface
        if (solidToMeshCoupling_ && coupled())
        {
            mapInterfaceSolidUToMeshMotionD();
        }

        // 6. Update the mesh motion residual, which now has the correct interface
        //    displacement
        refCast<foamPetscSnesHelper>(motionSolid()).formResidual
        (
            &f[motionFirstEqnID], &x[motionFirstEqnID]
        );

        // 7. Move the fluid mesh using the mesh motion field
        autoPtr<vectorField> oldPointsPtr_;
        if (meshToFluidCoupling_ && coupled())
        {
            // Store old points
            oldPointsPtr_.set(new vectorField(fluidMesh().points()));
            const vectorField& points0 =
                refCast<const meshMotionSolidModelFvMotionSolver>
                (
                    refCast<dynamicMotionSolverFvMesh>(fluidMesh()).motion()
                ).points0();

            // Move the mesh
            fluidMesh().movePoints(points0 + motionSolid().pointD());
            fluid().U().correctBoundaryConditions();
        }

        // 8. Update the fluid residual, which now has the correct interface
        //    velocity and mesh motion
        // Note that the fluid equations are first in the f (residual) and x
        // (solution) lists
        refCast<foamPetscSnesHelper>(fluid()).formResidual(f, x);

        // Reset the mesh
        if (oldPointsPtr_)
        {
            fluidMesh().movePoints(oldPointsPtr_.ref());
        }

        // 9. Apply scaling factor to solid and motion equations to preserve the
        // condition number of the monolithic system

        const label solidSystemEnd =
            solidFirstEqnID + solidMesh().nCells()*solidBlockSize;
        for (int i = solidFirstEqnID; i < solidSystemEnd; ++i)
        {
            f[i] *= solidSystemScaleFactor_;
        }

        const label motionSystemEnd =
            motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;
        for (int i = motionFirstEqnID; i < motionSystemEnd; ++i)
        {
            f[i] *= motionSystemScaleFactor_;
        }
    }
    else
    {
        FatalErrorInFunction
            << "nRegions must be 2 or 3!" << abort(FatalError);
    }

    return 0;
}


label newtonCouplingInterface::formJacobian
(
    Mat jac,
    const PetscScalar *x
)
{
    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    if (nRegions_ == 2)
    {
        // We will assembly the four sub-matrices:

        // Diagonal submatrices:
        //  - Aff: fluid equations (momentum and pressure)
        //  - Ass: solid equations
        //
        // Off-diagonal submatrices:
        //  - Afs: solid terms in the fluid equations
        //  - Asf: fluid terms (pressure on interface) in solid equations

        // Zero entries
        MatZeroEntries(jac);

        // Get access to the sub-matrices
        PetscInt nr, nc;
        Mat **subMats;
        CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

        if (nr != 2 || nc != 2)
        {
            FatalErrorInFunction
                << "The matrix has the wrong number of sub matrices: "
                << "nr = " << nr << ", nc = " << nc << abort(FatalError);
        }

        // The scalar row at which the solid equations start
        const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

        // Form diagonal submatrices
        //  - Aff: fluid equations (momentum and pressure)
        //  - Amm: mesh motion equations in the fluid domain
        //  - Ass: solid equations
        //

        // Aff
        refCast<foamPetscSnesHelper>(fluid()).formJacobian
        (
            subMats[0][0], x
        );

        // Ass
        refCast<foamPetscSnesHelper>(solid()).formJacobian
        (
            subMats[1][1], &x[solidFirstEqnID]
        );

        // Form off-diagonal submatrices:
        //  - Afs: solid terms in the fluid equations
        //  - Asf: fluid terms (pressure on interface) in solid equations

        // Afs
        if (solidToFluidCoupling_) // && coupled())
        {
            formAfs(subMats[0][1], fluidBlockSize, solidBlockSize, twoD);
        }

        // Asf
        if (fluidToSolidCoupling_) // && coupled())
        {
            formAsf(subMats[1][0], fluidBlockSize, solidBlockSize, twoD);
        }

        // Scale the solid matrix to preserve the condition number of the
        // monolithic system
        // We must assembly the matrix before we can use MatScale
        // Complete matrix assembly
        CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

        // Asf
        MatScale(subMats[1][0], solidSystemScaleFactor_);

        // Ass
        MatScale(subMats[1][1], solidSystemScaleFactor_);

        // Zero coupling matrices if not coupled
        // There is an issue with non-zeros if we do not insert the coefficients in
        // the first step, so we do it and zero them if not coupled
        // Todo: find a more elegant/efficient solution
        if (!coupled())
        {
            MatZeroEntries(subMats[0][1]);
            MatZeroEntries(subMats[1][0]);
        }
    }
    else if (nRegions_ == 3)
    {
        // We will assembly the seven sub-matrices:

        // Diagonal submatrices:
        //  - Aff: fluid equations (momentum and pressure)
        //  - Amm: mesh motion equations in the fluid domain
        //  - Ass: solid equations
        //
        // Off-diagonal submatrices:
        //  - Afm: mesh motion terms in the fluid equations
        //  - Afs: solid terms in the fluid equations
        //  - Ams: solid terms (interface motion) in the mesh motion equation
        //  - Asf: fluid terms (pressure on interface) in solid equations

        // Zero entries
        MatZeroEntries(jac);

        // Get access to the sub-matrices
        PetscInt nr, nc;
        Mat **subMats;
        CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

        if (nr != 3 || nc != 3)
        {
            FatalErrorInFunction
                << "The matrix has the wrong number of sub matrices: "
                << "nr = " << nr << ", nc = " << nc << abort(FatalError);
        }

        // Set the motion block size
        const label motionBlockSize = solidBlockSize;

        // The scalar row at which the motion equations start
        const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

        // The scalar row at which the solid equations start
        const label solidFirstEqnID =
            motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;

        // Form diagonal submatrices
        //  - Aff: fluid equations (momentum and pressure)
        //  - Amm: mesh motion equations in the fluid domain
        //  - Ass: solid equations
        //

        // Aff
        refCast<foamPetscSnesHelper>(fluid()).formJacobian
        (
            subMats[0][0], x
        );

        // Amm
        refCast<foamPetscSnesHelper>(motionSolid()).formJacobian
        (
            subMats[1][1], &x[motionFirstEqnID]
        );

        // Ass
        refCast<foamPetscSnesHelper>(solid()).formJacobian
        (
            subMats[2][2], &x[solidFirstEqnID]
        );

        // Form off-diagonal submatrices:
        //  - Afm: mesh motion terms in the fluid equations
        //  - Afs: solid terms in the fluid equations
        //  - Ams: solid terms (interface motion) in the mesh motion equation
        //  - Asf: fluid terms (pressure on interface) in solid equations

        // Afm
        if (meshToFluidCoupling_) // && coupled())
        {
            formAfm(subMats[0][1], fluidBlockSize, motionBlockSize, twoD);
        }

        // Afs
        if (solidToFluidCoupling_) // && coupled())
        {
            formAfs(subMats[0][2], fluidBlockSize, solidBlockSize, twoD);
        }

        // Ams
        if (solidToMeshCoupling_) // && coupled())
        {
            formAms(subMats[1][2], solidBlockSize, motionBlockSize, twoD);
        }

        // Asf
        if (fluidToSolidCoupling_) // && coupled())
        {
            formAsf(subMats[2][0], fluidBlockSize, solidBlockSize, twoD);
        }


        // Scale the solid matrix to preserve the condition number of the
        // monolithic system
        // We must assembly the matrix before we can use MatScale
        // Complete matrix assembly
        CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

        // Amm
        MatScale(subMats[1][1], motionSystemScaleFactor_);

        // Ams
        MatScale(subMats[1][2], motionSystemScaleFactor_);

        // Asf
        MatScale(subMats[2][0], solidSystemScaleFactor_);

        // Ass
        MatScale(subMats[2][2], solidSystemScaleFactor_);

        // Zero coupling matrices if not coupled
        // There is an issue with non-zeros if we do not insert the coefficients in
        // the first step, so we do it and zero them if not coupled
        // Todo: find a more elegant/efficient solution
        if (!coupled())
        {
            MatZeroEntries(subMats[0][1]);
            MatZeroEntries(subMats[0][2]);
            MatZeroEntries(subMats[1][2]);
            MatZeroEntries(subMats[2][0]);
        }
    }
    else
    {
        FatalErrorInFunction
            << "nRegions must be 2 or 3!" << abort(FatalError);
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
