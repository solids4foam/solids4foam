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

#ifdef USE_PETSC

#ifdef OPENFOAM_COM

#include "newtonQuasiMonolithicCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "directMapInterfaceToInterfaceMapping.H"
#include "fixedValueFvPatchFields.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"
#include "newtonIcoFluid.H"
#include "nonLinGeomTotalLagVelocitySolid.H"
#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "meshMotionSolidModelFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newtonQuasiMonolithicCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, newtonQuasiMonolithicCouplingInterface, dictionary
);


// * * * * * * * * * * * Private Member Functions* * * * * * * * * * * * * * //


void newtonQuasiMonolithicCouplingInterface::makeIsFluidIsMotionIsSolid() const
{
    if
    (
        isFluid_ != nullptr
     || isFluidVelocity_ != nullptr
     || isFluidPressure_ != nullptr
     // || isMotion_ != nullptr
     || isSolid_ != nullptr
    )
    {
        FatalErrorInFunction
            << "Pointer already set" << exit(FatalError);
    }

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Fluid
    label NFluid = -1;
    {
        // Set the number local block equations
        const label blockn = fluid().mesh().nCells();

        // Fluid block size
        const label fluidBlockSize = twoD ? 3 : 4;

        // Set the number local scalar equations
        const label n = blockn*fluidBlockSize;

        // Set the global system size
        NFluid = returnReduce(n, sumOp<label>());

        // Create the index sets, where it is assumed the motion equations come first
        ISCreateStride(PETSC_COMM_WORLD, NFluid, 0, 1, &isFluid_);
    }

    // Fluid velocity field
    {
        const label blockn = fluid().mesh().nCells();
        const label fluidBlockSize = twoD ? 3 : 4;

        const label nVelPerCell = twoD ? 2 : 3;
        const label nVel = blockn*nVelPerCell;

        List<PetscInt> idx(nVel);

        PetscInt k = 0;
        for (label c = 0; c < blockn; ++c)
        {
            const PetscInt base = c*fluidBlockSize;
            idx[k++] = base + 0;
            idx[k++] = base + 1;
            if (!twoD)
            {
                idx[k++] = base + 2;
            }
        }

        ISCreateGeneral
        (
            PETSC_COMM_WORLD,
            idx.size(),
            idx.begin(),
            PETSC_COPY_VALUES,
            &isFluidVelocity_
        );
    }

    // Fluid pressure field
    {
        const label blockn = fluid().mesh().nCells();
        const label fluidBlockSize = twoD ? 3 : 4;

        const label nP = blockn;
        List<PetscInt> idx(nP);

        for (label c = 0; c < blockn; ++c)
        {
            const PetscInt base = c * fluidBlockSize;
            idx[c] = base + (fluidBlockSize - 1);
        }

        ISCreateGeneral
        (
            PETSC_COMM_WORLD,
            idx.size(),
            idx.begin(),
            PETSC_COPY_VALUES,
            &isFluidPressure_
        );
    }

    // Fluid mesh motion
    // label NMotion = -1;
    // {
    //     // Set the number local block equations
    //     // Note: the motion mesh is the fluid mesh
    //     const label blockn = fluid().mesh().nCells();

    //     // Motion and solid block size
    //     const label motionBlockSize = twoD ? 2 : 3;

    //     // Set the number local scalar equations
    //     const label n = blockn*motionBlockSize;

    //     // Set the global system size
    //     NMotion = returnReduce(n, sumOp<label>());

    //     // Create the index sets, where it is assumed the motion equations come first
    //     ISCreateStride(PETSC_COMM_WORLD, NMotion, NFluid, 1, &isMotion_);
    // }

    // Solid
    {
        // Set the number local block equations
        const label blockn = solid().mesh().nCells();

        // Solid block size
        const label solidBlockSize = twoD ? 2 : 3;

        // Set the number local scalar equations
        const label n = blockn*solidBlockSize;

        // Set the global system size
        const label NSolid = returnReduce(n, sumOp<label>());

        // Create the solid index set assuming the solid equations are after the
        // motion equations
        // ISCreateStride(PETSC_COMM_WORLD, NSolid, NMotion, 1, &isSolid_);
        ISCreateStride(PETSC_COMM_WORLD, NSolid, NFluid, 1, &isSolid_);
    }
}


foamPetscSnesHelper& newtonQuasiMonolithicCouplingInterface::motion()
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


solidModel& newtonQuasiMonolithicCouplingInterface::motionSolid()
{
    return refCast<solidModel>(motion());
}


void newtonQuasiMonolithicCouplingInterface::createSubMatsAndMat
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
    for (label r = 0; r < nRegions; ++r)
    {
        ISDestroy(&isRow[r]);
        ISDestroy(&isCol[r]);
    }
}


label newtonQuasiMonolithicCouplingInterface::initialiseAfs
(
    Mat Afs,
    const fvMesh& fluidMesh,
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

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Afs matrix is equal to the number of cells at the interface

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
    std::vector<label> d_nnz(scalarRowN, 0);
    std::vector<label> o_nnz(scalarRowN, 0);

    // Set non-zeros for each interface fluid cells
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const labelList& fluidFaceCells =
        fluidMesh.boundary()[fluidPatchID].faceCells();
    forAll(fluidFaceMap, fluidFaceI)
    {
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Calculate the row index for this cells first scalar equation
        label rowI = fluidCellID*fluidBlockSize;

        // Set the number of non-zeros to be the number of solid equations
        d_nnz[rowI++] = solidBlockSize;
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
        MatMPIAIJSetPreallocation(Afs, 0, d_nnz.data(), 0, o_nnz.data())
    );

    // Raise an error if mallocs are required during matrix assembly
    CHKERRQ(MatSetOption(Afs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    //CHKERRQ(MatSetOption(Afs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

    return 0;
}


label newtonQuasiMonolithicCouplingInterface::initialiseAsf
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
    std::vector<label> d_nnz(scalarRowN, 0);
    std::vector<label> o_nnz(scalarRowN, 0);

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


label newtonQuasiMonolithicCouplingInterface::formAfs
(
    Mat Afs,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
)
{
    // This function incorporates the effect of the solid interface velocity on
    // the fluid interface velocity field

    // The fluid and mesh motion refer to the same mesh so they can use the same
    // global cells object
    const globalIndex& globalCells =
        refCast<foamPetscSnesHelper>(fluid()).globalCells();


    // The interface velocity equal the mesh motion interface
    // velocity

    // The fluid interface is a prescribed velocity (fixedValue) condition
    // where we assume the fluid wall velocity is equal to the mesh velocity of
    // the adjacent cell centre. This approximatation is sufficiently
    // accurate as a preconditioner for the matrix and will not affected the
    // converged solution (which is entirely governed by formResidual)

    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
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

    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatch = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatch.faceCells();

    // Lookup the fluid patch
    const fvPatchVectorField& fluidPatchU =
        fluid().U().boundaryField()[fluidPatchID];

    // Check the fluid U patch type is fixedValue (not just derived from it)
    // TODO: use the corrected version of fixedValue
    if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
    //if (fluidPatchU.type() != "fixedValue")
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

    // Lookup the motion ddt scheme
    // NOTE: we should match the solid's ddt scheme
    const word solidDdtScheme =
        word(solidMesh().ddtScheme("ddt(" + solid().D().name() +')'));

    // The known fluid boundary face value is now replaced by the adjacent
    // solid cell velocity
    scalarField fluidPatchDiffusionCoeffs(fluidPatch.size(), 0.0);
    vectorField fluidPatchPressureCoeffs(fluidPatch.size(), vector::zero);

    // Diffusion coefficient is nu*|Sf|/(|n & d|) for the momentum equation
    fluidPatchDiffusionCoeffs =
        fluidPatchNuEff*fluidPatch.magSf()*fluidPatch.deltaCoeffs();

    // Pressure coefficient from div(U)
    fluidPatchPressureCoeffs = -fluidPatchSf;

    // First we will insert the contribution to the fluid momentum equation
    // coming from the diffusion term

    // Second we will insert the contribution to the fluid continuity
    // (pressure) equation, where the div(U) term should use the adjacent
    // solid cell velocity instead of the known boundary face velocity

    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    forAll(fluidPatch, fluidFaceI)
    {
        // Fluid cell adjacent to the interface
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Neighbouring solid face index
        // TODO: fix for parallel! use zone face ID instead of face ID
        const label solidFaceI = fluidFaceMap[fluidFaceI];

        // Neighbouring solid cell
        const label solidCellID = solidFaceCells[solidFaceI];

        // We will add coefficients at block row "fluidCellID" and block
        // column "solidCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of fluid matrix
        const label globalBlockRowI = globalCells.toGlobal(fluidCellID);

        // Global block column ID of solid matrix
        const label globalBlockColI = globalCells.toGlobal(solidCellID);

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*fluidBlockSize;
        label globalColI = globalBlockColI*solidBlockSize;

        // Momentum coefficient for this face
        PetscScalar value = fluidPatchDiffusionCoeffs[fluidFaceI];

        // Manually insert the 3 scalar coefficients (2 in 2-D) for the momentum
        // coupling
        // Note that the pressure equation coupling is inserted after (further
        // below in this function)
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
        value = fluidPatchPressureCoeffs[fluidFaceI][vector::X];

        globalRowI++; // pressure equation
        globalColI = globalBlockColI*solidBlockSize; // x velocity
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        value = fluidPatchPressureCoeffs[fluidFaceI][vector::Y];
        globalColI++; // y velocity
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            value = fluidPatchPressureCoeffs[fluidFaceI][vector::Z];
            globalColI++; // z velocity
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


label newtonQuasiMonolithicCouplingInterface::formAsf
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
    // accurate as a preconditioner for the matrix and will not affect the
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


void newtonQuasiMonolithicCouplingInterface::mapInterfaceMotionUToFluidU()
{
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

    // The mesh motion patch is the same as the fluid patch
    const label motionPatchID = fluidPatchID;

    // Map the motion interface velocity to the fluid interface
    const fvPatchVectorField& motionPatchU =
        motionSolid().U().boundaryField()[motionPatchID];

    forAll(fluidPatchU, fluidFaceI)
    {
        fluidPatchU[fluidFaceI] = motionPatchU[fluidFaceI];
    }
}


void newtonQuasiMonolithicCouplingInterface::mapInterfaceSolidToMeshMotion()
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

    // Lookup the mesh motion velocity field
    fvPatchVectorField& motionPatchU =
        motionSolid().U().boundaryFieldRef()[fluidPatchID];

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];

    // Map the solid interface displacement to the motion interface
    if (extrapolateSolidInterfaceDisplacement_)
    {
        const fvPatchVectorField& solidPatchD =
            solid().D().boundaryField()[solidPatchID];
        const fvPatchVectorField& solidPatchU =
            solid().U().boundaryField()[solidPatchID];

        forAll(motionPatchD, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Extrapolated patch value (larger stencil)
            motionPatchD[fluidFaceI] = solidPatchD[solidFaceID];
            motionPatchU[fluidFaceI] = solidPatchU[solidFaceID];
        }
    }
    else
    {
        const labelList& solidFaceCells =
            solidMesh().boundary()[solidPatchID].faceCells();
        const vectorField& solidDI = solid().D();
        const vectorField& solidUI = solid().U();

        forAll(motionPatchD, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Adjacent cell value
            const label solidCellID = solidFaceCells[solidFaceID];
            motionPatchD[fluidFaceI] = solidDI[solidCellID];
            motionPatchU[fluidFaceI] = solidUI[solidCellID];
        }
    }

    if (interfaceToInterfaceList().size() != 1)
    {
        FatalErrorInFunction
            << "Only one interface allowed" << abort(FatalError);
    }

    // fixedDisplacement patches overwrite the value with an internally stored
    // totalDisp field, so we must set that
    if (isA<fixedDisplacementFvPatchVectorField>(motionPatchD))
    {
        refCast<fixedDisplacementFvPatchVectorField>(motionPatchD).totalDisp() =
            motionPatchD;
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

    // Check the motion pointD interface patch type is fixed value
    if
    (
        !isA<fixedValuePointPatchVectorField>
        (
            motionSolid().pointD().boundaryFieldRef()[fluidPatchID]
        )
    )
    {
        FatalErrorInFunction
            << "The meshMotionFluid pointD interface patch must be of type "
            << "'fixedValue'" << abort(FatalError);
    }

    // Set the mesh interface pointD
    // Use "==" to reassign fixedValue
    motionSolid().pointD().boundaryFieldRef()[fluidPatchID] == meshPatchPointD;

    // Correct boundary conditions to enforce the new patch values on the
    // internal field
    motionSolid().pointD().correctBoundaryConditions();
}


// void newtonQuasiMonolithicCouplingInterface::resetFieldsToOldTime()
// {
//     Info<< "Resetting primary fields to their old time values" << nl << endl;

//     // Velocity
//     fluid().U() = fluid().U().oldTime();

//     // Pressure
//     fluid().p() = fluid().p().oldTime();

//     // Displacement
//     solid().D() = solid().D().oldTime();

//     // Mesh motion
//     motionSolid().D() = motionSolid().D().oldTime();

//     // Update the interface fields
//     mapInterfaceSolidToMeshMotion();
//     mapInterfaceMotionUToFluidU();

//     // Insert the OpenFOAM fields into the PETSc solution vector

//     // Set twoD flag
//     const bool twoD = fluid().twoD();

//     // Set fluid and solid block sizes
//     const label fluidBlockSize = twoD ? 3 : 4;
//     const label solidBlockSize = twoD ? 2 : 3;

//     // The scalar row at which the solid equations start
//     const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

//     // The scalar row at which the motion equations start
//     const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

//     // Access the raw solution data
//     PetscScalar *xx;
//     VecGetArray(foamPetscSnesHelper::solution(), &xx);

//     // Insert the fluid velocity
//     foamPetscSnesHelper::InsertFieldComponents
//     (
//         fluid().U().primitiveFieldRef(),
//         xx,
//         0, // Location of U
//         fluidBlockSize,
//         twoD ? labelList({0,1}) : labelList({0,1,2})
//     );

//     // Insert the fluid pressure
//     foamPetscSnesHelper::InsertFieldComponents
//     (
//         fluid().p().primitiveFieldRef(),
//         xx,
//         fluidBlockSize - 1, // Location of p component
//         fluidBlockSize
//     );

//     // Insert the displacement
//     foamPetscSnesHelper::InsertFieldComponents
//     (
//         solid().D().primitiveFieldRef(),
//         &xx[solidFirstEqnID],
//         0, // Location of first component
//         solidBlockSize,
//         twoD ? labelList({0,1}) : labelList({0,1,2})
//     );

//     // Insert the motion displacement
//     foamPetscSnesHelper::InsertFieldComponents
//     (
//         motionSolid().D().primitiveFieldRef(),
//         &xx[motionFirstEqnID],
//         0, // Location of first component
//         solidBlockSize,
//         twoD ? labelList({0,1}) : labelList({0,1,2})
//     );

//     // Restore the solution vector
//     VecRestoreArray(foamPetscSnesHelper::solution(), &xx);
// }

void newtonQuasiMonolithicCouplingInterface::retrieveSolution
(
    Vec solution,
    volVectorField& UFluid,
    volScalarField& p,
    volVectorField& USolid,
    const label fluidBlockSize,
    const label solidBlockSize,
    const bool twoD
) const
{
    // Access the raw solution data
    const PetscScalar *xx;
    VecGetArrayRead(solution, &xx);

    // Retrieve the fluid velocity
    foamPetscSnesHelper::ExtractFieldComponents
    (
        xx,
        UFluid.primitiveFieldRef(),
        0, // Location of U
        fluidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );
    UFluid.correctBoundaryConditions();

    // Retrieve the fluid pressure
    foamPetscSnesHelper::ExtractFieldComponents
    (
        xx,
        p.primitiveFieldRef(),
        fluidBlockSize - 1, // Location of p component
        fluidBlockSize
    );
    p.correctBoundaryConditions();

    // Retrieve the solid velocity
    const label solidFirstEqnID = UFluid.size()*fluidBlockSize;
    foamPetscSnesHelper::ExtractFieldComponents
    (
        &xx[solidFirstEqnID],
        USolid.primitiveFieldRef(),
        0, // Location of first component
        solidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );
    USolid.correctBoundaryConditions();

    // Restore the x vector
    VecRestoreArrayRead(solution, &xx);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonQuasiMonolithicCouplingInterface::newtonQuasiMonolithicCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "UpU",
        fileName
        (
            fsiProperties().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        fluid().mesh(), // Should not be used
        solutionLocation::CELLS, // Used by helper functions
        fsiProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true // Will PETSc be used
    ),
    fluidSystemScaleFactor_
    (
        fsiProperties().lookupOrDefault<scalar>
        (
            "fluidSystemScaleFactor",
            1.0 //gAverage(1.0/fluid().nf().primitiveField())
        )
    ),
    solidSystemScaleFactor_
    (
        fsiProperties().lookupOrDefault<scalar>
        (
            "solidSystemScaleFactor",
            gAverage
            (
                1.0
               /(
                   solid().mechanical().shearModulus()().primitiveField()
                  *runTime.deltaTValue()
                )
            )
        )
    ),
    // fluidToSolidCoupling_(fsiProperties().lookup("fluidToSolidCoupling")),
    // meshToFluidCoupling_(fsiProperties().lookup("meshToFluidCoupling")),
    // solidToMeshCoupling_(fsiProperties().lookup("solidToMeshCoupling")),
    extrapolateSolidInterfaceDisplacement_
    (
        fsiProperties().lookupOrDefault<Switch>
        (
            "extrapolateSolidInterfaceDisplacement",
            true
        )
    ),
    passViscousStress_(fsiProperties().lookup("passViscousStress")),
    nRegions_(2),
    subMatsPtr_(nullptr),
    isFluid_(nullptr),
    isFluidVelocity_(nullptr),
    isFluidPressure_(nullptr),
    isSolid_(nullptr),
    configuredNestedFluidSplit_(false),
    tsLogPtr_()
    // oldTimeValue_(runTime.value()),
    // nConsecutiveFailedSolves_(0),
    // maxAllowedConsecutiveFailedSolves_
    // (
    //     fsiProperties().lookupOrDefault<label>
    //     (
    //         "maxAllowedConsecutiveFailedSolves", 5
    //     )
    // )
{
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

    Info<< "fluidSystemScaleFactor = " << fluidSystemScaleFactor_ << nl
        << "solidSystemScaleFactor = " << solidSystemScaleFactor_ << endl;

    // Store old time values
    fluid().U().storeOldTime();
    fluid().U().oldTime().storeOldTime();
    fluid().U().oldTime().oldTime().storeOldTime();
    fluid().p().storeOldTime();
    solid().U().storeOldTime();
    solid().U().oldTime().storeOldTime();
    solid().D().storeOldTime();
    motionSolid().D().storeOldTime();
    solid().sigma().storeOldTime();

    // Info<< "fluidToSolidCoupling: " << fluidToSolidCoupling_ << nl
    //     << "meshToFluidCoupling: " << meshToFluidCoupling_ << nl
    //     << "solidToMeshCoupling: " << solidToMeshCoupling_ << nl
    //     << "extrapolateSolidInterfaceDisplacement: "
    //     << extrapolateSolidInterfaceDisplacement_ << nl
    //     << "passViscousStress: " << passViscousStress_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void newtonQuasiMonolithicCouplingInterface::setDeltaT(Time& runTime)
{
    if
    (
        runTime.controlDict().getOrDefault("adjustTimeStep", false)
     && foamPetscSnesHelper::snesHasRun()
    )
    {
        const scalar maxDeltaT =
            runTime.controlDict().get<scalar>("maxDeltaT");
        const scalar minDeltaT =
            runTime.controlDict().get<scalar>("minDeltaT");

        const int minTargetNIter =
            runTime.controlDict().getOrDefault<int>("minTargetNIter", 3);
        const int maxTargetNIter =
            runTime.controlDict().getOrDefault<int>("maxTargetNIter", 6);

        const Switch enableTimeStepLog =
            runTime.controlDict().getOrDefault("logTimeStepAdjustments", true);

        PetscInt numIter;
        SNESGetIterationNumber(foamPetscSnesHelper::snes(), &numIter);

        SNESConvergedReason reason;
        SNESGetConvergedReason(foamPetscSnesHelper::snes(), &reason);

        const scalar currentDeltaT = runTime.deltaTValue();
        scalar newDeltaT = currentDeltaT;

        // if (reason == SNES_DIVERGED_FUNCTION_DOMAIN)
        if (reason < 0)
        {
            // SNES failed to converge
            newDeltaT = max(0.25*currentDeltaT, minDeltaT);
            Info<< nl << "SNES failed to converge: "
                << "reducing timestep to " << newDeltaT << endl;
        }
        else
        {
            // Guard against zero
            if (numIter <= 0)
            {
                numIter = 1;
            }

            scalar factor = 1.0;

            if (numIter > maxTargetNIter + 1)
            {
                factor = max(0.5, 0.9*scalar(maxTargetNIter)/numIter);
            }
            else if (numIter < minTargetNIter - 1)
            {
                factor = min(1.5, 1.1*scalar(maxTargetNIter)/numIter);
            }

            newDeltaT = clamp(factor*currentDeltaT, minDeltaT, maxDeltaT);
        }

        Info<< "Nonlinear iterations = " << numIter << nl
            << "Old time step        = " << currentDeltaT << nl
            << "New time step        = " << newDeltaT << nl << endl;

        runTime.setDeltaT(newDeltaT);

        if (enableTimeStepLog)
        {
            if (tsLogPtr_.empty())
            {
                const fileName timeStepLogFile =
                    runTime.controlDict().getOrDefault<fileName>
                    (
                        "timeStepLogFile", "timeStepLog.dat"
                    );

                tsLogPtr_.set(new OFstream(timeStepLogFile));

                tsLogPtr_()
                    << "Time currentDeltaT newDeltaT numIter reason" << endl;
            }

            tsLogPtr_()
                << runTime.timeName() << " "
                << currentDeltaT << " "
                << newDeltaT << " "
                << numIter << " "
                << reason << endl;
        }
    }
}


bool newtonQuasiMonolithicCouplingInterface::evolve()
{
    // Algorithm 2 from Jaiman and Joshi, adapted to a finite volume framework
    //
    // Time advancement from t^(n-1) → t^n
    //
    // 1. Start from known solutions at previous time steps:
    //      UFluid[n-1], UFluid[n-2], USolid[n-1], USolid[n-2], and D[n-1]
    //    where
    //      UFluid : fluid velocity
    //      USolid : solid velocity
    //      D      : solid displacement  (or we could have used position)
    //
    // 2. Advance to the new time level t^n:
    //
    //    (a) Explicitly update the solid displacement field D[n]:
    //          D[n] = D[n-1] + (3*dt/2)*USolid[n-1] - (dt/2)*USolid[n-2]
    //
    //    (b) Determine new positions of the vertices on the fluid domain
    //        mesh given D[n] on the fluid-solid interface, where D is known on
    //        the interface from the step (a)
    //
    //    (c) Update the fluid mesh motion based on the result from (b),
    //        i.e. fluidMesh.update()
    //
    //    (d) Evaluate the predicted fluid fluid (phi) and fluid mesh fluid
    //        (phiMotion): these are used for the linearised fluid convection
    //        term
    //
    //    (e) Compute the element-level stabilisation parameters : this step is
    //        not required in our finite volume implementation. Instead, we will
    //        incorporate typical finite volume stabilisation terms (e.g.,
    //        Rhie-Chow) when forming the monolithic system in step (f)
    //
    //    (f) Solve the coupled monolithic system for UFluid[n], pFluid[n], and
    //        USolid[n] using the fluid mesh from step (c) and the linear
    //        linearised mesh flux from step (d)
    //
    // Notes:
    //  - This formulation is second-order accurate in time and fully stabilised.
    //  - The algorithm couples the fluid and solid fields in a quasi-monolithic
    //    manner, ensuring consistent interface kinematics and dynamics.

    // Check if coupling switch needs to be updated
    if (!coupled())
    {
        updateCoupled();
    }

    // Preliminaries
    volVectorField& UFluid = fluid().U();
    surfaceScalarField& phi = fluid().phi();
    volScalarField& p = fluid().p();
    volVectorField& USolid = solid().U();
    volVectorField& D = solid().D();
    volVectorField& DMotion = motionSolid().D();
    const dimensionedScalar& deltaT = runTime().deltaT();
    const bool twoD = fluid().twoD();
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // Ensure boundary conditions are up-to-date
    UFluid.correctBoundaryConditions();
    p.correctBoundaryConditions();
    USolid.correctBoundaryConditions();
    D.correctBoundaryConditions();
    DMotion.correctBoundaryConditions();

    //    (a) Explicitly update the solid displacement field D[n]:
    //          D[n] = D[n-1] + (3*dt/2)*USolid[n-1] - (dt/2)*USolid[n-2]
    D = D.oldTime()
      + (3.0*deltaT/2.0)*USolid.oldTime()
      - (deltaT/2.0)*USolid.oldTime().oldTime();

    // Update displacement gradient
    solid().mechanical().grad(D, solid().gradD());

    // Interpolate cell displacements to vertices
    solid().mechanical().interpolate
    (
        D, solid().gradD(), solid().pointD()
    );
    solid().pointD().correctBoundaryConditions();


    //    (b) Determine new positions of the vertices on the fluid domain
    //        mesh given D[n] on the fluid-solid interface, where D is known on
    //        the interface from the step (a)

    // First, we will map DSolid to DMotion at the interface
    // This function maps DSolid to DMotion, pointDSolid to pointDMotion,
    // and USolid to UFluid at the interface
    // NOTE: mapping of USolid to UFluid is not needed here, but I already had
    // this function which did this - todo: we can remove that step
    mapInterfaceSolidToMeshMotion();


    //    (c) Update the fluid mesh motion based on the result from (b),
    //        i.e. fluidMesh.update()
    // This updates DMotion based on the interface D which was set in step (b)
    // In addition, this function moves the fluid mesh using this D field
    fluidMesh().update();

    // Print the mesh Courant number
    {
        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(fluidMesh().phi()))().primitiveField()
        );

        const scalar meshCoNum =
            0.5*gMax(sumPhi/fluidMesh().V().field())*deltaT.value();

        const scalar meanMeshCoNum =
            0.5*(gSum(sumPhi)/gSum(fluidMesh().V().field()))*deltaT.value();

        Info<< "Fluid mesh Courant number mean: " << meanMeshCoNum
            << " max: " << meshCoNum << endl;
    }


    //    (d) Evaluate the predicted fluid fluid (phi) and fluid mesh fluid
    //        (phiMotion): these are used for the linearised fluid convection
    //        term

    // Extrapolate the fluid velocity
    if
    (
        Switch
        (
            fluid().fluidProperties().lookup("fluidFluxExtrapolationAlgorithm1")
        )
    )
    {
        // Equation 6.10
        phi = fvc::interpolate
              (
                  2.0*UFluid.oldTime() - UFluid.oldTime().oldTime()
              ) & fluidMesh().Sf();
    }
    else
    {
        // Equation 6.30
        phi = fvc::interpolate
              (
                  2.25*UFluid.oldTime()
                - 1.5*UFluid.oldTime().oldTime()
                + 0.25*UFluid.oldTime().oldTime().oldTime()
              ) & fluidMesh().Sf();
    }

    // Print the Courant number
    {
        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        const scalar CoNum =
            0.5*gMax(sumPhi/fluidMesh().V().field())*deltaT.value();

        const scalar meanCoNum =
            0.5*(gSum(sumPhi)/gSum(fluidMesh().V().field()))*deltaT.value();

        Info<< "Fluid Courant number mean: " << meanCoNum
            << " max: " << CoNum << endl;
    }

    // In Jamain and Joshi, they update the motionU (fluid mesh velocity) based
    // on the current, old and old-old fluid mesh positions (or equivalently the
    // displacements) using a 2nd order method. These mesh velocities are used
    // in the linearised convection term.
    // In OpenFOAM, mesh.phi() will already return a 2nd order approximation of
    // the fluid mesh fluid (velocity times area at the faces), which is also
    // consistent with the geometric conservation law, assuming we use a 2nd
    // order time discretisation.
    // So we don't need to do anything here; we just need to use mesh.phi()


    //    (e) Compute the element-level stabilisation parameters : this step is
    //        not required in our finite volume implementation. Instead, we will
    //        incorporate typical finite volume stabilisation terms (e.g.,
    //        Rhie-Chow) when forming the monolithic system in step (f)
    //
    // Nothing to do here.


    //    (f) Solve the coupled monolithic system for UFluid[n], pFluid[n], and
    //        USolid[n] using the fluid mesh from step (c) and the linear
    //        linearised mesh flux from step (d)

    // We will call PETSc SNES to form and solve the coupled monolithic system
    // Note that this will call formResidual, formJacobian, initialiseResidual,
    // initialiseJacobian as defined above
    // This could system will be linear so the PETSc option should be selected
    // with this in mind
    Info<< "Solving the monolithic momentum-continuity system for Up using "
        << "PETSc SNES" << endl;
    foamPetscSnesHelper::solve();

    // Retrieve the solution from PETSc and copy into UFluid, p, and USolid
    retrieveSolution
    (
        foamPetscSnesHelper::solution(),
        UFluid, p, USolid,
        fluidBlockSize, solidBlockSize,
        twoD
    );

    return 0;
}


label newtonQuasiMonolithicCouplingInterface::initialiseJacobian(Mat& jac)
{
    // The monolithic fluid-solid interaction system can be expressed in
    // linearised Ax=b form, where the matrix A is
    //
    //     *-----------*
    //     | Aff | Afs |
    //     |-----------|
    //     | Asf | Ass |
    //     *-----------*
    //
    // The diagonal submatrices are
    //  - Aff: fluid equations (momentum and pressure)
    //  - Ass: solid equations (momentum)
    //
    // The off diagonal submatrices, which represent coupling between
    // equations, are
    //  - Afs: solid terms (interface velocity) in the fluid equations
    //  - Asf: fluid terms (pressure on interface) in solid equation

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // For brevity and convenience, we will store the size and blockSize of the
    // regions in a labelPairList
    const labelPairList nBlocksAndBlockSize
    (
        {
            {fluidMesh().nCells(), fluidBlockSize},
            {solidMesh().nCells(), solidBlockSize}
        }
    );

    // Set the submatrix sizes
    createSubMatsAndMat(jac, subMatsPtr_, nBlocksAndBlockSize);

    // Get access to the sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    // Initialise the diagonal submatrices:
    //  - Aff: fluid equations (momentum and pressure)
    //  - Ass: solid equations (momentum)

    // Aff
    foamPetscSnesHelper::initialiseJacobian
    (
        subMats[0][0], fluid().mesh(), fluidBlockSize, false
    );
    PetscObjectSetName((PetscObject)subMats[0][0], "Aff");

    // Ass
    foamPetscSnesHelper::initialiseJacobian
    (
        subMats[1][1], solid().mesh(), solidBlockSize, false
    );
    PetscObjectSetName((PetscObject)subMats[1][1], "Ass");

    // Initialise the diagonal submatrices:
    //  - Afs: solid velocity terms (interface) in the fluid equations
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

    return 0;
}


label newtonQuasiMonolithicCouplingInterface::initialiseSolution(Vec& x)
{
    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // Count the number of unknowns
    label n = 0;
    label N = 0;

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


void newtonQuasiMonolithicCouplingInterface::customiseSolver()
{
    KSP ksp; PC pc; const char* pct = nullptr;
    SNESGetKSP(snes(), &ksp);
    KSPGetPC(ksp, &pc);
    PCGetType(pc, &pct);

    // If a split preconditioner is used, let is know where the fluid unknowns
    // and solid unknowns are located
    // CHECK: is it ok to only call this for timeIndex == 1?
    if (pct && !std::strcmp(pct, PCFIELDSPLIT) && runTime().timeIndex() == 1)
    {
        Info<< "    Defining the fluid-solid field indices for the fluid-solid"
            << nl << "    split preconditioner" << endl;
        PCFieldSplitSetIS(pc, "fluid", isFluid());
        PCFieldSplitSetIS(pc, "solid", isSolid());
    }
}


label newtonQuasiMonolithicCouplingInterface::formResidual
(
    Vec f,
    const Vec x
)
{
    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // Preliminaries
    volVectorField& UFluid = fluid().U();
    volScalarField& p = fluid().p();
    volVectorField& USolid = solid().U();
    const label fluidSystemEnd = fluidMesh().nCells()*fluidBlockSize;
    const label solidFirstEqnID = fluidSystemEnd;


    // Considerations on the order of updating the fluid, solid, and mesh motion
    // residuals
    //  - solid depends on the fluid interface traction
    //  - fluid depends on the solid interface velocity
    //
    // Approach
    // 1. Retrieve UFluid, p and USolid from the PETSc solution vector
    // 2. Map the fluid interface traction to the solid interface
    // 3. Update the solid residual, which now has the correct interface
    //    traction, as a function of the current solid velocity (we do not yet
    //    have displacement)
    // 4. Map the solid interface velocity to the fluid interface
    // 5. Update the fluid residual, which now has the correct interface
    //    velocity, using the current (linearised) mesh motion
    // 6. Apply scaling factors to solid and fluid equations to preserve the
    //    condition number of the monolithic system


    // 1. Retrieve UFluid, p and USolid from the PETSc solution vector
    retrieveSolution
    (
        x,
        UFluid, p, USolid,
        fluidBlockSize, solidBlockSize,
        twoD
    );

    // 1b. When passViscousStress is enabled, update the fluid interface
    //     velocity from the current USolid BEFORE computing the viscous
    //     traction. Without this, patchViscousForce computes grad(U) using the
    //     stale fixedValue boundary from the previous SNES iteration, making
    //     the residual an inconsistent function of x.
    if (passViscousStress_)
    {
        mapInterfaceSolidToMeshMotion();
        mapInterfaceMotionUToFluidU();
    }

    // 2. Map the fluid interface traction to the solid interface

    // Fluid interface traction
    const label fluidPatchID =
        fluidSolidInterface::fluidPatchIndices()[0];

    // Fluid interface unit normals
    const vectorField fluidNf
    (
        fluidMesh().boundary()[fluidPatchID].nf()
    );

    // Fluid interface traction
    vectorField fluidTraction
    (
      - fluid().patchPressureForce(fluidPatchID)*fluidNf
    );

    if (passViscousStress_)
    {
        fluidTraction +=
            fluid().patchViscousForce(fluidPatchID);
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


    // 3. Update the solid residual, which now has the correct interface
    //    traction, as a function of the current solid velocity (we do not yet
    //    have displacement)

    // We must use a velocity formulation
    if (!isA<solidModels::nonLinGeomTotalLagVelocitySolid>(solid()))
    {
        FatalErrorInFunction
            << "Currently, the solid model must be of type "
            << solidModels::nonLinGeomTotalLagVelocitySolid::typeName
            << abort(FatalError);
    }

    // We create a temporary no-copy xSolid and fSolid Vec pointers, which are a
    // "view" of the solid equations in the full solution vectors x and f
    {
        Vec xSolid = nullptr;
        Vec fSolid = nullptr;
        VecGetSubVector(x, isSolid(), &xSolid);
        VecGetSubVector(f, isSolid(), &fSolid);
        VecSetBlockSize(xSolid, solidBlockSize);
        VecSetBlockSize(fSolid, solidBlockSize);
        refCast<foamPetscSnesHelper>(solid()).formResidual
        (
            fSolid, xSolid
        );
        VecRestoreSubVector(x, isSolid(), &xSolid);
        VecRestoreSubVector(f, isSolid(), &fSolid);
    }


    // 4. Map the solid interface velocity to the fluid interface via motionU
    mapInterfaceSolidToMeshMotion();
    mapInterfaceMotionUToFluidU();


    // 5. Update the fluid residual, which now has the correct interface
    //    velocity, using the current (linearised) mesh motion

    // We create a temporary no-copy xFluid and fFluid Vec pointers, which are a
    // "view" of the fluid equations in the full solution vectors x and f
    {
        Vec xFluid = nullptr;
        Vec fFluid = nullptr;
        VecGetSubVector(x, isFluid(), &xFluid);
        VecGetSubVector(f, isFluid(), &fFluid);
        VecSetBlockSize(xFluid, fluidBlockSize);
        VecSetBlockSize(fFluid, fluidBlockSize);
        refCast<fluidModels::newtonIcoFluid>(fluid()).formResidual
        (
            fFluid, xFluid, true // extrapolaed flux
        );
        VecRestoreSubVector(x, isFluid(), &xFluid);
        VecRestoreSubVector(f, isFluid(), &fFluid);
    }


    // 6. Apply scaling factors to solid and fluid equations to preserve the
    //    condition number of the monolithic system
    {
        PetscScalar* ff;
        VecGetArray(f, &ff);

        for (int i = 0; i < fluidSystemEnd; ++i)
        {
            ff[i] *= fluidSystemScaleFactor_;
        }

        const label solidSystemEnd =
            solidFirstEqnID + solidMesh().nCells()*solidBlockSize;
        for (int i = solidFirstEqnID; i < solidSystemEnd; ++i)
        {
            ff[i] *= solidSystemScaleFactor_;
        }

        VecRestoreArray(f, &ff);
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


label newtonQuasiMonolithicCouplingInterface::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Our monolithic system matrix will take the form:
    //
    //     *-----------*
    //     | Aff | Afs |
    //     |-----------|
    //     | Asf | Ass |
    //     *-----------*
    //

    // We will assembly the four sub-matrices:

    // Diagonal submatrices:
    //  - Aff: fluid equations (momentum and pressure)
    //  - Ass: solid equations (momentum)
    //
    // Off-diagonal submatrices:
    //  - Afs: solid motion terms in the fluid equations: solid velocity at the
    //         interface
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

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

    // Zero the entries
    MatZeroEntries(jac);

    // Form diagonal submatrices
    //  - Aff: fluid equations (momentum and pressure)
    //  - Ass: solid equations (momentum)

    // Aff
    {
        Vec xFluid = nullptr;
        VecGetSubVector(x, isFluid(), &xFluid);
        VecSetBlockSize(xFluid, fluidBlockSize);
        refCast<fluidModels::newtonIcoFluid>(fluid()).formJacobian
        (
            subMats[0][0], xFluid, true // use explicit flux
        );
        VecRestoreSubVector(x, isFluid(), &xFluid);
    }

    // Ass
    {
        Vec xSolid = nullptr;
        VecGetSubVector(x, isSolid(), &xSolid);
        VecSetBlockSize(xSolid, solidBlockSize);
        refCast<foamPetscSnesHelper>(solid()).formJacobian
        (
            subMats[1][1], xSolid
        );
        VecRestoreSubVector(x, isSolid(), &xSolid);
    }

    // Form off-diagonal submatrices:
    //  - Afs: solid terms in the fluid equations
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Afs
    formAfs(subMats[0][1], fluidBlockSize, solidBlockSize, twoD);

    // Asf
    formAsf(subMats[1][0], fluidBlockSize, solidBlockSize, twoD);


    // Scale the su-matrices to preserve the condition number of the
    // monolithic system

    // We must assembly the matrix before we can use MatScale
    // Complete matrix assembly
    CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

    // Aff
    MatScale(subMats[0][0], fluidSystemScaleFactor_);

    // Afs
    MatScale(subMats[0][1], fluidSystemScaleFactor_);

    // Asf
    MatScale(subMats[1][0], solidSystemScaleFactor_);

    // Ass
    MatScale(subMats[1][1], solidSystemScaleFactor_);

    // Configure the field split for velocity and pressure in the fluid
    if (!configuredNestedFluidSplit_)
    {
        KSP ksp;
        PC pc;
        const char* pct = nullptr;

        SNESGetKSP(snes(), &ksp);
        KSPGetPC(ksp, &pc);

        SNESGetKSP(snes(), &ksp);
        KSPGetPC(ksp, &pc);
        PCGetType(pc, &pct);

        // Only proceed if the top-level PC is a field split
        if (pct && std::strcmp(pct, PCFIELDSPLIT) == 0)
        {
            configuredNestedFluidSplit_ = true;

            // Ensure PC/KSP is set up *now that J is assembled*
            KSPSetUp(ksp);

            PetscInt nsplit;
            KSP* subksps;
            PCFieldSplitGetSubKSP(pc, &nsplit, &subksps);

            for (PetscInt i = 0; i < nsplit; ++i)
            {
                PC subpc;
                const char* subpct = nullptr;
                const char* prefix = nullptr;
                KSPGetPC(subksps[i], &subpc);
                PCGetType(subpc, &subpct);
                PCGetOptionsPrefix(subpc, &prefix);

                if (prefix && std::strstr(prefix, "fluid"))
                {
                    if (subpct && !std::strcmp(subpct, PCFIELDSPLIT))
                    {
                        Info<< "    Defining the fluid velocity-pressure field indices"
                            << "for the fluid sub-problem" << nl
                            << "    split preconditioner" << endl;

                        // local-to-fluid indices (0..NFluid-1)

                        PCFieldSplitSetIS(subpc, "velocity", isFluidVelocity());
                        PCFieldSplitSetIS(subpc, "pressure", isFluidPressure());

                        // ISDestroy(&isV);
                        // ISDestroy(&isP);

                        break;
                    }
                }
            }
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

#endif // OPENFOAM_NOT_EXTEND

#endif // ifdef USE_PETSC

// ************************************************************************* //
