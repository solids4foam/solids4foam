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
#include "linGeomVelocitySolid.H"
#include "nonLinGeomTotalLagVelocitySolid.H"
#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "meshMotionSolidModelFvMotionSolver.H"
#include "globalIndex.H"
#include "compatibilityFunctions.H"

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

namespace
{

PetscBool petscOptionEnabled(const char* name)
{
    PetscBool enabled = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, name, &enabled, nullptr);

    return enabled;
}


class retryTimeStateBuilder
:
    public Foam::TimeState
{
public:

    static Foam::TimeState rollbackState(const Foam::Time& runTime)
    {
        retryTimeStateBuilder state;

        static_cast<Foam::TimeState&>(state) =
            static_cast<const Foam::TimeState&>(runTime);

        // Restore the previous successful time-step length as the "saved"
        // value so the retry does not treat the failed step as the last
        // accepted one when operator++() updates deltaT0.
        state.deltaTSave_ = runTime.deltaT0Value();
        state.writeTime_ = false;

        return state;
    }
};

}


// * * * * * * * * * * * Private Member Functions* * * * * * * * * * * * * * //


void newtonQuasiMonolithicCouplingInterface::makeIsFluidIsMotionIsSolid() const
{
    if
    (
        isFluid_ != nullptr
     || isFluidVelocity_ != nullptr
     || isFluidPressure_ != nullptr
     || isSolid_ != nullptr
    )
    {
        FatalErrorInFunction
            << "Pointer already set" << exit(FatalError);
    }

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Fluid block parameters
    const label fluidBlockn = fluid().mesh().nCells();
    const label fluidBlockSize = twoD ? 3 : 4;
    const label fluidLocalSize = fluidBlockn*fluidBlockSize;

    // Solid block parameters (needed early for global start calculation)
    const label solidBlockn = solid().mesh().nCells();
    const label solidBlockSize = twoD ? 2 : 3;
    const label solidLocalSize = solidBlockn*solidBlockSize;

    // The Vec layout is per-process interleaved: each process owns
    // [proc_fluid, proc_solid] as a contiguous chunk. We need the IS
    // global indices to match this layout, NOT block-wise numbering.
    const label totalLocalSize = fluidLocalSize + solidLocalSize;
    const globalIndex vecGI(totalLocalSize);
    const label vecGlobalStart = vecGI.localStart();
    const label localFluidStart = vecGlobalStart;

    // Fluid IS: each process owns its local portion of the fluid indices
    ISCreateStride
    (
        PETSC_COMM_WORLD, fluidLocalSize, localFluidStart, 1, &isFluid_
    );

    // Fluid velocity field IS: velocity DOF indices within the fluid block
    {
        const label nVelPerCell = twoD ? 2 : 3;
        const label nVel = fluidBlockn*nVelPerCell;

        List<PetscInt> idx(nVel);

        PetscInt k = 0;
        for (label c = 0; c < fluidBlockn; ++c)
        {
            const PetscInt base = localFluidStart + c*fluidBlockSize;
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

    // Fluid pressure field IS: pressure DOF indices within the fluid block
    {
        List<PetscInt> idx(fluidBlockn);

        for (label c = 0; c < fluidBlockn; ++c)
        {
            const PetscInt base = localFluidStart + c*fluidBlockSize;
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

    // Solid IS: each process owns its local portion, offset by the local
    // fluid size within this process's chunk of the Vec
    {
        ISCreateStride
        (
            PETSC_COMM_WORLD,
            solidLocalSize,
            vecGlobalStart + fluidLocalSize,
            1,
            &isSolid_
        );
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
    // located.
    // Each process must specify its own global index range within each
    // region, accounting for DOFs on preceding processes (parallel fix).
    std::vector<IS> isRow(nRegions), isCol(nRegions);

    // The Vec layout is per-process interleaved: each process owns
    // [region0_local, region1_local, ...] as a contiguous chunk.
    // Compute total local size to get the Vec global start for this process.
    label totalLocalSize = 0;
    for (label r = 0; r < nRegions; ++r)
    {
        totalLocalSize +=
            nBlocksAndBlockSize[r].first()*nBlocksAndBlockSize[r].second();
    }
    const globalIndex vecGI(totalLocalSize);
    const label vecGlobalStart = vecGI.localStart();

    label cumulativeLocalSize = 0;
    for (label r = 0; r < nRegions; ++r)
    {
        const label nBlocks = nBlocksAndBlockSize[r].first();
        const label blockSize = nBlocksAndBlockSize[r].second();
        const label regionLocalSize = nBlocks*blockSize;

        // Global start index for this process's region r DOFs
        const label first = vecGlobalStart + cumulativeLocalSize;

        // Create an IS that covers the indices for region r
        ISCreateStride
        (
            PETSC_COMM_WORLD, regionLocalSize, first, 1, &isRow[r]
        );
        ISCreateStride
        (
            PETSC_COMM_WORLD, regionLocalSize, first, 1, &isCol[r]
        );
        cumulativeLocalSize += regionLocalSize;
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

    // Set non-zeros for each interface fluid cell
    // In parallel, we need to determine if the coupled solid cell is on
    // this processor (d_nnz) or another processor (o_nnz)
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatchRef = fluidMesh.boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatchRef.faceCells();

    // Build zone-level global solid cell IDs
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();

    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatchRef = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatchRef.faceCells();

    const globalIndex& solidGlobalCells =
        refCast<const foamPetscSnesHelper>(solid()).globalCells();

    // Local field: global cell ID for each local solid interface face
    scalarField localSolidGlobalCellIDs(solidPatchRef.size(), 0);
    forAll(solidPatchRef, faceI)
    {
        localSolidGlobalCellIDs[faceI] =
            solidGlobalCells.toGlobal(solidFaceCells[faceI]);
    }

    // Broadcast to zone level (all processors get the full array)
    const scalarField zoneSolidGlobalCellIDs
    (
        solidGlobalPatch.patchFaceToGlobal(localSolidGlobalCellIDs)
    );

    // Zone-level face map: fluid zone face -> solid zone face
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();
    const labelList& fluidFaceToZone = fluidGlobalPatch.faceToGlobalAddr();

    // Loop over LOCAL fluid interface faces
    forAll(fluidPatchRef, fluidFaceI)
    {
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Map local fluid face -> zone fluid face -> zone solid face
        const label fluidZoneFaceI = fluidFaceToZone[fluidFaceI];
        const label solidZoneFaceI = fluidFaceMap[fluidZoneFaceI];
        const label globalSolidCellID =
            label(zoneSolidGlobalCellIDs[solidZoneFaceI]);

        // Check if the coupled solid cell is on this processor
        const bool solidCellLocal =
            solidGlobalCells.isLocal(globalSolidCellID);

        // Calculate the row index for this cells first scalar equation
        label rowI = fluidCellID*fluidBlockSize;

        // Accumulate the number of non-zeros (a corner cell may have
        // multiple faces on the interface, each coupling to a different
        // solid cell)
        if (solidCellLocal)
        {
            d_nnz[rowI++] += solidBlockSize;
            d_nnz[rowI++] += solidBlockSize;
            d_nnz[rowI++] += solidBlockSize;

            if (!twoD)
            {
                d_nnz[rowI++] += solidBlockSize;
            }
        }
        else
        {
            o_nnz[rowI++] += solidBlockSize;
            o_nnz[rowI++] += solidBlockSize;
            o_nnz[rowI++] += solidBlockSize;

            if (!twoD)
            {
                o_nnz[rowI++] += solidBlockSize;
            }
        }
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

    // Set non-zeros for each interface solid cell
    // In parallel, we need to determine if the coupled fluid cell is on
    // this processor (d_nnz) or another processor (o_nnz)
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatchRef = solidMesh.boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatchRef.faceCells();

    // Build zone-level global fluid cell IDs
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();

    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatchRef = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatchRef.faceCells();

    const globalIndex& fluidGlobalCells =
        refCast<const foamPetscSnesHelper>(fluid()).globalCells();

    // Local field: global cell ID for each local fluid interface face
    scalarField localFluidGlobalCellIDs(fluidPatchRef.size(), 0);
    forAll(fluidPatchRef, faceI)
    {
        localFluidGlobalCellIDs[faceI] =
            fluidGlobalCells.toGlobal(fluidFaceCells[faceI]);
    }

    // Broadcast to zone level (all processors get the full array)
    const scalarField zoneFluidGlobalCellIDs
    (
        fluidGlobalPatch.patchFaceToGlobal(localFluidGlobalCellIDs)
    );

    // Zone-level face map: solid zone face -> fluid zone face
    const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
    const labelList& solidFaceToZone = solidGlobalPatch.faceToGlobalAddr();

    // Loop over LOCAL solid interface faces
    forAll(solidPatchRef, solidFaceI)
    {
        const label solidCellID = solidFaceCells[solidFaceI];

        // Map local solid face -> zone solid face -> zone fluid face
        const label solidZoneFaceI = solidFaceToZone[solidFaceI];
        const label fluidZoneFaceI = solidFaceMap[solidZoneFaceI];
        const label globalFluidCellID =
            label(zoneFluidGlobalCellIDs[fluidZoneFaceI]);

        // Check if the coupled fluid cell is on this processor
        const bool fluidCellLocal =
            fluidGlobalCells.isLocal(globalFluidCellID);

        // Calculate the row index for this cells first scalar equation
        label rowI = solidCellID*solidBlockSize;

        // Accumulate the number of non-zeros (a corner cell may have
        // multiple faces on the interface, each coupling to a different
        // fluid cell)
        // e.g., The x-momentum equation could have a coefficient for the
        // fluid x/y/z velocity and fluid pressure
        if (fluidCellLocal)
        {
            d_nnz[rowI++] += fluidBlockSize;
            d_nnz[rowI++] += fluidBlockSize;

            if (!twoD)
            {
                d_nnz[rowI++] += fluidBlockSize;
            }
        }
        else
        {
            o_nnz[rowI++] += fluidBlockSize;
            o_nnz[rowI++] += fluidBlockSize;

            if (!twoD)
            {
                o_nnz[rowI++] += fluidBlockSize;
            }
        }
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

    // Global cells for the fluid and solid regions
    const globalIndex& fluidGlobalCells =
        refCast<foamPetscSnesHelper>(fluid()).globalCells();
    const globalIndex& solidGlobalCells =
        refCast<foamPetscSnesHelper>(solid()).globalCells();

    // The fluid interface is a prescribed velocity (fixedValue) condition
    // where we assume the fluid wall velocity is equal to the mesh velocity of
    // the adjacent cell centre. This approximation is sufficiently
    // accurate as a preconditioner for the matrix and will not affect the
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

    // No scaling: the Afs Jacobian uses the standard coupling d(U_f)/d(U_s) = 1.
    // (The Liu interface condition modifies only the residual; the Jacobian
    // approximation is kept as the direct-mapping baseline for consistency.)
    const scalar liuScale = 1.0;
    const scalar pressureScale =
        refCast<const fluidModels::newtonIcoFluid>(fluid()).pressureScaleFactor();

    // First we will insert the contribution to the fluid momentum equation
    // coming from the diffusion term

    // Second we will insert the contribution to the fluid continuity
    // (pressure) equation, where the div(U) term should use the adjacent
    // solid cell velocity instead of the known boundary face velocity

    // Build zone-level global solid cell IDs for parallel face mapping
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();

    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const fvPatch& solidPatch = solidMesh().boundary()[solidPatchID];
    const labelList& solidFaceCells = solidPatch.faceCells();

    // Local field: global cell ID for each local solid interface face
    scalarField localSolidGlobalCellIDs(solidPatch.size(), 0);
    forAll(solidPatch, faceI)
    {
        localSolidGlobalCellIDs[faceI] =
            solidGlobalCells.toGlobal(solidFaceCells[faceI]);
    }

    // Broadcast to zone level (all processors get the full array)
    const scalarField zoneSolidGlobalCellIDs
    (
        solidGlobalPatch.patchFaceToGlobal(localSolidGlobalCellIDs)
    );

    // Zone-level face map and local-to-zone addressing
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();
    const labelList& fluidFaceToZone = fluidGlobalPatch.faceToGlobalAddr();

    forAll(fluidPatch, fluidFaceI)
    {
        // Fluid cell adjacent to the interface
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // Map local fluid face -> zone fluid face -> zone solid face
        const label fluidZoneFaceI = fluidFaceToZone[fluidFaceI];
        const label solidZoneFaceI = fluidFaceMap[fluidZoneFaceI];

        // Get global solid cell ID from the zone-level array
        const label globalSolidCellID =
            label(zoneSolidGlobalCellIDs[solidZoneFaceI]);

        // Global block row ID of fluid matrix
        const label globalBlockRowI =
            fluidGlobalCells.toGlobal(fluidCellID);

        // Global block column ID of solid matrix
        const label globalBlockColI = globalSolidCellID;

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*fluidBlockSize;
        label globalColI = globalBlockColI*solidBlockSize;

        // Momentum coefficient for this face (scaled by liuScale for the
        // Liu Eq.31 interface condition)
        PetscScalar value = liuScale*fluidPatchDiffusionCoeffs[fluidFaceI];

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
        // (scaled by liuScale for the Liu Eq.31 interface condition)
        value =
            pressureScale*liuScale
           *fluidPatchPressureCoeffs[fluidFaceI][vector::X];

        globalRowI++; // pressure equation
        globalColI = globalBlockColI*solidBlockSize; // x velocity
        CHKERRQ
        (
            MatSetValues
            (
                Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        value =
            pressureScale*liuScale
           *fluidPatchPressureCoeffs[fluidFaceI][vector::Y];
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
            value =
                pressureScale*liuScale
               *fluidPatchPressureCoeffs[fluidFaceI][vector::Z];
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
    // compact stencil. For pressure, we assume the traction on a fluid
    // interface face is equal to the pressure at the centre of the adjacent
    // fluid cell. For the viscous stress (when passViscousStress is enabled),
    // we use the snGrad approximation: the viscous traction depends on the
    // adjacent fluid cell velocity through the boundary face gradient.
    // This approximation is sufficiently accurate as a preconditioner for
    // the matrix and will not affect the converged solution (which is
    // entirely governed by formResidual)

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

    // Zone-level face map and local-to-zone addressing for parallel
    const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();
    const labelList& solidFaceToZone = solidGlobalPatch.faceToGlobalAddr();

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
    fluidModels::newtonIcoFluid& newtonFluid =
        refCast<fluidModels::newtonIcoFluid>(fluid());

    // Build zone-level global fluid cell IDs for parallel face mapping
    const globalIndex& fluidGlobalCells =
        refCast<foamPetscSnesHelper>(fluid()).globalCells();

    scalarField localFluidGlobalCellIDs(fluidPatch.size(), 0);
    forAll(fluidPatch, faceI)
    {
        localFluidGlobalCellIDs[faceI] =
            fluidGlobalCells.toGlobal(fluidFaceCells[faceI]);
    }

    const scalarField zoneFluidGlobalCellIDs
    (
        fluidGlobalPatch.patchFaceToGlobal(localFluidGlobalCellIDs)
    );

    // Pressure coupling coefficients: solidSf * rho (local to solid, OK)
    const vectorField pressureCoeffs
    (
        solidPatch.Sf()*newtonFluid.rho().value()
    );

    // Viscous coupling coefficients (only when passViscousStress is enabled)
    // The viscous traction from the fluid on the solid is:
    //   t_viscous = -rho * nf & (nuEff * (gradU + gradU^T))
    // Using the snGrad approximation for the boundary face gradient:
    //   d(t_viscous_i)/d(U_cell_k) = rho*nuEff*deltaCoeffs*(delta_ik + nf_i*nf_k)
    // So the force contribution is:
    //   d(F_i)/d(U_cell_k) = rho*nuEff*deltaCoeffs*|Sf|*(delta_ik + nf_i*nf_k)
    // In parallel, the viscous coefficients are on the fluid side so we
    // broadcast them to the zone level for cross-processor access
    scalarField zoneViscDiffCoeffs;
    vectorField zoneFluidNf;
    if (passViscousStress_)
    {
        const scalarField fluidPatchNuEff
        (
            newtonFluid.turbulence().nuEff(fluidPatchID)
        );

        scalarField localViscDiffCoeffs
        (
            newtonFluid.rho().value()
           *fluidPatchNuEff
           *fluidPatch.deltaCoeffs()
           *fluidPatch.magSf()
        );

        vectorField localFluidNf(fluidPatch.nf());

        // Broadcast to zone level for parallel access
        zoneViscDiffCoeffs =
            fluidGlobalPatch.patchFaceToGlobal(localViscDiffCoeffs);
        zoneFluidNf =
            fluidGlobalPatch.patchFaceToGlobal(localFluidNf);
    }

    const globalIndex& solidGlobalCells =
        refCast<foamPetscSnesHelper>(solid()).globalCells();

    forAll(solidPatch, solidFaceI)
    {
        // Map local solid face -> zone solid face -> zone fluid face
        const label solidZoneFaceI = solidFaceToZone[solidFaceI];
        const label fluidZoneFaceI = solidFaceMap[solidZoneFaceI];

        // Solid cell is local
        const label solidCellID = solidFaceCells[solidFaceI];

        // Get global fluid cell ID from the zone-level array
        const label globalFluidCellID =
            label(zoneFluidGlobalCellIDs[fluidZoneFaceI]);

        // Global block row ID of solid matrix
        const label globalBlockRowI =
            solidGlobalCells.toGlobal(solidCellID);

        // Global block column ID of fluid matrix
        const label globalBlockColI = globalFluidCellID;

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // --- Pressure coupling ---
        // The column corresponds to the pressure equation in the fluid cell
        {
            label globalRowI = globalBlockRowI*solidBlockSize;
            const label globalColI =
                globalBlockColI*fluidBlockSize + (fluidBlockSize - 1);

            // Manually insert the 3 scalar coefficients (2 in 2-D)
            PetscScalar value = -pressureCoeffs[solidFaceI][vector::X];
            CHKERRQ
            (
                MatSetValues
                (
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            globalRowI++;
            value = -pressureCoeffs[solidFaceI][vector::Y];
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
                value = -pressureCoeffs[solidFaceI][vector::Z];
                CHKERRQ
                (
                    MatSetValues
                    (
                        Asf, 1, &globalRowI, 1, &globalColI,
                        &value, ADD_VALUES
                    )
                );
            }
        }

        // --- Viscous stress coupling ---
        // The columns correspond to the velocity equations in the fluid cell
        if (passViscousStress_)
        {
            // Use zone-level coefficients for parallel access
            const scalar D = zoneViscDiffCoeffs[fluidZoneFaceI];
            const vector& nf = zoneFluidNf[fluidZoneFaceI];

            // Number of velocity components
            const label nVelCmpts = twoD ? 2 : 3;

            for (label i = 0; i < nVelCmpts; ++i)
            {
                for (label k = 0; k < nVelCmpts; ++k)
                {
                    // d(F_i)/d(U_k) = D * (delta_ik + nf_i*nf_k)
                    PetscScalar value = D*((i == k ? 1.0 : 0.0) + nf[i]*nf[k]);

                    PetscInt globalRowI =
                        globalBlockRowI*solidBlockSize + i;
                    PetscInt globalColI =
                        globalBlockColI*fluidBlockSize + k;

                    CHKERRQ
                    (
                        MatSetValues
                        (
                            Asf, 1, &globalRowI, 1, &globalColI,
                            &value, ADD_VALUES
                        )
                    );
                }
            }
        }
    }

    return 0;
}


void newtonQuasiMonolithicCouplingInterface::mapInterfaceMotionUToFluidU()
{
    // Lookup the fluid interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    fvPatchVectorField& fluidPatchU =
        boundaryFieldRef(fluid().U())[fluidPatchID];
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


void newtonQuasiMonolithicCouplingInterface::setLiuInterfaceVelocity()
{
    // Set the fluid interface velocity using Liu (2014) Eq.31:
    //   U_fluid = (3/4)*U_solid + (1/2)*U_solid_old - (1/4)*U_solid_oldold
    //
    // This is required for energy-stable temporal coupling. The (3/4, 1/2,
    // -1/4) coefficients arise from the compatibility between the BDF2
    // temporal discretisation and the trapezoidal stress averaging in the
    // solid equations. See Liu et al., JCP 270, 2014, Eq. 31 and Section 4.

    // Lookup interface mapping
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

    // Global patches for parallel face data transfer
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();
    const labelList& fluidFaceToZone = fluidGlobalPatch.faceToGlobalAddr();

    // Solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
    const labelList& solidFaceCells =
        solidMesh().boundary()[solidPatchID].faceCells();

    // Access solid velocity at current, old, and old-old time levels
    const vectorField& solidUI = solid().U();
    const vectorField& solidUOld = solid().U().oldTime();
    const vectorField& solidUOldOld = solid().U().oldTime().oldTime();

    // Build local solid interface velocities at each time level using
    // adjacent cell-centre values
    const label nSolidFaces = solidMesh().boundary()[solidPatchID].size();
    vectorField localSolidU(nSolidFaces);
    vectorField localSolidUOld(nSolidFaces);
    vectorField localSolidUOldOld(nSolidFaces);

    forAll(localSolidU, faceI)
    {
        const label cellI = solidFaceCells[faceI];
        localSolidU[faceI] = solidUI[cellI];
        localSolidUOld[faceI] = solidUOld[cellI];
        localSolidUOldOld[faceI] = solidUOldOld[cellI];
    }

    // Broadcast to zone level for parallel access
    const vectorField zoneSolidU
    (
        solidGlobalPatch.patchFaceToGlobal(localSolidU)
    );
    const vectorField zoneSolidUOld
    (
        solidGlobalPatch.patchFaceToGlobal(localSolidUOld)
    );
    const vectorField zoneSolidUOldOld
    (
        solidGlobalPatch.patchFaceToGlobal(localSolidUOldOld)
    );

    // Set the fluid interface velocity with Liu Eq.31 coefficients
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    fvPatchVectorField& fluidPatchU =
        boundaryFieldRef(fluid().U())[fluidPatchID];

    forAll(fluidPatchU, fluidFaceI)
    {
        const label fluidZoneFaceI = fluidFaceToZone[fluidFaceI];
        const label solidZoneFaceI = fluidFaceMap[fluidZoneFaceI];

        fluidPatchU[fluidFaceI] =
            0.75*zoneSolidU[solidZoneFaceI]
          + 0.50*zoneSolidUOld[solidZoneFaceI]
          - 0.25*zoneSolidUOldOld[solidZoneFaceI];
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

    // Zone-level face map: fluidZoneFace -> solidZoneFace
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Global patches for parallel face data transfer
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();
    const labelList& fluidFaceToZone = fluidGlobalPatch.faceToGlobalAddr();

    // Lookup the fluid mesh interface patch
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];

    // Lookup the mesh motion displacement field
    fvPatchVectorField& motionPatchD =
        boundaryFieldRef(motionSolid().D())[fluidPatchID];
    if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
    {
        FatalErrorInFunction
            << "The meshMotionFluid interface patch must be of type "
            << "'fixedValue'" << abort(FatalError);
    }

    // Lookup the mesh motion velocity field
    fvPatchVectorField& motionPatchU =
        boundaryFieldRef(motionSolid().U())[fluidPatchID];

    // Lookup the solid interface patch
    const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];

    // Map the solid interface displacement to the motion interface
    // In parallel, we broadcast solid data to zone level, map solid->fluid
    // zone faces, then extract local fluid values
    if (extrapolateSolidInterfaceDisplacement_)
    {
        const fvPatchVectorField& solidPatchD =
            solid().D().boundaryField()[solidPatchID];
        const fvPatchVectorField& solidPatchU =
            solid().U().boundaryField()[solidPatchID];

        // Broadcast to zone level
        const vectorField zoneSolidD
        (
            solidGlobalPatch.patchFaceToGlobal(solidPatchD)
        );
        const vectorField zoneSolidU
        (
            solidGlobalPatch.patchFaceToGlobal(solidPatchU)
        );

        // Map from solid zone to fluid local via zone face mapping
        forAll(motionPatchD, fluidFaceI)
        {
            const label fluidZoneFaceI = fluidFaceToZone[fluidFaceI];
            const label solidZoneFaceI = fluidFaceMap[fluidZoneFaceI];

            motionPatchD[fluidFaceI] = zoneSolidD[solidZoneFaceI];
            motionPatchU[fluidFaceI] = zoneSolidU[solidZoneFaceI];
        }
    }
    else
    {
        const labelList& solidFaceCells =
            solidMesh().boundary()[solidPatchID].faceCells();
        const vectorField& solidDI = solid().D();
        const vectorField& solidUI = solid().U();

        // Build local solid cell-centre D and U, then broadcast to zone level
        vectorField localSolidCellD(solidMesh().boundary()[solidPatchID].size());
        vectorField localSolidCellU(solidMesh().boundary()[solidPatchID].size());
        forAll(localSolidCellD, faceI)
        {
            localSolidCellD[faceI] = solidDI[solidFaceCells[faceI]];
            localSolidCellU[faceI] = solidUI[solidFaceCells[faceI]];
        }

        const vectorField zoneSolidCellD
        (
            solidGlobalPatch.patchFaceToGlobal(localSolidCellD)
        );
        const vectorField zoneSolidCellU
        (
            solidGlobalPatch.patchFaceToGlobal(localSolidCellU)
        );

        forAll(motionPatchD, fluidFaceI)
        {
            const label fluidZoneFaceI = fluidFaceToZone[fluidFaceI];
            const label solidZoneFaceI = fluidFaceMap[fluidZoneFaceI];

            motionPatchD[fluidFaceI] = zoneSolidCellD[solidZoneFaceI];
            motionPatchU[fluidFaceI] = zoneSolidCellU[solidZoneFaceI];
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
            boundaryFieldRef(motionSolid().pointD())[fluidPatchID]
        )
    )
    {
        FatalErrorInFunction
            << "The meshMotionFluid pointD interface patch must be of type "
            << "'fixedValue'" << abort(FatalError);
    }

    // Set the mesh interface pointD
    // Use "==" to reassign fixedValue
    boundaryFieldRef(motionSolid().pointD())[fluidPatchID] == meshPatchPointD;

    // Correct boundary conditions to enforce the new patch values on the
    // internal field
    motionSolid().pointD().correctBoundaryConditions();
}


void newtonQuasiMonolithicCouplingInterface::restoreOldTimeState
(
    const pointField& oldFluidPoints
)
{
    dynamicFvMesh& fMesh = fluidMesh();
    volVectorField& UFluid = fluid().U();
    volScalarField& p = fluid().p();
    surfaceScalarField& phi = fluid().phi();
    volVectorField& USolid = solid().U();
    volVectorField& D = solid().D();
    volVectorField& DMotion = motionSolid().D();

    UFluid = UFluid.oldTime();
    UFluid.correctBoundaryConditions();

    p = p.oldTime();
    p.correctBoundaryConditions();

    USolid = USolid.oldTime();
    USolid.correctBoundaryConditions();

    D = D.oldTime();
    D.correctBoundaryConditions();

    solid().mechanical().grad(D, solid().gradD());
    solid().mechanical().interpolate(D, solid().gradD(), solid().pointD());
    solid().pointD().correctBoundaryConditions();

    solid().sigma() = solid().sigma().oldTime();
    solid().sigma().correctBoundaryConditions();

    DMotion = DMotion.oldTime();
    DMotion.correctBoundaryConditions();

    motionSolid().pointD() = motionSolid().pointD().oldTime();
    motionSolid().pointD().correctBoundaryConditions();

    motionSolid().sigma() = motionSolid().sigma().oldTime();
    motionSolid().sigma().correctBoundaryConditions();

    fMesh.movePoints(oldFluidPoints);

    phi = fvc::interpolate(UFluid) & fMesh.Sf();

    if (fMesh.changing())
    {
        fvc::makeRelative(phi, UFluid);
    }
}


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
        primitiveFieldRef(UFluid),
        0, // Location of U
        fluidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );
    UFluid.correctBoundaryConditions();

    // Retrieve the fluid pressure
    foamPetscSnesHelper::ExtractFieldComponents
    (
        xx,
        primitiveFieldRef(p),
        fluidBlockSize - 1, // Location of p component
        fluidBlockSize
    );
    p.correctBoundaryConditions();

    // Retrieve the solid velocity
    const label solidFirstEqnID = UFluid.size()*fluidBlockSize;
    foamPetscSnesHelper::ExtractFieldComponents
    (
        &xx[solidFirstEqnID],
        primitiveFieldRef(USolid),
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
        fluid().mesh(), // Used by helper utilities; fvSolution is overridden
        solutionLocation::CELLS, // Used by helper functions
        fsiProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true, // Will PETSc be used
        runTime.system()
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
    liuInterfaceCondition_
    (
        fsiProperties().lookupOrDefault<Switch>
        (
            "liuInterfaceCondition",
            false
        )
    ),
    nRegions_(2),
    subMatsPtr_(nullptr),
    isFluid_(nullptr),
    isFluidVelocity_(nullptr),
    isFluidPressure_(nullptr),
    isSolid_(nullptr),
    configuredNestedFluidSplit_(false),
    nestMat_(nullptr),
    Pmat_(nullptr),
    tsLogPtr_()
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
        << "solidSystemScaleFactor = " << solidSystemScaleFactor_ << nl
        << "liuInterfaceCondition = " << liuInterfaceCondition_ << endl;

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
    Time& time = physicsModel::runTime();
    volVectorField& UFluid = fluid().U();
    surfaceScalarField& phi = fluid().phi();
    volScalarField& p = fluid().p();
    volVectorField& USolid = solid().U();
    volVectorField& D = solid().D();
    volVectorField& DMotion = motionSolid().D();
    const bool twoD = fluid().twoD();
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    const Switch adjustTimeStep
    (
        time.controlDict().lookupOrDefault("adjustTimeStep", false)
    );

    const label maxTimeStepRetries
    (
        fsiProperties().lookupOrDefault<label>("maxTimeStepRetries", 10)
    );

    label timeStepRetry = 0;

    while (true)
    {
        const dimensionedScalar& deltaT = time.deltaT();
        const scalar failedTimeValue = time.value();
        const scalar failedDeltaT = time.deltaTValue();
        const label oldTimeIndex = time.timeIndex() - 1;
        const scalar oldTimeValue = failedTimeValue - failedDeltaT;
        const pointField oldFluidPoints(fluidMesh().points());
        const TimeState retryTimeState =
            retryTimeStateBuilder::rollbackState(time);

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
        //        mesh given D[n] on the fluid-solid interface, where D is known
        //        on the interface from the step (a)

        // First, we will map DSolid to DMotion at the interface
        // This function maps DSolid to DMotion, pointDSolid to pointDMotion,
        // and USolid to UFluid at the interface
        // NOTE: mapping of USolid to UFluid is not needed here, but I already
        // had this function which did this - todo: we can remove that step
        mapInterfaceSolidToMeshMotion();


        //    (c) Update the fluid mesh motion based on the result from (b),
        //        i.e. fluidMesh.update()
        // This updates DMotion based on the interface D which was set in
        // step (b). In addition, this function moves the fluid mesh using this
        // D field
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
                0.5
               *(gSum(sumPhi)/gSum(fluidMesh().V().field()))*deltaT.value();

            Info<< "Fluid mesh Courant number mean: " << meanMeshCoNum
                << " max: " << meshCoNum << endl;
        }


        //    (d) Evaluate the predicted fluid fluid (phi) and fluid mesh fluid
        //        (phiMotion): these are used for the linearised fluid
        //        convection term

        // Extrapolate the fluid velocity
        if
        (
            Switch
            (
                fluid().fluidProperties().lookup
                (
                    "fluidFluxExtrapolationAlgorithm1"
                )
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
                0.5
               *(gSum(sumPhi)/gSum(fluidMesh().V().field()))*deltaT.value();

            Info<< "Fluid Courant number mean: " << meanCoNum
                << " max: " << CoNum << endl;
        }

        // In Jamain and Joshi, they update the motionU (fluid mesh velocity)
        // based on the current, old and old-old fluid mesh positions (or
        // equivalently the displacements) using a 2nd order method. These mesh
        // velocities are used in the linearised convection term.
        // In OpenFOAM, mesh.phi() will already return a 2nd order
        // approximation of the fluid mesh fluid (velocity times area at the
        // faces), which is also consistent with the geometric conservation law,
        // assuming we use a 2nd order time discretisation.
        // So we don't need to do anything here; we just need to use mesh.phi()


        //    (e) Compute the element-level stabilisation parameters : this step
        //        is not required in our finite volume implementation. Instead,
        //        we will incorporate typical finite volume stabilisation terms
        //        (e.g., Rhie-Chow) when forming the monolithic system in
        //        step (f)
        //
        // Nothing to do here.


        //    (f) Solve the coupled monolithic system for UFluid[n], pFluid[n],
        //        and USolid[n] using the fluid mesh from step (c) and the
        //        linear linearised mesh flux from step (d)

        // Keep the pre-solve state so a failed SNES solve can be retried.
        foamPetscSnesHelper::storeSolutionBackup();

        // We will call PETSc SNES to form and solve the coupled monolithic
        // system. Note that this will call formResidual, formJacobian,
        // initialiseResidual, initialiseJacobian as defined above.
        // This coupled system will be linear so the PETSc options should be
        // selected with this in mind.
        Info<< "Solving the monolithic momentum-continuity system for Up using "
            << "PETSc SNES" << endl;
        const int solveStatus = foamPetscSnesHelper::solve(true);

        if (solveStatus >= 0)
        {
            break;
        }

        VecCopy
        (
            foamPetscSnesHelper::solutionBackup(),
            foamPetscSnesHelper::solution()
        );

        restoreOldTimeState(oldFluidPoints);

        if (!adjustTimeStep)
        {
            FatalErrorInFunction
                << "PETSc SNES failed to converge and the previous "
                << "fluid-solid time-step state has been restored, but "
                << "`adjustTimeStep` is disabled." << nl
                << "Enable `adjustTimeStep` to retry the failed time step "
                << "with a reduced deltaT."
                << abort(FatalError);
        }

        ++timeStepRetry;

        if (timeStepRetry > maxTimeStepRetries)
        {
            FatalErrorInFunction
                << "Exceeded the maximum number of failed PETSc retries ("
                << maxTimeStepRetries << ") at time " << failedTimeValue
                << " with deltaT = " << failedDeltaT << nl
                << "Set a larger `maxTimeStepRetries` if you would like more "
                << "recovery attempts."
                << abort(FatalError);
        }

        static_cast<TimeState&>(time) = retryTimeState;
        time.setTime(oldTimeValue, oldTimeIndex);
        setDeltaT(time);

        if (time.deltaTValue() >= failedDeltaT*(1.0 - SMALL))
        {
            FatalErrorInFunction
                << "PETSc SNES failed to converge at the minimum allowed time "
                << "step. The old-time state has been restored, but deltaT "
                << "could not be reduced below " << failedDeltaT
                << " for a retry."
                << abort(FatalError);
        }

        ++time;

        Info<< "Retrying the failed PETSc time step with deltaT = "
            << time.deltaTValue() << " at Time = "
            << time.timeName() << nl << endl;
    }

    // Retrieve the solution from PETSc and copy into UFluid, p, and USolid
    retrieveSolution
    (
        foamPetscSnesHelper::solution(),
        UFluid, p, USolid,
        fluidBlockSize, solidBlockSize,
        twoD
    );

    return true;
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

    // If a split preconditioner is used, let it know where the fluid unknowns
    // and solid unknowns are located. Guard with configuredNestedFluidSplit_ to
    // avoid duplicate IS registration (formJacobian also sets IS).
    if
    (
        pct
     && !std::strcmp(pct, PCFIELDSPLIT)
     && !configuredNestedFluidSplit_
    )
    {
        Info<< "    Defining the fluid-solid field indices for the fluid-solid"
            << nl << "    split preconditioner" << endl;
        PCFieldSplitSetIS(pc, "fluid", isFluid());
        PCFieldSplitSetIS(pc, "solid", isSolid());
        configuredNestedFluidSplit_ = true;
    }
}


label newtonQuasiMonolithicCouplingInterface::formResidual
(
    Vec f,
    const Vec x
)
{
    const PetscBool reportBlockResiduals =
        petscOptionEnabled("-s4f_petsc_report_block_residuals");

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // Preliminaries
    volVectorField& UFluid = fluid().U();
    volScalarField& p = fluid().p();
    volVectorField& USolid = solid().U();
    // const label fluidSystemEnd = fluidMesh().nCells()*fluidBlockSize;
    // const label solidFirstEqnID = fluidSystemEnd;


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

        if (liuInterfaceCondition_)
        {
            setLiuInterfaceVelocity();
        }
        else
        {
            mapInterfaceMotionUToFluidU();
        }
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

    // Zone-level face map: fluidZoneFace -> solidZoneFace
    // zoneBToZoneAFaceMap: size = zoneA (fluid), map[fluidZoneFaceI] = solidZoneFaceI
    const labelList& fluidFaceMap = interfaceMap.zoneBToZoneAFaceMap();

    // Global patch for parallel face data transfer
    const globalPolyPatch& fluidGlobalPatch = interfaceMap.globalPatchA();
    const globalPolyPatch& solidGlobalPatch = interfaceMap.globalPatchB();

    // Calculate the solid interface traction
    // Flip the sign as the solid normals point in the opposite direction to
    // the fluid normals
    const label solidPatchID =
        fluidSolidInterface::solidPatchIndices()[0];

    // Broadcast local fluid traction to zone level
    const vectorField negFluidTraction(-fluidTraction);
    const vectorField zoneFluidTraction
    (
        fluidGlobalPatch.patchFaceToGlobal(negFluidTraction)
    );

    // Map from fluid zone faces to solid zone faces at zone level
    vectorField zoneSolidTraction(solidGlobalPatch.globalPatch().size(), vector::zero);
    forAll(fluidFaceMap, fluidZoneFaceI)
    {
        zoneSolidTraction[fluidFaceMap[fluidZoneFaceI]] =
            zoneFluidTraction[fluidZoneFaceI];
    }

    // Extract local solid traction from zone level
    const vectorField solidTraction
    (
        solidGlobalPatch.globalFaceToPatch(zoneSolidTraction)
    );

    // Lookup the displacement interface traction patch and set the traction
    fvPatchVectorField& solidPatchD =
        boundaryFieldRef(solid().D())[solidPatchID];
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
    if
    (
        !isA<solidModels::nonLinGeomTotalLagVelocitySolid>(solid())
     && !isA<solidModels::linGeomVelocitySolid>(solid())
    )
    {
        FatalErrorInFunction
            << "Currently, the solid model must be of type "
            << solidModels::nonLinGeomTotalLagVelocitySolid::typeName
            << " or "
            << solidModels::linGeomVelocitySolid::typeName
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


    // 4. Map the solid interface velocity to the fluid interface
    //    When liuInterfaceCondition is enabled, the fluid interface velocity
    //    uses the temporally consistent condition from Liu (2014) Eq.31:
    //      U_fluid = (3/4)*U_solid + (1/2)*U_solid_old - (1/4)*U_solid_oldold
    //    This is O(dt^2) accurate and ensures energy stability at the FSI
    //    interface for the BDF2/trapezoidal temporal discretisation.
    mapInterfaceSolidToMeshMotion();

    if (liuInterfaceCondition_)
    {
        setLiuInterfaceVelocity();
    }
    else
    {
        mapInterfaceMotionUToFluidU();
    }


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
        Vec fFluid = nullptr;
        VecGetSubVector(f, isFluid(), &fFluid);
        VecScale(fFluid, fluidSystemScaleFactor_);
        VecRestoreSubVector(f, isFluid(), &fFluid);

        Vec fSolid = nullptr;
        VecGetSubVector(f, isSolid(), &fSolid);
        VecScale(fSolid, solidSystemScaleFactor_);
        VecRestoreSubVector(f, isSolid(), &fSolid);
    }

    // 7. Report per-block residual norms (diagnostic)
    if (reportBlockResiduals)
    {
        PetscReal fluidNormScaled = 0;
        PetscReal solidNormScaled = 0;

        Vec fFluid = nullptr;
        VecGetSubVector(f, isFluid(), &fFluid);
        VecNorm(fFluid, NORM_2, &fluidNormScaled);
        VecRestoreSubVector(f, isFluid(), &fFluid);

        Vec fSolid = nullptr;
        VecGetSubVector(f, isSolid(), &fSolid);
        VecNorm(fSolid, NORM_2, &solidNormScaled);
        VecRestoreSubVector(f, isSolid(), &fSolid);

        const PetscReal fluidNorm =
            fluidNormScaled/mag(fluidSystemScaleFactor_);
        const PetscReal solidNorm =
            solidNormScaled/mag(solidSystemScaleFactor_);

        Info<< "    Per-block residual norms:"
            << " fluid(scaled)=" << fluidNormScaled
            << " solid(scaled)=" << solidNormScaled
            << " fluid(unscaled)=" << fluidNorm
            << " solid(unscaled)=" << solidNorm
            << endl;
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

    // Determine the assembly matrix: always the MatNest.
    // On the first call, jac IS the MatNest and we store a pointer to it.
    // On subsequent calls (when Pmat_ is active), jac may be the converted
    // MPIAIJ, so we use the stored nestMat_ for assembly.
    Mat assemblyMat = jac;
    if (nestMat_ != nullptr)
    {
        // Use the stored MatNest for assembly
        assemblyMat = nestMat_;
        MatZeroEntries(assemblyMat);
    }
    else
    {
        // First call: store the MatNest pointer
        nestMat_ = jac;
    }

    // Get access to the sub-matrices from the MatNest
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(assemblyMat, &nr, &nc, &subMats));

    if (nr != 2 || nc != 2)
    {
        FatalErrorInFunction
            << "The matrix has the wrong number of sub matrices: "
            << "nr = " << nr << ", nc = " << nc << abort(FatalError);
    }

    // Zero the entries of the MatNest (may already be zeroed by callback)
    MatZeroEntries(assemblyMat);

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
            subMats[0][0], xFluid, true // allow explicit flux
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

    // Explicitly assemble each sub-matrix before assembling the nest.
    // This ensures off-process entries (from parallel coupling blocks)
    // are properly flushed.
    CHKERRQ(MatAssemblyBegin(subMats[0][0], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(subMats[0][0], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyBegin(subMats[0][1], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(subMats[0][1], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyBegin(subMats[1][0], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(subMats[1][0], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyBegin(subMats[1][1], MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(subMats[1][1], MAT_FINAL_ASSEMBLY));

    // Complete nest matrix assembly
    CHKERRQ(MatAssemblyBegin(assemblyMat, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(assemblyMat, MAT_FINAL_ASSEMBLY));

    // Scale the sub-matrices
    MatScale(subMats[0][0], fluidSystemScaleFactor_);
    MatScale(subMats[0][1], fluidSystemScaleFactor_);
    MatScale(subMats[1][0], solidSystemScaleFactor_);
    MatScale(subMats[1][1], solidSystemScaleFactor_);

    // Convert the assembled MatNest to monolithic AIJ for preconditioning.
    // This allows standard PCs (LU, ILU, ASM, bjacobi, etc.) to work with the
    // block-structured system. Required in both serial and parallel because
    // MatNest does not support operations like MatIncreaseOverlap (ASM).
    const bool firstMatConvert = (Pmat_ == nullptr);
    {
        if (Pmat_ == nullptr)
        {
            // First call: create the MPIAIJ by conversion
            CHKERRQ
            (
                MatConvert(assemblyMat, MATAIJ, MAT_INITIAL_MATRIX, &Pmat_)
            );

            // Update the SNES to use Pmat_ for preconditioning
            Mat currentJac;
            SNESJacobianFn *jacFunc;
            void *jacCtx;
            CHKERRQ
            (
                SNESGetJacobian
                (
                    snes(), &currentJac, nullptr, &jacFunc, &jacCtx
                )
            );
            CHKERRQ
            (
                SNESSetJacobian(snes(), currentJac, Pmat_, jacFunc, jacCtx)
            );
        }
        else
        {
            // Subsequent calls: reuse the existing AIJ structure
            CHKERRQ
            (
                MatConvert(assemblyMat, MATAIJ, MAT_REUSE_MATRIX, &Pmat_)
            );
        }
    }

    // Configure field split index sets when the top-level PC is fieldsplit.
    // Two-level split:
    //   Level 1: fluid vs solid (top-level PC)
    //   Level 2: velocity vs pressure within fluid (nested PC, optional)
    if (!configuredNestedFluidSplit_)
    {
        KSP ksp;
        PC pc;
        const char* pct = nullptr;

        SNESGetKSP(snes(), &ksp);
        KSPGetPC(ksp, &pc);
        PCGetType(pc, &pct);

        // Only proceed if the top-level PC is a field split
        if (pct && std::strcmp(pct, PCFIELDSPLIT) == 0)
        {
            configuredNestedFluidSplit_ = true;

            // Set top-level IS for fluid/solid split
            Info<< "    Setting fluid/solid field split index sets" << endl;
            PCFieldSplitSetIS(pc, "fluid", isFluid());
            PCFieldSplitSetIS(pc, "solid", isSolid());

            // Ensure the KSP has the correct preconditioning matrix
            // before setup. The SNES may not have updated the KSP
            // operators yet (this is inside formJacobian, before
            // SNES processes the returned matrices).
            {
                Mat amat;
                SNESGetJacobian(snes(), &amat, nullptr, nullptr, nullptr);
                Mat pmat = (Pstream::parRun() && Pmat_ != nullptr)
                    ? Pmat_ : assemblyMat;
                KSPSetOperators(ksp, amat, pmat);
            }

            // Note: nested velocity/pressure split within the fluid
            // sub-block is configured via PETSc command-line options
            // using block_size (e.g. -fieldsplit_fluid_pc_fieldsplit_
            // block_size 3 for 2D, 4 for 3D).
        }
    }

    // One-time check: in parallel, a direct solver (LU/Cholesky) on the
    // solid sub-block crashes with SIGBUS because the fieldsplit sub-matrix
    // is distributed.  Use "redundant" to gather the small solid block to
    // one rank instead.
    if (firstMatConvert && Pstream::parRun())
    {
        KSP ksp;
        PC pc;
        const char* pct = nullptr;

        SNESGetKSP(snes(), &ksp);
        KSPGetPC(ksp, &pc);
        PCGetType(pc, &pct);

        if (pct && std::strcmp(pct, PCFIELDSPLIT) == 0)
        {
            KSPSetOperators
            (
                ksp,
                assemblyMat,
                Pmat_ ? Pmat_ : assemblyMat
            );
            PCSetUp(pc);

            PetscInt nSub = 0;
            KSP* subKSPs = nullptr;
            PCFieldSplitGetSubKSP(pc, &nSub, &subKSPs);

            // Solid is the last split (index nSub - 1)
            if (nSub >= 2)
            {
                PC solidPC;
                KSPGetPC(subKSPs[nSub - 1], &solidPC);

                const char* solidPCType = nullptr;
                PCGetType(solidPC, &solidPCType);

                if
                (
                    solidPCType
                 && (
                        !std::strcmp(solidPCType, PCLU)
                     || !std::strcmp(solidPCType, PCCHOLESKY)
                    )
                )
                {
                    FatalErrorInFunction
                        << "The solid sub-block uses direct solver '"
                        << solidPCType
                        << "' in parallel, which causes bus errors." << nl
                        << "Use 'redundant' to gather the solid block"
                        << " to one rank:" << nl << nl
                        << "    fieldsplit_solid_pc_type redundant;" << nl
                        << "    fieldsplit_solid_redundant_pc_type lu;"
                        << nl << nl
                        << exit(FatalError);
                }
            }

            PetscFree(subKSPs);
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

#endif // OPENFOAM_COM

#endif // ifdef USE_PETSC

// ************************************************************************* //
