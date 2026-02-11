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

#include "newtonMonolithicCouplingInterface.H"
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

defineTypeNameAndDebug(newtonMonolithicCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, newtonMonolithicCouplingInterface, dictionary
);


// * * * * * * * * * * * Private Member Functions* * * * * * * * * * * * * * //


void newtonMonolithicCouplingInterface::makeIsFluidIsMotionIsSolid() const
{
    if (isFluid_ != nullptr || isMotion_ != nullptr || isSolid_ != nullptr)
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

        // Fluid and solid block size
        const label fluidBlockSize = twoD ? 3 : 4;

        // Set the number local scalar equations
        const label n = blockn*fluidBlockSize;

        // Set the global system size
        NFluid = returnReduce(n, sumOp<label>());

        // Create the index sets, where it is assumed the motion equations come first
        ISCreateStride(PETSC_COMM_WORLD, NFluid, 0, 1, &isFluid_);
    }

    // Fluid mesh motion
    label NMotion = -1;
    {
        // Set the number local block equations
        // Note: the motion mesh is the fluid mesh
        const label blockn = fluid().mesh().nCells();

        // Motion and solid block size
        const label motionBlockSize = twoD ? 2 : 3;

        // Set the number local scalar equations
        const label n = blockn*motionBlockSize;

        // Set the global system size
        NMotion = returnReduce(n, sumOp<label>());

        // Create the index sets, where it is assumed the motion equations come first
        ISCreateStride(PETSC_COMM_WORLD, NMotion, NFluid, 1, &isMotion_);
    }

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
        ISCreateStride(PETSC_COMM_WORLD, NSolid, NMotion, 1, &isSolid_);
    }
}


foamPetscSnesHelper& newtonMonolithicCouplingInterface::motion()
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


solidModel& newtonMonolithicCouplingInterface::motionSolid()
{
    return refCast<solidModel>(motion());
}


void newtonMonolithicCouplingInterface::createSubMatsAndMat
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


label newtonMonolithicCouplingInterface::initialiseAfm
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
    std::vector<label> d_nnz(scalarRowN, fluidBlockSize);
    std::vector<label> o_nnz(scalarRowN, 0);

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


label newtonMonolithicCouplingInterface::initialiseAms
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


label newtonMonolithicCouplingInterface::initialiseAsf
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


label newtonMonolithicCouplingInterface::formAfm
(
    Mat Afm,
    const label fluidBlockSize,
    const label motionBlockSize,
    const bool twoD
)
{
    // We will include two forms of mesh-in-fluid coupling:
    //   i.   The interface velocity equal the mesh motion interface velocity
    //   ii.  The mesh flux in the advection term of the fluid momentum equation

    // TODO: CHECK! I am adding pressure equation coeffs; do I need to?

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


    // Coupling i. The interface velocity equal the mesh motion interface
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

    // Lookup the fluid interface patch
    // The mesh motion is the same patch as it is the same mesh
    const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
    const fvPatch& fluidPatch = fluidMesh().boundary()[fluidPatchID];
    const labelList& fluidFaceCells = fluidPatch.faceCells();

    // Lookup the fluid patch
    const fvPatchVectorField& fluidPatchU =
        fluid().U().boundaryField()[fluidPatchID];

    // Check the fluid U patch type is fixedValue (not just derived from it)
    // if (!isA<fixedValueFvPatchVectorField>(fluidPatchU))
    if (fluidPatchU.type() != "fixedValue")
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
    // NOTE: we should match the solid's ddt scheme, not the mesh motion ddt
    // scheme
    const word motionDdtScheme =
        // word(fluidMesh().ddtScheme("ddt(" + motionSolid().D().name() +')'));
        word(solidMesh().ddtScheme("ddt(" + solid().D().name() +')'));

    // The known fluid boundary face value is now replaced by the adjacent
    // motion cell velocity
    scalarField fluidPatchDiffusionCoeffs(fluidPatch.size(), 0.0);
    vectorField fluidPatchPressureCoeffs(fluidPatch.size(), vector::zero);
    if (motionDdtScheme == "Euler")
    {
        // Diffusion coefficient is nu*|Sf|/(|n & d|*dt), where dt comes
        // from converting the motion displacement to velocity
        const scalar deltaT = motionSolid().time().deltaTValue();
        fluidPatchDiffusionCoeffs =
            fluidPatchNuEff*fluidPatch.magSf()*fluidPatch.deltaCoeffs()/deltaT;

        // Pressure coefficient from div(U)
        fluidPatchPressureCoeffs = -fluidPatchSf/deltaT;
    }
    else if (motionDdtScheme == "backward")
    {
        // Diffusion coefficient is 3.0*nu*|Sf|/(|n & d|*2.0*dt), where (3/2)*dt
        // comes from converting the motion displacement to velocity
        const scalar deltaT = motionSolid().time().deltaTValue();
        fluidPatchDiffusionCoeffs =
            3.0*fluidPatchNuEff*fluidPatch.magSf()*fluidPatch.deltaCoeffs()
           /(2.0*deltaT);

        // Pressure coefficient from div(U)
        fluidPatchPressureCoeffs = -3.0*fluidPatchSf/(2.0*deltaT);
    }
    else
    {
        // steady-state does not really make sense as this would mean the fluid
        // velocity at the interface would be zero... or maybe it is OK.. to
        // check!
        FatalErrorInFunction
            << "Unknown motion ddtScheme " << motionDdtScheme << ": only "
            << "Euler and backward currently allowed"
            << exit(FatalError);
    }

    // First we will insert the contribution to the fluid momentum equation
    // coming from the diffusion term

    // Second we will insert the contribution to the fluid continuity
    // (pressure) equation, where the div(U) term should use the adjacent
    // solid cell velocity instead of the known boundary face velocity

    forAll(fluidPatch, fluidFaceI)
    {
        // Fluid cell adjacent to the interface
        const label fluidCellID = fluidFaceCells[fluidFaceI];

        // To help readability
        const label motionCellID = fluidCellID;

        // We will add coefficients at block row "fluidCellID" and block
        // column "motionCellID"
        // Note that the nested monolithic matrix takes care of offsetting
        // the rows/columns when the submatrices are inserted into the
        // nested matrix

        // Global block row ID of fluid matrix
        const label globalBlockRowI = globalCells.toGlobal(fluidCellID);

        // Global block column ID of motion matrix
        const label globalBlockColI = globalCells.toGlobal(motionCellID);

        // CAREFUL: we cannot use MatSetValuesBlocked as it only works with
        // square block coefficients, so we will insert the scalar
        // coefficients manually

        // Calculate the scalar global row ID (not the block row ID)
        label globalRowI = globalBlockRowI*fluidBlockSize;
        label globalColI = globalBlockColI*motionBlockSize;

        // Momentum coefficient for this face
        PetscScalar value = fluidPatchDiffusionCoeffs[fluidFaceI];

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        CHKERRQ
        (
            MatSetValues
            (
                Afm, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        globalRowI++;
        globalColI++;
        CHKERRQ
        (
            MatSetValues
            (
                Afm, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
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
                    Afm, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }

        // Secondly we will insert the contributions for the pressure
        // equation

        // Manually insert the 3 scalar coefficients (2 in 2-D)
        value = fluidPatchPressureCoeffs[fluidFaceI][vector::X];

        globalRowI++; // pressure equation
        globalColI = globalBlockColI*motionBlockSize; // x motion displacement
        CHKERRQ
        (
            MatSetValues
            (
                Afm, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        value = fluidPatchPressureCoeffs[fluidFaceI][vector::Y];
        globalColI++; // z motion displacement
        CHKERRQ
        (
            MatSetValues
            (
                Afm, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
            )
        );

        if (!twoD)
        {
            value = fluidPatchPressureCoeffs[fluidFaceI][vector::Z];
            globalColI++; // z motion displacement
            CHKERRQ
            (
                MatSetValues
                (
                    Afm, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );
        }
    }


    // Coupling ii. The mesh flux in the advection term of the fluid momentum
    // equation

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


label newtonMonolithicCouplingInterface::formAms
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
    // Diffusion coefficient is impK*|Sf|/(n & d)
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


label newtonMonolithicCouplingInterface::formAsf
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


void newtonMonolithicCouplingInterface::mapInterfaceMotionUToFluidU()
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


void newtonMonolithicCouplingInterface::mapInterfaceSolidToMeshMotion()
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
    motionSolid().pointD().boundaryFieldRef()[fluidPatchID] ==
        meshPatchPointD;

    // Correct boundary conditions to enforce the new patch values on the
    // internal field
    motionSolid().pointD().correctBoundaryConditions();
}


void newtonMonolithicCouplingInterface::predict()
{
    if (runTime().timeIndex() == 0)
    {
        return;
    }

    Info<< "Linear predictor" << endl;

    // Predict solution using previous time steps

    // Velocity
    // fluid().U() =
    //     fluid().U().oldTime() + fluid().A()*runTime().deltaT();

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

    // Mesh motion
    motionSolid().D() =
        motionSolid().D().oldTime() + motionSolid().U()*runTime().deltaT();

    // Update the mesh motion interface
    if (solidToMeshCoupling_)
    {
        mapInterfaceSolidToMeshMotion();
    }

    // Insert the OpenFOAM fields into the PETSc solution vector

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // The scalar row at which the solid equations start
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

    // The scalar row at which the motion equations start
    const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

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

    // Insert the motion displacement
    foamPetscSnesHelper::InsertFieldComponents
    (
        motionSolid().D().primitiveFieldRef(),
        &xx[motionFirstEqnID],
        0, // Location of first component
        solidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );

    // Restore the solution vector
    VecRestoreArray(foamPetscSnesHelper::solution(), &xx);
}


void newtonMonolithicCouplingInterface::resetFieldsToOldTime()
{
    Info<< "Resetting primary fields to their old time values" << nl << endl;

    // Velocity
    fluid().U() = fluid().U().oldTime();

    // Pressure
    fluid().p() = fluid().p().oldTime();

    // Displacement
    solid().D() = solid().D().oldTime();

    // Mesh motion
    motionSolid().D() = motionSolid().D().oldTime();

    // Update the interface fields
    mapInterfaceSolidToMeshMotion();
    mapInterfaceMotionUToFluidU();

    // Insert the OpenFOAM fields into the PETSc solution vector

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // The scalar row at which the solid equations start
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

    // The scalar row at which the motion equations start
    const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

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

    // Insert the motion displacement
    foamPetscSnesHelper::InsertFieldComponents
    (
        motionSolid().D().primitiveFieldRef(),
        &xx[motionFirstEqnID],
        0, // Location of first component
        solidBlockSize,
        twoD ? labelList({0,1}) : labelList({0,1,2})
    );

    // Restore the solution vector
    VecRestoreArray(foamPetscSnesHelper::solution(), &xx);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonMonolithicCouplingInterface::newtonMonolithicCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "UpD",
        fileName
        (
            fsiProperties().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        fluid().mesh(), // fluid/fvSolution will be used
        solutionLocation::NONE,
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
    nRegions_(3),
    subMatsPtr_(nullptr),
    isFluid_(nullptr),
    isSolid_(nullptr),
    tsLogPtr_(),
    oldTimeValue_(runTime.value()),
    nConsecutiveFailedSolves_(0),
    maxAllowedConsecutiveFailedSolves_
    (
        fsiProperties().lookupOrDefault<label>
        (
            "maxAllowedConsecutiveFailedSolves", 5
        )
    )
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

    // Store old time values
    fluid().U().storeOldTime();
    fluid().p().storeOldTime();
    solid().D().storeOldTime();
    motionSolid().D().storeOldTime();

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
        << "meshToFluidCoupling: " << meshToFluidCoupling_ << nl
        << "solidToMeshCoupling: " << solidToMeshCoupling_ << nl
        << "extrapolateSolidInterfaceDisplacement: "
        << extrapolateSolidInterfaceDisplacement_ << nl
        << "passViscousStress: " << passViscousStress_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void newtonMonolithicCouplingInterface::setDeltaT(Time& runTime)
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


bool newtonMonolithicCouplingInterface::evolve()
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

    // Ensure boundary conditions are up-to-date
    fluid().U().correctBoundaryConditions();
    fluid().p().correctBoundaryConditions();
    solid().D().correctBoundaryConditions();
    motionSolid().D().correctBoundaryConditions();

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

    // Mesh motion block size
    const label motionBlockSize = solidBlockSize;

    // The scalar row at which the motion equations start
    const label motionFirstEqnID = fluidMesh().nCells()*fluidBlockSize;

    // The scalar row at which the solid equations start
    const label solidFirstEqnID =
        motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;

    // Store the SNES solution
    foamPetscSnesHelper::storeSolutionBackup();

    // Solve the nonlinear system and check the convergence
    // If unsuccesful, reduce the time, reset the fields and try again
    bool failed = false;
    do
    {
        // Solve the monolithic system
        Info<< "Calling snes::solve" << endl;
        failed = false;
        foamPetscSnesHelper::solve(true);

        SNESConvergedReason reason;
        SNESGetConvergedReason(foamPetscSnesHelper::snes(), &reason);

        if (reason < 0)
        {
            Info<< "Solution fail: reducing the time step and trying again"
                << endl;

            failed = true;

            // Reduce the time-step
            setDeltaT(const_cast<Time&>(runTime()));

            // Reset the time
            const_cast<Time&>(runTime()).setTime
            (
                oldTimeValue_ + runTime().deltaTValue(), runTime().timeIndex()
            );

            // Reset the fields to the old time fields
            resetFieldsToOldTime();

            // Reset the SNES solution
            VecCopy
            (
                foamPetscSnesHelper::solutionBackup(),
                foamPetscSnesHelper::solution()
            );

            // Reset the SNES object
            foamPetscSnesHelper::resetSnes();

            // Clear rAUf field
            refCast<fluidModels::newtonIcoFluid>(fluid()).clearRAUf();
        }

        // Keep track of the number of failed solves and quit if the maximum
        // allowed is reached
        if (nConsecutiveFailedSolves_++ == maxAllowedConsecutiveFailedSolves_)
        {
            FatalErrorInFunction
                << nl
                << "The monolithic nonlinear solver (SNES) failed to converge "
                << "after "
                << maxAllowedConsecutiveFailedSolves_ << " consecutive attempts."
                << nl << "This may indicate:\n"
                << "  - A time step that is too large for stability," << nl
                << "  - Highly nonlinear or poorly scaled physics," << nl
                << "  - Insufficient or inappropriate preconditioning." << nl
                << "Suggestions:" << nl
                << "  - Try reducing the initial time step in 'controlDict'," << nl
                << "  - Check solver settings or residuals," << nl
                << "  - Enable SNES monitoring for diagnostics."
                << exit(FatalError);
        }
    }
    while (failed);

    // Reset the nConsecutiveFailedSolves counter
    nConsecutiveFailedSolves_ = 0;

    // Update the old time value
    oldTimeValue_ = runTime().value();

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
    solid().mechanical().interpolate
    (
        solid().D(), solid().gradD(), solid().pointD()
    );
    solid().pointD().correctBoundaryConditions();
    motionSolid().mechanical().interpolate
    (
        motionSolid().D(), motionSolid().gradD(), motionSolid().pointD()
    );
    motionSolid().pointD().correctBoundaryConditions();

    // Increment of displacement
    solid().DD() = solid().D() - solid().D().oldTime();

    // Increment of point displacement
    solid().pointDD() = solid().pointD() - solid().pointD().oldTime();

    // Velocity
    solid().U() = fvc::ddt(solid().D());
    motionSolid().U() = fvc::ddt(motionSolid().D());

    // Update the mesh motion interface
    if (solidToMeshCoupling_)
    {
        mapInterfaceSolidToMeshMotion();
    }

    // Update the fluid interface velocity
    if (meshToFluidCoupling_)
    {
        mapInterfaceMotionUToFluidU();
    }

    // Update phi
    // Is there a neater way to do this? Maybe ask the fluid model...
    {
        const volTensorField Fm(I + motionSolid().gradD().T());
        const volScalarField Jm(det(Fm));
        const volTensorField invFm(inv(Fm));
        const surfaceTensorField Fmf(fvc::interpolate(Fm));
        const surfaceScalarField Jmf(det(Fmf));
        const surfaceTensorField invFmf(inv(Fmf));
        const surfaceVectorField deformedSf
        (
            Jmf*invFmf.T() & fluidMesh().Sf()
        );
        fluid().phi() =
            deformedSf
          & (fvc::interpolate(fluid().U() - motionSolid().U()));
    }

    // Update fluid velocity acceleration
    fluid().A() = fvc::ddt(fluid().U());

    // Update fluid rate of change of pressure
    fluid().dpdt() = fvc::ddt(fluid().p());

    return 0;
}


label newtonMonolithicCouplingInterface::initialiseJacobian(Mat& jac)
{
    // A fluid-solid interaction problem with a moving mesh (arbitrary
    // Lagrangian Eulerian) fluid can be expressed in Ax=b form, where the
    // matrix A is
    //
    //     *-----------------*
    //     | Aff | Afm |  0  |
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
    //  - Afm: mesh motion terms in the fluid equations, including the
    //         interface motion
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations
    //
    // Note
    //  - Afs is empty as the fluid equation do not directly depend on the solid
    //    displacement; instead, the solid interface motion indirectly comes from
    //    the fluid mesh interface motion
    //  - Amf is empty as the the mesh motion does not directly
    //    depend on the fluid solution
    //  - Asm is empty as the solid equations do not depend on the fluid
    //    mesh motion

    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

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

    // Create the empty submatrices and skip Afs, Amf and Asm, which will be set
    // to null ptrs
    labelPairHashSet nullSubMats;
    nullSubMats.insert(labelPair(0, 2)); // Afs
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
    foamPetscSnesHelper::initialiseJacobian
    (
        subMats[0][0], fluid().mesh(), fluidBlockSize, false
    );
    PetscObjectSetName((PetscObject)subMats[0][0], "Aff");

    // Amm
    foamPetscSnesHelper::initialiseJacobian
    (
        subMats[1][1], fluid().mesh(), motionBlockSize, false
    );
    PetscObjectSetName((PetscObject)subMats[1][1], "Amm");

    // Ass
    foamPetscSnesHelper::initialiseJacobian
    (
        subMats[2][2], solid().mesh(), solidBlockSize, false
    );
    PetscObjectSetName((PetscObject)subMats[2][2], "Ass");

    // Initialise the four off diagonal submatrices:
    //  - Afm: mesh motion terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Afm
    initialiseAfm
    (
        subMats[0][1], fluidMesh(), fluidBlockSize, motionBlockSize, twoD
    );
    PetscObjectSetName((PetscObject)subMats[0][1], "Afm");

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

    return 0;
}


label newtonMonolithicCouplingInterface::initialiseSolution(Vec& x)
{
    // Set twoD flag
    const bool twoD = fluid().twoD();

    // Set fluid and solid block sizes
    const label fluidBlockSize = twoD ? 3 : 4;
    const label solidBlockSize = twoD ? 2 : 3;

    // Mesh motion block size
    const label motionBlockSize = solidBlockSize;

    // Count the number of unknowns
    label n = 0;
    label N = 0;

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


label newtonMonolithicCouplingInterface::formResidual
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

    // Considerations on the order of updating the fluid, solid, and mesh motion
    // residuals
    //  - solid depends on the fluid interface traction
    //  - mesh motion depends on the solid interface displacement
    //  - fluid depends on the entire mesh motion flux field
    //
    // Approach
    // 1. Update the fluid velocity and pressure fields and calculate the
    //    traction at the fluid interface
    // 2. Map the fluid interface traction to the solid interface
    // 3. Update the solid residual, which now has the correct interface
    //    traction
    // 4. Map the solid interface displacement to the mesh motion interface
    // 5. Update the mesh motion residual, which now has the correct interface
    //    displacement
    // 6. Map the mesh motion interface velocity to the fluid interface
    // 7. Update the fluid residual, which now has the correct interface
    //    velocity and mesh motion deformation fields (Fm, Jm)
    // 8. Apply scaling factor to solid and motion equations to preserve the
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
        Vec xFluid = nullptr;
        VecGetSubVector(x, isFluid(), &xFluid);
        VecSetBlockSize(xFluid, fluidBlockSize);
        foamPetscSnesHelper::ExtractFieldComponents<vector>
        (
            xFluid,
            U.primitiveFieldRef(),
            0, // Location of U
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
            xFluid,
            p.primitiveFieldRef(),
            fluidBlockSize - 1 // Location of p component
        );

        p.correctBoundaryConditions();

        VecRestoreSubVector(x, isFluid(), &xFluid);
    }

    // Extract D field (needed for motion interface)
    {
        Vec xSolid = nullptr;
        VecGetSubVector(x, isSolid(), &xSolid);
        VecSetBlockSize(xSolid, solidBlockSize);
        foamPetscSnesHelper::ExtractFieldComponents<vector>
        (
            xSolid,
            solid().D().primitiveFieldRef(),
            0, // Location of first component
            twoD ? labelList({0,1}) : labelList({0,1,2})
        );

        // Enforce the boundary conditions
        solid().D().correctBoundaryConditions();

        VecRestoreSubVector(x, isSolid(), &xSolid);
    }

    // Update gradient of displacement
    solid().mechanical().grad
    (
        solid().D(), solid().gradD()
    );

    // Map D to motionD
    mapInterfaceSolidToMeshMotion();


    // Check for unphysical values for J as an indicator for divergence
    {
        const volScalarField J(det(I + solid().gradD().T()));
        const scalar maxJ = max(J).value();
        const scalar minJ = min(J).value();

        if (minJ < 1e-8 || maxJ > 1e6)
        {
            Info<< "The solution is diverging as minJ = " << minJ
                << " and maxJ = " << maxJ
                << endl;

            return -1;
        }
    }

    // 2. Map the fluid interface traction to the solid interface
    if (fluidToSolidCoupling_ && coupled())
    {
        // Fluid interface traction
        const label fluidPatchID =
            fluidSolidInterface::fluidPatchIndices()[0];

        {
            Vec xMotion = nullptr;
            VecGetSubVector(x, isMotion(), &xMotion);
            VecSetBlockSize(xMotion, motionBlockSize);
            foamPetscSnesHelper::ExtractFieldComponents<vector>
            (
                xMotion,
                motionSolid().D().primitiveFieldRef(),
                0, // Location of first component
                twoD ? labelList({0,1}) : labelList({0,1,2})
            );

            // Enforce the boundary conditions
            motionSolid().D().correctBoundaryConditions();

            VecRestoreSubVector(x, isMotion(), &xMotion);
        }

        // Map D to motionD
        mapInterfaceSolidToMeshMotion();

        // Update gradient of displacement
        motionSolid().mechanical().grad
        (
            motionSolid().D(), motionSolid().gradD()
        );

        // Map D to motionD
        mapInterfaceSolidToMeshMotion();

        const tensorField& gradD =
            motionSolid().gradD().boundaryField()[fluidPatchID];
        const tensorField Fm(I + gradD.T());
        const scalarField Jm(det(Fm));
        const tensorField invFm(inv(Fm));
        const vectorField& Sf = fluidMesh().boundary()[fluidPatchID].Sf();
        const vectorField deformedSf(Jm*invFm.T() & Sf);
        // const tensorField gradU
        // (
        //     fvc::grad(fluid().U())().boundaryField()[patchID]
        // ); // do I need to calculate the whole field?
        // tvF.ref() = rho_.value()*deformedSf & (nuEff*invFm.T() & gradU);
        const vectorField fluidNf
        (
            fluidMesh().boundary()[fluidPatchID].nf()
        );
        const vectorField fluidDeformedNf(deformedSf/mag(deformedSf));
        vectorField fluidTraction
        (
          - fluid().patchPressureForce(fluidPatchID)*fluidDeformedNf
          // - fluid().patchPressureForce(fluidPatchID)*fluidNf
        );
        if (passViscousStress_)
        {
            const scalar rho = 1;
            const scalar nuEff = 1;
            tensorField gradU
            (
                fvc::grad(fluid().U())().boundaryField()[fluidPatchID]
            );
            gradU +=
                fluidNf*fluid().U().boundaryField()[fluidPatchID].snGrad()
              - ((I - sqr(fluidNf)) & gradU);
            fluidTraction +=
                // refCast<fluidModels::newtonIcoFluid>
                // (
                //     fluid()
                // ).patchViscousForce
                // (
                //     fluidPatchID,
                //     motionSolid()
                // );
                rho*nuEff*fluidDeformedNf
              & (
                  //gradU // not rotated!
                  //invFm & gradU
                  Jm*invFm.T() & gradU // why does this work?
                  //invFm.T() & gradU // should be right
                  //gradU & invFm.T()
                  //gradU & invFm
                );
                //rho*fluidNf & (nuEff*invFm.T() & gradU);
                //rho*nuEff*(fluidDeformedNf & (gradU & invFm));
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

    // Check for unphysical values for Jm as an indicator for divergence
    {
        const volScalarField J(det(I + motionSolid().gradD().T()));
        const scalar maxJ = max(J).value();
        const scalar minJ = min(J).value();

        if (minJ < 1e-8 || maxJ > 1e6)
        {
            Info<< "The solution is diverging as minJm = " << minJ
                << " and maxJm = " << maxJ
                << endl;

            return -1;
        }
    }

    // 3. Update the solid residual, which now has the correct interface
    //    traction
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

    // 4. Map the solid interface displacement and velocity to the mesh
    //    motion interface
    if (solidToMeshCoupling_ && coupled())
    {
        mapInterfaceSolidToMeshMotion();
    }

    // 5. Update the mesh motion residual, which now has the correct interface
    //    displacement
    {
        Vec xMotion = nullptr;
        Vec fMotion = nullptr;
        VecGetSubVector(x, isMotion(), &xMotion);
        VecGetSubVector(f, isMotion(), &fMotion);
        VecSetBlockSize(xMotion, motionBlockSize);
        VecSetBlockSize(fMotion, motionBlockSize);
        refCast<foamPetscSnesHelper>(motionSolid()).formResidual
        (
            fMotion, xMotion
        );
        VecRestoreSubVector(x, isMotion(), &xMotion);
        VecRestoreSubVector(f, isMotion(), &fMotion);
    }

    // 6. Map the mesh motion interface velocity to the fluid interface
    // Note: it is assumed the motion.U() is dDm/dt even for steady-state cases
    if (meshToFluidCoupling_ && coupled())
    {
        mapInterfaceMotionUToFluidU();
    }

    // 7. Update the fluid residual, which now has the correct interface
    //    velocity and mesh motion
    // Note that the fluid equations are first in the f (residual) and x
    // (solution) lists
    //refCast<foamPetscSnesHelper>(fluid()).formResidual(f, x);
    // The fluid residual should be calculated over the reference
    // configuration
    {
        Vec xFluid = nullptr;
        Vec fFluid = nullptr;
        VecGetSubVector(x, isFluid(), &xFluid);
        VecGetSubVector(f, isFluid(), &fFluid);
        VecSetBlockSize(xFluid, fluidBlockSize);
        VecSetBlockSize(fFluid, fluidBlockSize);
        refCast<fluidModels::newtonIcoFluid>(fluid()).formResidual
        (
            fFluid, xFluid, motionSolid()
        );
        VecRestoreSubVector(x, isFluid(), &xFluid);
        VecRestoreSubVector(f, isFluid(), &fFluid);
    }

    // 8. Apply scaling factor to solid and motion equations to preserve the
    // condition number of the monolithic system
    // TODO: what about scaling the fluid equations?

    PetscScalar* ff;
    VecGetArray(f, &ff);

    const label solidSystemEnd =
        solidFirstEqnID + solidMesh().nCells()*solidBlockSize;
    for (int i = solidFirstEqnID; i < solidSystemEnd; ++i)
    {
        ff[i] *= solidSystemScaleFactor_;
    }

    const label motionSystemEnd =
        motionFirstEqnID + fluidMesh().nCells()*motionBlockSize;
    for (int i = motionFirstEqnID; i < motionSystemEnd; ++i)
    {
        ff[i] *= motionSystemScaleFactor_;
    }

    VecRestoreArray(f, &ff);

    PetscFunctionReturn(0);
}


label newtonMonolithicCouplingInterface::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Our monolithic system matrix will take the form:
    //
    //     *-----------------*
    //     | Aff | Afm |  0  |
    //     |-----------------|
    //     |  0  | Amm | Ams |
    //     |-----------------|
    //     | Asf |  0  | Ass |
    //     *-----------------*
    //

    // We will assembly the six sub-matrices:

    // Diagonal submatrices:
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations
    //
    // Off-diagonal submatrices:
    //  - Afm: mesh motion terms in the fluid equations: mesh flux in advection
    //    term and interface motion
    //  - Ams: solid terms (interface motion) in the mesh motion equation
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

    if (nr != 3 || nc != 3)
    {
        FatalErrorInFunction
            << "The matrix has the wrong number of sub matrices: "
            << "nr = " << nr << ", nc = " << nc << abort(FatalError);
    }

    // Zero the entries
    MatZeroEntries(jac);

    // Set the motion block size
    const label motionBlockSize = solidBlockSize;

    // Form diagonal submatrices
    //  - Aff: fluid equations (momentum and pressure)
    //  - Amm: mesh motion equations in the fluid domain
    //  - Ass: solid equations

    // Aff
    {
        Vec xFluid = nullptr;
        VecGetSubVector(x, isFluid(), &xFluid);
        VecSetBlockSize(xFluid, fluidBlockSize);
        refCast<foamPetscSnesHelper>(fluid()).formJacobian
        (
            subMats[0][0], xFluid
        );
        VecRestoreSubVector(x, isFluid(), &xFluid);
    }

    // Amm
    {
        Vec xMotion = nullptr;
        VecGetSubVector(x, isMotion(), &xMotion);
        VecSetBlockSize(xMotion, motionBlockSize);
        refCast<foamPetscSnesHelper>(motionSolid()).formJacobian
        (
            subMats[1][1], xMotion
        );
        VecRestoreSubVector(x, isMotion(), &xMotion);
    }

    // Ass
    {
        Vec xSolid = nullptr;
        VecGetSubVector(x, isSolid(), &xSolid);
        VecSetBlockSize(xSolid, solidBlockSize);
        refCast<foamPetscSnesHelper>(solid()).formJacobian
        (
            subMats[2][2], xSolid
        );
        VecRestoreSubVector(x, isSolid(), &xSolid);
    }

    // Form off-diagonal submatrices:
    //  - Afm: mesh motion terms in the fluid equations
    //  - Ams: solid terms (interface motion) in the mesh motion equation
    //  - Asf: fluid terms (pressure on interface) in solid equations

    // Afm
    if (meshToFluidCoupling_) // && coupled())
    {
        formAfm(subMats[0][1], fluidBlockSize, motionBlockSize, twoD);
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
        MatZeroEntries(subMats[1][2]);
        MatZeroEntries(subMats[2][0]);
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

#endif // OPENFOAM_NOT_EXTEND

#endif // ifdef USE_PETSC

// ************************************************************************* //
