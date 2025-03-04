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
        fluid().mesh().nCells() + solid().mesh().nCells(),
        fsiProperties().lookupOrDefault<Switch>("stopOnPetscError", true),
        true // Will PETSc be used
    ),
    meshMotionPtr_(),
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
    nRegions_(2),
    subMatsPtr_(nullptr)
{
    // Create the fluid mesh motion solver
    meshMotionPtr_ = solidModel::New(runTime, "meshMotionFluid");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool newtonCouplingInterface::evolve()
{
    // Steps
    // 1. Optional: predict solution field using old-time fields
    // 2. Map foam solution fields to PETSc
    // 3. Solve PETSc system
    // 4. Map PETSc solution to foam fields
    // 5. Update an secondary foam fields

    // Take references
    volVectorField& U = fluid().U();
    volScalarField& p = fluid().p();
    surfaceScalarField& phi = fluid().phi();
    volVectorField& D = solid().D();

    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 3;

    // WarningInFunction
    //     << "Todo: Add linear predictor" << endl;

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Access the raw solution data
    const PetscScalar *xx;
    VecGetArrayRead(foamPetscSnesHelper::solution(), &xx);

    // Retrieve the solution
    // Map the PETSc solution to the U field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        //foamPetscSnesHelper::solution(),
        xx,
        U.primitiveFieldRef(),
        0, // Location of U
        fluidBlockSize,
        fluid().twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    U.correctBoundaryConditions();

    // Update the flux
    // CHECK
    // Make relative?
    phi = fvc::interpolate(U) & fluidMesh().Sf();

    // Map the PETSc solution to the p field
    // p is located in the final component
    foamPetscSnesHelper::ExtractFieldComponents<scalar>
    (
        //foamPetscSnesHelper::solution(),
        xx,
        p.primitiveFieldRef(),
        fluidBlockSize - 1, // Location of p component
        fluidBlockSize
    );

    p.correctBoundaryConditions();

    // Extract the displacement, which starts at row fluid.nCells*fluidBlockSize
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;
    ExtractFieldComponents
    (
        &xx[solidFirstEqnID],
        D.primitiveFieldRef(),
        0, // Location of first component
        solidBlockSize,
        solid().twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Restore the x vector
    VecRestoreArrayRead(foamPetscSnesHelper::solution(), &xx);

    // Interpolate cell displacements to vertices
    solid().mechanical().interpolate(solid().D(), solid().gradD(), solid().pointD());

    // Increment of displacement
    solid().DD() = solid().D() - solid().D().oldTime();

    // Increment of point displacement
    solid().pointDD() = solid().pointD() - solid().pointD().oldTime();

    // Velocity
    solid().U() = fvc::ddt(solid().D());




    // Test mesh motion solver

    // Map solid interface motion to the mesh motion interface
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

        // Lookup the fluid mesh interface patch
        const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];

        // Lookup the mesh motion displacement field
        fvPatchVectorField& motionPatchD =
            meshMotionPtr_->D().boundaryFieldRef()[fluidPatchID];
        if (!isA<fixedValueFvPatchVectorField>(motionPatchD))
        {
            FatalErrorInFunction
                << "The meshMotionFluid interface patch must be of type "
                << "'fixedValue'" << abort(FatalError);
        }

        // Lookup the solid interface patch
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        const fvPatchVectorField& solidPatchD =
            solid().D().boundaryField()[solidPatchID];

        // Map the solid interface displacement to the motion interface
        forAll(motionPatchD, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Extrapolated patch value (larger stencil)
            motionPatchD[fluidFaceI] = solidPatchD[solidFaceID];
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
        meshMotionPtr_->pointD().boundaryFieldRef()[fluidPatchID] ==
            meshPatchPointD;
    }

    // Solve mesh motion
    Info<< nl << "Solving mesh motion" << endl;
    meshMotionPtr_->evolve();

    return 0;
}


label newtonCouplingInterface::initialiseJacobian(Mat& jac)
{
    // Create a nest matrix.

    // We will initialise the four sub-matrices:
    //  - Afluid: fluid
    //  - Asolid: solid
    //  - Asf: solid-in-fluid-coupling
    //  - Afs: fluid-in-solid coupling

    if (nRegions_ != 2)
    {
        FatalErrorInFunction
            << "nRegions must be 2!" << abort(FatalError);
    }

    // Set fluid and solid block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 3;

    // For brevity and convenience, we will store the size and blockSize of the
    // regions in a labelPairList
    const labelPairList nBlocksAndBlockSize
    (
        {
            {fluidMesh().nCells(), fluidBlockSize},
            {solidMesh().nCells(), solidBlockSize}
        }
    );

    // Create arrays (vectors) of IS objects for rows and columns.
    // These will indicate where in the matrix the different regions are
    // located
    std::vector<IS> isRow(nRegions_), isCol(nRegions_);
    label globalOffset = 0;
    for (label r = 0; r < nRegions_; ++r)
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
    // For example, subMat[i*nRegions_ + i] might be the matrix for
    // region i, while subMat[i*nRegions_ + j] (i != j) are the coupling
    // matrices.
    subMatsPtr_ = new Mat[nRegions_*nRegions_];
    for (label i = 0; i < nRegions_; ++i)
    {
        for (label j = 0; j < nRegions_; ++j)
        {
            // Create the submatrix for regions i and j.
            Mat subA;
            MatCreate(PETSC_COMM_WORLD, &subA);

            // Number of rows
            const label nRowBlocks = nBlocksAndBlockSize[i].first();
            const label rowBlockSize = nBlocksAndBlockSize[i].second();
            const label nRowsLocal = nRowBlocks*rowBlockSize;
            const label nRowsGlobal =
                returnReduce(nRowsLocal, sumOp<label>());

            // Number of columns
            const label nColBlocks = nBlocksAndBlockSize[j].first();
            const label colBlockSize = nBlocksAndBlockSize[j].second();
            const label nColsLocal = nColBlocks*colBlockSize;
            const label nColsGlobal =
                returnReduce(nColsLocal, sumOp<label>());

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

            // Store subA in row-major order at index [i*nRegions_ + j]
            subMatsPtr_[i*nRegions_ + j] = subA;
        }
    }

    // Create the nest matrix using the index sets and submatrices.
    jac = Mat();
    MatCreateNest
    (
        PETSC_COMM_WORLD,
        nRegions_,
        isRow.data(),
        nRegions_,
        isCol.data(),
        (Mat*)subMatsPtr_,
        &jac
    );

    // Cleanup IS objects.
    for (label r = 0; r < nRegions_; ++r)
    {
        ISDestroy(&isRow[r]);
        ISDestroy(&isCol[r]);
    }


    // Get access to the sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    // Initialise Afluid based on compact stencil fvMesh
    {
        if (debug)
        {
            InfoInFunction
                << "Initialising Afluid" << endl;
        }

        Mat Afluid = subMats[0][0];
        Foam::initialiseJacobian(Afluid, fluid().mesh(), fluidBlockSize, false);
    }

    // Initialise Asolid based on compact stencil fvMesh
    {
        if (debug)
        {
            InfoInFunction
                << "Initialising Asolid" << endl;
        }

        Mat Asolid = subMats[1][1];
        Foam::initialiseJacobian(Asolid, solid().mesh(), solidBlockSize, false);
    }

    // Initially we assume a conformal FSI interface, where each fluid cell
    // shares a face with a solid cell. So we assume the number of blocks in
    // the Asf (and Afs) matrix is equal to the number of cells at the
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

    // Initialise solid-in-fluid-coupling matrix
    {
        if (debug)
        {
            InfoInFunction
                << "Initialising Asf" << endl;
        }

        Mat Asf = subMats[0][1];

        // CAREFUL: we are setting non-zeros here based on the scalar rows, not
        // the block rows

        // Set matrix type to AIJ (since BAIJ does not support non-square
        // blocks)
        CHKERRQ(MatSetType(Asf, MATMPIAIJ));

        // Total number of scalar rows in the fluid region
        const label scalarRowN = fluidMesh().nCells()*fluidBlockSize;

        // Allocate per-scalar-row nonzeros, initialised to 0
        std::vector<int> d_nnz(scalarRowN, 0);
        std::vector<int> o_nnz(scalarRowN, 0);

        // Set non-zeros for each interface fluid cells
        const labelList& fluidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
        const label fluidPatchID = fluidSolidInterface::fluidPatchIndices()[0];
        const labelList& fluidFaceCells =
            fluidMesh().boundary()[fluidPatchID].faceCells();
        forAll(fluidFaceMap, fluidFaceI)
        {
            const label fluidCellID = fluidFaceCells[fluidFaceI];

            // Calculate the row index for this cells first scalar equation
            label rowI = fluidCellID*fluidBlockSize;

            // Set the number of non-zeros to be the number of solid equations
            // e.g., The x-momentum equation could have a coefficient for the
            // solid x/y/z displacement
            d_nnz[rowI++] = solidBlockSize;
            d_nnz[rowI++] = solidBlockSize;
            d_nnz[rowI++] = solidBlockSize;

            if (fluidBlockSize > 3)
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
            MatMPIAIJSetPreallocation(Asf, 0, d_nnz.data(), 0, o_nnz.data())
        );

        // Raise an error if mallocs are required during matrix assembly
        CHKERRQ(MatSetOption(Asf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    }


    {
        if (debug)
        {
            InfoInFunction
                << "Initialising Afs" << endl;
        }

        Mat Afs = subMats[1][0];

        // CAREFUL: we are setting non-zeros her based on the scalar rows, not
        // the block rows

        // Set matrix type to AIJ (since BAIJ does not support non-square
        // blocks)
        CHKERRQ(MatSetType(Afs, MATMPIAIJ));

        // Total number of scalar rows in the solid region
        const label scalarRowN = solidMesh().nCells()*solidBlockSize;

        // Allocate per-scalar-row nonzeros, initialised to 0
        std::vector<int> d_nnz(scalarRowN, 0);
        std::vector<int> o_nnz(scalarRowN, 0);

        // Set non-zeros for each interface solid cells
        const labelList& solidFaceMap = interfaceMap.zoneAToZoneBFaceMap();
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        const labelList& solidFaceCells =
            solidMesh().boundary()[solidPatchID].faceCells();
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

            if (solidBlockSize > 2)
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
            MatMPIAIJSetPreallocation(Afs, 0, d_nnz.data(), 0, o_nnz.data())
        );

        // Raise an error if mallocs are required during matrix assembly
        CHKERRQ(MatSetOption(Afs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    }

    if (debug)
    {
        InfoInFunction
            << "End" << endl;
    }

    return 0;
}


label newtonCouplingInterface::initialiseSolution(Vec& x)
{
    // Set fluid and solid block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 3;

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
    label n = 0;
    forAll(nBlocksAndBlockSize, regionI)
    {
        const label nBlocksRegionI = nBlocksAndBlockSize[regionI].first();
        const label blockSizeRegionI = nBlocksAndBlockSize[regionI].second();
        n += nBlocksRegionI*blockSizeRegionI;
    }

    // Global system size: total number of scalar equation across all
    // processors
    const label N = returnReduce(n, sumOp<label>());

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
    // Approach
    // 1. Update the fluid velocity and pressure fields and calculate the
    //     traction at the fluid interface
    // 2. Map the fluid interface traction to the solid interface
    // 3. Update the solid residual, which now has the correct interface
    //    traction
    // 4. Map the solid interface velocity to the fluid interface
    // 5. Update the fluid residual, which now has the correct interface
    //    velocity

    // Block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 3;

    // Currently limited to one interface: it should be straight-forward to add
    // a loop over multiple interface => todo
    if (fluidSolidInterface::fluidPatchIndices().size() != 1)
    {
        FatalError
            << "Only one interface patch is currently allowed"
            << abort(FatalError);
    }

    // 1. Update the fluid velocity and pressure fields and calculate the
    //     traction at the fluid interface
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

        if (fluidToSolidCoupling_)
        {
            // Fluid interface traction
            const label fluidPatchID =
                fluidSolidInterface::fluidPatchIndices()[0];
            const vectorField fluidNf
            (
                fluidMesh().boundary()[fluidPatchID].nf()
            );
            const vectorField fluidTraction
            (
                fluid().patchViscousForce(fluidPatchID)
              - fluid().patchPressureForce(fluidPatchID)*fluidNf
            );

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
    }


    // Calculate the solid residual first
    // The first solid scalar equation is at row fluidMesh.nCells*fluidBlockSize
    const label solidFirstEqnID = fluidMesh().nCells()*fluidBlockSize;
    refCast<foamPetscSnesHelper>(solid()).formResidual
    (
        &f[solidFirstEqnID], &x[solidFirstEqnID]
    );

    {
        // Scale the solid residual to preserve the condition number of the
        // monolithic system
        const label solidSystemEnd =
            solidFirstEqnID + solidMesh().nCells()*solidBlockSize;
        for (int i = solidFirstEqnID; i < solidSystemEnd; ++i)
        {
            f[i] *= solidSystemScaleFactor_;
        }
    }

    // We update the velocity at the fluid interface with the solid interface
    // velocity
    // Note: it is assumed the solid.U is dD/dt even for steady-state cases
    if (solidToFluidCoupling_)
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

        // Lookup the solid interface patch
        const label solidPatchID = fluidSolidInterface::solidPatchIndices()[0];
        const fvPatchVectorField& solidPatchU =
            solid().U().boundaryField()[solidPatchID];

        // Map the solid interface velocity to the fluid interface
        // const labelList& solidFaceCells =
        //     solidMesh().boundary()[solidPatchID].faceCells();
        // const vectorField& solidUI = solid().U();
        forAll(fluidPatchU, fluidFaceI)
        {
            const label solidFaceID = fluidFaceMap[fluidFaceI];

            // Extrapolated patch value (larger stencil)
            fluidPatchU[fluidFaceI] = solidPatchU[solidFaceID];

            // Adjacent cell value
            // const label solidCellID = solidFaceCells[solidFaceID];
            // fluidPatchU[fluidFaceI] = solidUI[solidCellID];
        }

        fluid().phi() = fvc::interpolate(fluid().U()) & fluidMesh().Sf();
    }

    // Then, calculate the fluid residual
    // Note that the fluid equations are first in the f (residual) and x
    // (solution) lists
    refCast<foamPetscSnesHelper>(fluid()).formResidual(f, x);

    return 0;
}


label newtonCouplingInterface::formJacobian
(
    Mat jac,
    const PetscScalar *x
)
{
    // We will assembly the four sub-matrices:
    //  - Afluid: fluid
    //  - Asolid: solid
    //  - Asf: solid-in-fluid-coupling
    //  - Afs: fluid-in-solid coupling

    // Zero entries
    MatZeroEntries(jac);

    // Get access to the two sub-matrices
    PetscInt nr, nc;
    Mat **subMats;
    CHKERRQ(MatNestGetSubMats(jac, &nr, &nc, &subMats));

    if (nr != 2 || nc != 2)
    {
        FatalErrorInFunction
            << "The matrix has the wrong number of sub matrices: "
            << "nr = " << nr << ", nc = " << nc << abort(FatalError);
    }

    // Assemble the fluid matrix
    {
        if (debug)
        {
            Info<< "Forming Afluid" << endl;
        }

        // Mat Afluid;
        // MatNestGetSubMat(jac, 0, 0, Afluid);
        Mat Afluid = subMats[0][0];

        if (!isA<foamPetscSnesHelper>(fluid()))
        {
            FatalErrorInFunction
                << "A fluid model derived from foamPetscSnesHelper must be used!"
                << abort(FatalError);
        }

        // The FSI interface is a prescribed velocity condition
        refCast<foamPetscSnesHelper>(fluid()).formJacobian(Afluid, x);
    }

    // Assemble the solid matrix
    {
        if (debug)
        {
            Info<< "Forming Asolid" << endl;
        }

        Mat Asolid = subMats[1][1];

        if (!isA<foamPetscSnesHelper>(solid()))
        {
            FatalErrorInFunction
                << "A solid model derived from foamPetscSnesHelper must be used!"
                << abort(FatalError);
        }

        // The FSI interface is a prescribed velocity condition
        refCast<foamPetscSnesHelper>(solid()).formJacobian(Asolid, x);
    }

    // Set fluid and solid block sizes
    const label fluidBlockSize = fluid().twoD() ? 3 : 4;
    const label solidBlockSize = solid().twoD() ? 2 : 3;

    // const label nFluidCells = fluidMesh().nCells();

    // Assemble the solid-in-fluid coupling matrix
    if (solidToFluidCoupling_)
    {
        if (debug)
        {
            Info<< "Forming Asf" << endl;
        }

        Mat Asf = subMats[0][1];

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
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            globalRowI++;
            globalColI++;
            CHKERRQ
            (
                MatSetValues
                (
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            if (fluidBlockSize > 3)
            {
                globalRowI++;
                globalColI++;
                CHKERRQ
                (
                    MatSetValues
                    (
                        Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
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
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            value = -fluidPatchSf[fluidFaceI][vector::Y]/deltaT;
            globalColI++; // y solid displacement
            CHKERRQ
            (
                MatSetValues
                (
                    Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            if (fluidBlockSize > 3)
            {
                value = -fluidPatchSf[fluidFaceI][vector::Z]/deltaT;
                globalColI++; // y solid displacement
                CHKERRQ
                (
                    MatSetValues
                    (
                        Asf, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                    )
                );
            }
        }
    }

    // Assemble the fluid-in-solid coupling matrix
    if (fluidToSolidCoupling_)
    {
        if (debug)
        {
            Info<< "Forming Afs" << endl;
        }

        Mat Afs = subMats[1][0];

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
                    Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            globalRowI++;
            //globalColI++;
            value = -patchCoeffs[solidFaceI][vector::Y];
            CHKERRQ
            (
                MatSetValues
                (
                    Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                )
            );

            if (solidBlockSize > 2)
            {
                globalRowI++;
                //globalColI++;
                value = -patchCoeffs[solidFaceI][vector::Z];
                CHKERRQ
                (
                    MatSetValues
                    (
                        Afs, 1, &globalRowI, 1, &globalColI, &value, ADD_VALUES
                    )
                );
            }
        }
    }

    {
        // Scale the solid matrix to preserve the condition number of the
        // monolithic system

        // We must assembly the matrix before we can use MatScale
        // Complete matrix assembly
        CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

        MatScale(subMats[1][1], solidSystemScaleFactor_);
        MatScale(subMats[1][0], solidSystemScaleFactor_);
    }

    if (debug)
    {
        Info<< "End of newtonCouplingInterface::formJacobian" << endl;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
