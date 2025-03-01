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

#include "foamPetscSnesHelper.H"
#include "SparseMatrixTemplate.H"
#include "processorFvPatch.H"
#include "symmetryFvPatchFields.H"
#include "symmetryPlaneFvPatchFields.H"
#include "leastSquaresVectors.H"
#include "fvm.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

PetscErrorCode formResidualFoamPetscSnesHelper
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Vec f,        // residual
    void *ctx     // user context
)
{
    const PetscScalar *xx;
    PetscScalar       *ff;
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;

    // Access x and f data
    CHKERRQ(VecGetArrayRead(x, &xx));
    CHKERRQ(VecGetArray(f, &ff));

    // Compute the residual
    if (user->solMod_.formResidual(ff, xx) != 0)
    {
        Foam::FatalError
            << "formResidual(ff, xx) returned an error code!"
            << Foam::abort(Foam::FatalError);
    }

    // Restore the solution and residual vectors
    CHKERRQ(VecRestoreArrayRead(x, &xx));
    CHKERRQ(VecRestoreArray(f, &ff));

    return 0;
}


PetscErrorCode formJacobianFoamPetscSnesHelper
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Mat jac,      // Jacobian
    Mat B,        // Preconditioner matrix (can be jac)
    void *ctx     // user context
)
{
    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Access the OpenFOAM data
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;

    // Check if the matrix is a nested matrix
    PetscBool isNest;
    PetscObjectTypeCompare((PetscObject)B, MATNEST, &isNest);

    // Lookup the matrix info
    MatInfo info;
    if (isNest)
    {
        // For a nested matrix, we will only check the top left sub-matrix - if
        // it is initialised we will assume all sub-matrices are initialised
        Mat sub;
        MatNestGetSubMat(B, 0, 0, &sub);
        MatGetInfo(sub, MAT_LOCAL, &info);
    }
    else
    {
        MatGetInfo(B, MAT_LOCAL, &info);
    }

    // Initialise the matrix if it has yet to be allocated; otherwise zero all
    // entries
    if (info.nz_used)
    {
        // Zero the matrix but do not reallocate the space
        // The "-snes_lag_jacobian -2" PETSc option can be used to avoid
        // re-building the matrix
        // For a nested matrix, this will zero all sub-matrices
        CHKERRQ(MatZeroEntries(B));
    }
    else
    {
        // Set the non-zero structure of the matrix B
        Foam::Info<< "    Initialising the matrix..." << Foam::flush;
        if (user->solMod_.initialiseJacobian(B) != 0)
        {
            Foam::FatalError
                << "initialiseJacobian(B) returned an error code!"
                << Foam::abort(Foam::FatalError);
        }
        Foam::Info<< "done" << Foam::endl;
    }

    // Populate the Jacobian => implemented by the solid model
    if (user->solMod_.formJacobian(B, xx) != 0)
    {
        Foam::FatalError
            << "formJacobian(B, xx) returned an error code!"
            << Foam::abort(Foam::FatalError);
    }

    // Restore solution vector
    CHKERRQ(VecRestoreArrayRead(x, &xx));

    // Complete matrix assembly
    CHKERRQ(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));

    if (jac != B)
    {
        CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
    }

    return 0;
}



Foam::label Foam::initialiseJacobian
(
    Mat jac,
    const fvMesh& mesh,
    const label blockSize
)
{
    // Set the block size
    CHKERRQ(MatSetBlockSize(jac, blockSize));

    // Number of vector unknowns on this processor
    //const int blockn = user->solMod_.globalCells().localSize();
    const int blockn = mesh.nCells();

    // Count the number of non-zeros in the matrix
    // Note: we assume a compact stencil, i.e. face only face neighbours

    // Number of on-processor non-zeros per row
    int* d_nnz = (int*)malloc(blockn*sizeof(int));

    // Number of off-processor non-zeros per row
    int* o_nnz = (int*)malloc(blockn*sizeof(int));

    // Initialise d_nnz and o_nnz to zero
    for (int i = 0; i < blockn; ++i)
    {
        d_nnz[i] = 1; // count diagonal cell
        o_nnz[i] = 0;
    }

    // Take a reference to the mesh
    //const Foam::fvMesh& mesh = user->solMod_.fmesh();

    // Count neighbours sharing an internal face
    const Foam::labelUList& own = mesh.owner();
    const Foam::labelUList& nei = mesh.neighbour();
    forAll(own, faceI)
    {
        const Foam::label ownCellID = own[faceI];
        const Foam::label neiCellID = nei[faceI];
        d_nnz[ownCellID]++;
        d_nnz[neiCellID]++;
    }

    // Count off-processor neighbour cells
    forAll(mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].type() == "processor")
        {
            const Foam::labelUList& faceCells =
                mesh.boundary()[patchI].faceCells();

            forAll(faceCells, fcI)
            {
                const Foam::label cellID = faceCells[fcI];
                o_nnz[cellID]++;
            }
        }
        else if (mesh.boundary()[patchI].coupled())
        {
            // Other coupled boundaries are not implemented
            Foam::FatalError
                << "Coupled boundary are not implemented, except for"
                << " processor boundaries" << Foam::abort(Foam::FatalError);
        }
    }

    // Allocate parallel matrix
    //CHKERRQ(MatMPIAIJSetPreallocation(jac, 0, d_nnz, 0, o_nnz));
    // Allocate parallel matrix with the same conservative stencil per node
    //CHKERRQ(MatMPIAIJSetPreallocation(jac, d_nz, NULL, 0, NULL));
    CHKERRQ(MatMPIBAIJSetPreallocation(jac, blockSize, 0, d_nnz, 0, o_nnz));

    // Raise an error if mallocs are required during matrix assembly
    MatSetOption(jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    // Free memory
    free(d_nnz);
    free(o_nnz);

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(foamPetscSnesHelper, 0);


// * * * * * * * * * * * * * * * Private Function  * * * * * * * * * * * * * //

void foamPetscSnesHelper::makeNeiProcFields(const fvMesh& mesh) const
{
    if
    (
        neiProcGlobalIDs_.size() != 0
     || neiProcVolumes_.size() != 0
    )
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    neiProcGlobalIDs_.setSize(mesh.boundaryMesh().size());
    neiProcVolumes_.setSize(mesh.boundaryMesh().size());

    PtrList<labelList>& neiProcGlobalIDs = neiProcGlobalIDs_;
    PtrList<scalarField>& neiProcVolumes = neiProcVolumes_;

    const scalarField& VI = mesh.V();

    // Send the data
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];

        if (fp.type() == "processor")
        {
            // Take a copy of the faceCells (local IDs) and convert them to
            // global IDs
            labelList globalFaceCells(fp.faceCells());
            foamPetscSnesHelper::globalCells().inplaceToGlobal(globalFaceCells);

            // Send global IDs to the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(fp);
            procPatch.send
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );

            // Construct a list of the volumes of cells at the patch
            scalarField patchVols(fp.size());
            const labelUList& faceCells = fp.faceCells();
            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                patchVols[faceI] = VI[cellID];
            }

            // Send patch cell volumes to the neighbour proc
            procPatch.send
            (
                Pstream::commsTypes::blocking, patchVols
            );
        }
    }

    // Receive the data
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];

        if (fp.type() == "processor")
        {
            neiProcGlobalIDs.set(patchI, new labelList(fp.size()));
            labelList& globalFaceCells = neiProcGlobalIDs[patchI];

            // Receive global IDs from the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(fp);
            procPatch.receive
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );

            // Receive patch cell volmes from the neighbour proc
            neiProcVolumes.set(patchI, new scalarField(fp.size()));
            scalarField& patchNeiVols = neiProcVolumes_[patchI];
            procPatch.receive
            (
                Pstream::commsTypes::blocking, patchNeiVols
            );
        }
    }
}


const PtrList<labelList>& foamPetscSnesHelper::neiProcGlobalIDs
(
    const fvMesh& mesh
) const
{
    if (neiProcGlobalIDs_.size() == 0)
    {
        makeNeiProcFields(mesh);
    }

    return neiProcGlobalIDs_;
}


const PtrList<scalarField>& foamPetscSnesHelper::neiProcVolumes
(
    const fvMesh& mesh
) const
{
    if (neiProcVolumes_.size() == 0)
    {
        makeNeiProcFields(mesh);
    }

    return neiProcVolumes_;
}


const leastSquaresS4fVectors& foamPetscSnesHelper::lsVectors
(
    const volScalarField& p
) const
{
    // Lookup and return vectors if they exist for this field
    // Otherwise create and return them
    const word lsName("leastSquaresVectors" + p.name());
    if (p.mesh().foundObject<leastSquaresS4fVectors>(lsName))
    {
        return p.mesh().lookupObject<leastSquaresS4fVectors>(lsName);
    }

    boolList useBoundaryFaceValues(p.mesh().boundary().size(), false);
    forAll(p.mesh().boundary(), patchI)
    {
        if (p.boundaryField()[patchI].fixesValue())
        {
            useBoundaryFaceValues[patchI] = true;
        }
    }

    return leastSquaresS4fVectors::New(lsName, p.mesh(), useBoundaryFaceValues);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamPetscSnesHelper::foamPetscSnesHelper
(
    fileName optionsFile,
    const label nBlocks,
    const label blockSize,
    //const labelPairList& nBlocksAndBlockSize,
    const labelListList& fieldDefs,
    const Switch stopOnPetscError,
    const Switch initialise
)
:
    foamPetscSnesHelper // Delegate constructor
    (
        optionsFile,
        labelPairList({{nBlocks, blockSize}}),
        fieldDefs,
        stopOnPetscError,
        initialise
    )
{}


foamPetscSnesHelper::foamPetscSnesHelper
(
    fileName optionsFile,
    const labelPairList& nBlocksAndBlockSize,
    const labelListList& fieldDefs,
    const Switch stopOnPetscError,
    const Switch initialise
)
:
    initialised_(initialise),
    globalCellsPtr_(),
    stopOnPetscError_(stopOnPetscError),
    snes_(nullptr),
    x_(nullptr),
    A_(nullptr),
    nRegions_(nBlocksAndBlockSize.size()),
    subMatsPtr_(nullptr),
    snesUserPtr_(),
    neiProcGlobalIDs_(),
    neiProcVolumes_()
{
    if (initialise)
    {
        // Expand the options file name
        optionsFile.expand();

        // Check the options file exists
        {
            IFstream is(optionsFile);
            if (!is.good())
            {
                FatalErrorInFunction
                    << "Cannot find the PETSc options file: " << optionsFile
                    << abort(FatalError);
            }
        }

        // Initialise PETSc with an options file
        PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

        // Create the PETSc SNES object
        snes_ = SNES();
        SNESCreate(PETSC_COMM_WORLD, &snes_);

        // Create user data context
        snesUserPtr_.set(new appCtxfoamPetscSnesHelper(*this));
        appCtxfoamPetscSnesHelper& user = snesUserPtr_();

        // Set the user context
        SNESSetApplicationContext(snes_, &user);

        // Set the residual function
        SNESSetFunction(snes_, NULL, formResidualFoamPetscSnesHelper, &user);

        if (nBlocksAndBlockSize.size() == 0)
        {
            FatalErrorInFunction
                << "nBlocksAndBlockSize should have a size greater than 0"
                << abort(FatalError);
        }

        // Count number of local blocks and local scalar equations
        label blockn = 0;
        label n = 0;
        forAll(nBlocksAndBlockSize, regionI)
        {
            const label nBlocksRegionI = nBlocksAndBlockSize[regionI].first();
            const label blockSizeRegionI = nBlocksAndBlockSize[regionI].second();
            blockn += nBlocksRegionI;
            n += nBlocksRegionI*blockSizeRegionI;
        }

        // Create globalCells object
        globalCellsPtr_.set(new globalIndex(blockn));

        // Global system size: total number of blocks across all processors
        //const int blockN = globalCellsPtr_->size();

        // Global system size: total number of scalar equation across all
        // processors
        const label N = returnReduce(n, sumOp<label>());

        // If there is only one region, create a normal matrix.
        if (nRegions_ == 1)
        {
            // Create the Jacobian matrix
            A_ = Mat();
            MatCreate(PETSC_COMM_WORLD, &A_);
            MatSetSizes(A_, n, n, N, N);
            //MatSetSizes(A, n, n, PETSC_DECIDE, PETSC_DECIDE);
            MatSetFromOptions(A_);
            MatSetType(A_, MATMPIAIJ);
            //MatSetType(A, MATMPIBAIJ);
        }
        else
        {
            // We have more than one region; create a nest matrix.

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
                    // Info<< "region = " << i << " " << j << endl;

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
            A_ = Mat();
            MatCreateNest
            (
                PETSC_COMM_WORLD,
                nRegions_,
                isRow.data(),
                nRegions_,
                isCol.data(),
                (Mat*)subMatsPtr_,
                &A_
            );

            // Cleanup: destroy the IS objects.
            for (label r = 0; r < nRegions_; ++r)
            {
                ISDestroy(&isRow[r]);
                ISDestroy(&isCol[r]);
            }
        }

        // Set the Jacobian function
        SNESSetJacobian(snes_, A_, A_, formJacobianFoamPetscSnesHelper, &user);

        // Set solver options
        // Uses default options, can be overridden by command line options
        SNESSetFromOptions(snes_);

        // Create the solution vector
        x_ = Vec();
        VecCreate(PETSC_COMM_WORLD, &x_);
        VecSetSizes(x_, n, N);
        if (nBlocksAndBlockSize.size() == 1)
        {
            // Only set the Vec block size for one region
            // Otherwise, it is not needed, as the Mat will be nested and get
            // its block sizes from its sub-matrices
            VecSetBlockSize(x_, nBlocksAndBlockSize[0].second());
        }
        VecSetType(x_, VECMPI);
        PetscObjectSetName((PetscObject) x_, "Solution");
        VecSetFromOptions(x_);
        VecZeroEntries(x_);

        // TO BE FIXED
        // Set up the fieldsplit index sets if fieldDefs is non-empty.
        // After setting the fieldsplit index sets, SNES/KSP will use them via
        // the PC.
        // We assume the ordering of the global vector is in blocks:
        // For each cell, the entries are:
        //     global_index = cell*blockSize + local_index,
        // where local_index runs from 0 to blockSize - 1.
        // The user provides fieldDefs such as, e.g.
        //     fieldDefs = { {0,1,2}, {3} }
        // to indicate that indices 0-2 in each block belong to field 0
        // (momentum) and index 3 belongs to field 1 (pressure).
        if (fieldDefs.size() > 0)
        {
            WarningInFunction
                << "fieldDefs not used: to be fixed!" << endl;

            // // Get the KSP from SNES, then retrieve the PC
            // KSP ksp;
            // SNESGetKSP(snes_, &ksp);
            // PC pc;
            // KSPGetPC(ksp, &pc);

            // // Obtain the local ownership range of the solution vector
            // PetscInt rstart, rend;
            // VecGetOwnershipRange(x, &rstart, &rend);
            // // Determine which cells (blocks) belong locally
            // const PetscInt localCellStart = rstart/blockSize;
            // const PetscInt localCellEnd = rend/blockSize;

            // // Loop over each field as defined in fieldDefs.
            // for
            // (
            //     label fieldIndex = 0;
            //     fieldIndex < fieldDefs.size();
            //     ++fieldIndex
            // )
            // {
            //     // Get the list of local indices for this field.
            //     const labelList& field = fieldDefs[fieldIndex];

            //     // Create an empty OpenFOAM List for the PETSc indices.
            //     DynamicList<PetscInt> indices(rend - rstart);

            //     // Loop over each locally owned cell
            //     for
            //     (
            //         PetscInt cell = localCellStart; cell < localCellEnd; ++cell
            //     )
            //     {
            //         // For each degree-of-freedom in this field...
            //         for (label j = 0; j < field.size(); j++)
            //         {
            //             PetscInt idx = cell*blockSize + field[j];

            //             // Only add the index if it is within our local ownership
            //             if (idx >= rstart && idx < rend)
            //             {
            //                 indices.append(idx);
            //             }
            //         }
            //     }

            //     // Create a PETSc index set from the indices.
            //     // Note: indices.begin() returns a pointer to the first element
            //     IS is;
            //     ISCreateGeneral
            //     (
            //         PETSC_COMM_WORLD,
            //         indices.size(),
            //         indices.begin(),
            //         PETSC_COPY_VALUES,
            //         &is
            //     );

            //     // Build a field name, "field0", "field1", etc.
            //     std::string fieldName = "field" + std::to_string(fieldIndex);
            //     PCFieldSplitSetIS(pc, fieldName.c_str(), is);
            //     ISDestroy(&is);
            // }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

foamPetscSnesHelper::~foamPetscSnesHelper()
{
    if (initialised_)
    {
        SNESDestroy(&snes_);

        VecDestroy(&x_);

        // Check if A is a nested matrix
        if (nRegions_ > 1)
        {
            // Destroy each submatrix
            for (label idx = 0; idx < nRegions_*nRegions_; ++idx)
            {
                MatDestroy(&subMatsPtr_[idx]);
            }

            // Free the array of Mat handles
            delete[] subMatsPtr_;
            subMatsPtr_ = nullptr;
        }

        MatDestroy(&A_);

        snesUserPtr_.clear();

        WarningInFunction
            << "Find a neat solution to finalise PETSc only once"
                << endl;
        //PetscFinalize();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
(
    const volScalarField& p,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const bool flipSign
) const
{
    // Get reference to least square vectors
    const fvMesh& mesh = p.mesh();
    const leastSquaresS4fVectors& lsv = lsVectors(p);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const scalarField& VI = mesh.V();

    const scalar sign = flipSign ? -1.0 : 1.0;

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise the block coefficient
    const label nCoeffCmpts = blockSize*blockSize;
    PetscScalar values[nCoeffCmpts];
    std::memset(values, 0, sizeof(values));

    forAll(own, faceI)
    {
        // Local block row ID
        const label ownCellID = own[faceI];

        // Local block column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            foamPetscSnesHelper::globalCells().toGlobal(neiCellID);

        // Explicit gradient
        // lsGrad[ownCellID] += ownLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);

        // Insert coefficients for block row globalBlockRowI

        // lsGrad[ownCellID] += ownLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // mat(ownCellID, ownCellID) -= VI[ownCellID]*ownLs[faceI];
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            values[cmptI*blockSize + colOffset] =
                -sign*VI[ownCellID]*ownLs[faceI][cmptI];
        }
        CHKERRQ
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                ADD_VALUES
            )
        );

        // lsGrad[ownCellID] += ownLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // mat(ownCellID, neIFaceI) += VI[ownCellID]*ownLs[faceI];
        // Flip the sign of the values
        for (label i = 0; i < nCoeffCmpts; ++i)
        {
            values[i] = -values[i];
        }
        CHKERRQ
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockRowI, 1, &globalBlockColI, values,
                ADD_VALUES
            )
        );

        // Insert coefficients for block row globalBlockColI

        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        //mat(neIFaceI, neIFaceI) += VI[neiCellID]*neiLs[faceI];
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            values[cmptI*blockSize + colOffset] =
                -sign*VI[neiCellID]*neiLs[faceI][cmptI];
        }
        CHKERRQ
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockColI, 1, &globalBlockColI, values,
                ADD_VALUES
            )
        );

        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        //mat(neIFaceI, ownCellID) -= VI[neiCellID]*neiLs[faceI];
        // Flip the sign of the values
        for (label i = 0; i < nCoeffCmpts; ++i)
        {
            values[i] = -values[i];
        }
        CHKERRQ
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockColI, 1, &globalBlockRowI, values,
                ADD_VALUES
            )
        );
    }

    // Boundary face contributions
    //const boolList& useBoundaryFaceValues = lsv.useBoundaryFaceValues();
    const boolList useBoundaryFaceValues(mesh.boundary().size(), true);
    forAll(mesh.boundary(), patchI)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchI];
        const labelUList& faceCells = mesh.boundary()[patchI].faceCells();
        const fvPatch& fp = mesh.boundary()[patchI];
        if (fp.type() == "processor")
        {
            const labelList& neiGlobalFaceCells =
                neiProcGlobalIDs(mesh)[patchI];
            const scalarField& patchNeiVols = neiProcVolumes(mesh)[patchI];
            const vectorField& patchNeiLs = neiLs.boundaryField()[patchI];
            // const vectorField patchNeiLs(... patchNeighbourField());

            forAll(fp, patchFaceI)
            {
                // Local block row ID
                const label ownCellID = faceCells[patchFaceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                // On-proc diagonal coefficient
                // mat(ownCellID, ownCellID) -= VI[ownCellID]*ownLs[faceI];
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize + colOffset] =
                        -sign*VI[ownCellID]*patchOwnLs[patchFaceI][cmptI];
                }
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                        ADD_VALUES
                    );
                );

                // Neighbour global cell ID
                const label globalBlockColI = neiGlobalFaceCells[patchFaceI];

                // Off-proc off-diagonal coefficient
                // mat(ownCellID, neiCellID) += VI[ownCellID]*neiLs[faceI];
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize + colOffset] =
                        sign*patchNeiVols[patchFaceI]
                       *patchNeiLs[patchFaceI][cmptI];
                }
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockColI, values,
                        ADD_VALUES
                    );
                );
            }
        }
        else if (fp.coupled()) // coupled but not processor
        {
            FatalErrorInFunction
                << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        else if
        (
            isA<symmetryPolyPatch>(fp.patch())
         || isA<symmetryPlanePolyPatch>(fp.patch())
        )
        {
            // The delta in scalar p across the symmetry is zero by definition
            // so the symmetry plane does not contribution coefficients
        }
        else
        {
            if (useBoundaryFaceValues[patchI])
            {
                forAll(faceCells, patchFaceI)
                {
                    // Explicit calculation
                    // lsGrad[faceCells[patchFaceI]] +=
                    //      patchOwnLs[patchFaceI]
                    //     *(patchVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);

                    // Subtract patchOwnLs[patchFaceI] from (faceCellI, faceCellI)
                    // Nothing else to do as patch value is a known value

                    // Local block row ID
                    const label ownCellID = faceCells[patchFaceI];

                    // Global block row ID
                    const label globalBlockRowI =
                        foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                    // mat(ownCellID, ownCellID) -=
                    //     VI[ownCellID]*patchOwnLs[patchFaceI];
                    for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                    {
                        values[cmptI*blockSize + colOffset] =
                            -sign*VI[ownCellID]*patchOwnLs[patchFaceI][cmptI];
                    }
                    CHKERRQ
                    (
                        MatSetValuesBlocked
                        (
                            jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                            ADD_VALUES
                        )
                    );
                }
            }
        }
    }

    return 0;
}


label foamPetscSnesHelper::InsertFvmDivPhiUIntoPETScMatrix
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const bool flipSign
) const
{
    // Take references for efficiency and brevity
    const fvMesh& mesh = U.mesh();
    const scalarField& phiI = phi;
    const vectorField& UI = U;
    const vectorField& SfI = mesh.Sf();
    const scalarField& wI = mesh.weights();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const scalar sign = flipSign ? -1 : 1;
    tensor coeff = tensor::zero;

    // Add segregated component of convection term
    {
        fvVectorMatrix UEqn
        (
          - fvm::div(phi, U, "jacobian-div(phi,U)")
        );

        UEqn.relax();

        // Convert fvMatrix matrix to PETSc matrix
        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
        (
            UEqn, jac, 0, 0, nScalarEqns
        );
    }

    // Add additional term from div(phi,U) upwind linearisation
    // (d/dU)(\sum phi*U = phi*I + w*Sf*Up
    // for phi > 0, where phi = Sf*(w*Up + (1-w)*Un)
    // The w*Sf*Up term is missing from fvm::div(phi(), U) as it requires a
    // coupled solver. So we will add the w*Sf*Up term here
    // Similary, for phi < 0, the neighbour coefficient is (1 - w)*Sf*Un

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Prepare coeff array
    const label nCoeffCmpts = blockSize*blockSize;
    PetscScalar values[nCoeffCmpts];
    std::memset(values, 0, sizeof(values));

    forAll(phiI, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            foamPetscSnesHelper::globalCells().toGlobal(neiCellID);

        if (phiI[faceI] > 0)
        {
            // Add w*Sf*Up to owner eqn
            // coeff = sign*wI[faceI]*SfI[faceI]*UI[ownCellID];
            coeff = sign*wI[faceI]*UI[ownCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                    ADD_VALUES
                )
            );

            // Add (1 - w)*Sf*Up as nei to contribution to own eqn
            // coeff = sign*(1.0 - wI[faceI])*SfI[faceI]*UI[ownCellID];
            coeff = sign*(1.0 - wI[faceI])*UI[ownCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI, values,
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Sf*Up to the neighbour diagonal
            // Flip the sign
            // coeff = -sign*(1.0 - wI[faceI])*SfI[faceI]*UI[ownCellID];
            coeff = -sign*(1.0 - wI[faceI])*UI[ownCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockColI, values,
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Sf*Up as the own contribution to the nei eqn
            // Flip the sign
            // coeff = -sign*(1.0 - wI[faceI])*SfI[faceI]*UI[ownCellID];
            coeff = -sign*(1.0 - wI[faceI])*UI[ownCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI, values,
                    ADD_VALUES
                )
            );
        }
        else
        {
            // Add w*Sf*Un to owner diagonal
            // coeff = sign*wI[faceI]*SfI[faceI]*UI[neiCellID];
            coeff = sign*wI[faceI]*UI[neiCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                    ADD_VALUES
                )
            );

            // Add w*Sf*Un as nei to contribution to own eqn
            // coeff = sign*wI[faceI]*SfI[faceI]*UI[neiCellID];
            coeff = sign*wI[faceI]*UI[neiCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI, values,
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Sf*Un to neighbour diagonal
            // coeff = -sign*(1.0 - wI[faceI])*SfI[faceI]*UI[neiCellID];
            coeff = -sign*(1.0 - wI[faceI])*UI[neiCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockColI, values,
                    ADD_VALUES
                )
            );

            // Add -w*Sf*Un as own to contribution to nei eqn
            // coeff = -sign*wI[faceI]*SfI[faceI]*UI[neiCellID];
            coeff = -sign*wI[faceI]*UI[neiCellID]*SfI[faceI];
            for (label i = 0; i < (blockSize - 1); ++i)
            {
                for (label j = 0; j < (blockSize - 1); ++j)
                {
                    // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                    // the 4x4 (or 3x3 in 2-D) values matrix
                    values[(i + rowOffset)*blockSize + j + colOffset] =
                        coeff[i*3 + j];
                }
            }
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI, values,
                    ADD_VALUES
                )
            );
        }
    }

    // Boundary contribution for outlets
    forAll(mesh.boundary(), patchI)
    {
        const fvPatchVectorField& pU = U.boundaryField()[patchI];
        const fvPatch& patch = pU.patch();
        const vectorField& pSf = patch.Sf();
        //const scalarField& pphi = phi().boundaryField()[patchI];
        const scalarField& pw = patch.weights();
        const labelUList& fc = patch.faceCells();

        if (patch.coupled())
        {
            notImplemented("div(phi,U) additional term for processors");
        }
        else if (patch.type() != "empty" && !pU.fixesValue())
        {
            forAll(fc, faceI)
            {
                // Local row ID
                const label ownCellID = fc[faceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                // Add w*Sf*Up to owner eqn
                coeff = sign*pw[faceI]*UI[ownCellID]*pSf[faceI];
                for (label i = 0; i < (blockSize - 1); ++i)
                {
                    for (label j = 0; j < (blockSize - 1); ++j)
                    {
                        // Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
                        // the 4x4 (or 3x3 in 2-D) values matrix
                        values[(i + rowOffset)*blockSize + j + colOffset] =
                            coeff[i*3 + j];
                    }
                }
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                        ADD_VALUES
                    )
                );
            }
        }
    }

    return 0;
}


label foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
(
    const volScalarField& p,
    const volVectorField& U,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const bool flipSign
) const
{
    // Take references for efficiency and brevity
    const fvMesh& mesh = p.mesh();
    const scalar sign = flipSign ? -1 : 1;

    for (label cmptI = 0; cmptI < 3; ++cmptI)
    {
        if (nScalarEqns == 2 && cmptI == 2)
        {
            break;
        }

        fvScalarMatrix divUCoeffs(p, dimArea*dimPressure);

        const vectorField& Sf = mesh.Sf();
        const surfaceScalarField& weights = mesh.weights();
        const scalarField& w = weights;

        scalarField& upper = divUCoeffs.upper();
        scalarField& lower = divUCoeffs.lower();

        lower = w*Sf.component(cmptI);
        upper = lower - Sf.component(cmptI);

        divUCoeffs.negSumDiag();

        // Boundary contributions
        scalarField& diag = divUCoeffs.diag();
        forAll(mesh.boundary(), patchI)
        {
            const fvPatchVectorField& pU = U.boundaryField()[patchI];
            const fvPatch& patch = pU.patch();
            const vectorField& Sf = patch.Sf();
            const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
            const labelUList& fc = patch.faceCells();

            const vectorField internalCoeffs(pU.valueInternalCoeffs(pw));

            // Diag contribution
            forAll(pU, faceI)
            {
                diag[fc[faceI]] -=
                    sign*internalCoeffs[faceI][cmptI]*Sf[faceI][cmptI];
            }

            if (patch.coupled())
            {
                // Todo: add off-core coeffs
                notImplemented("patch.coupled(): D-in-p");

                // CoeffField<vector>::linearTypeField& pcoupleUpper =
                //     bs.coupleUpper()[patchI].asLinear();
                // CoeffField<vector>::linearTypeField& pcoupleLower =
                //     bs.coupleLower()[patchI].asLinear();

                // const vectorField pcl = -pw*Sf;
                // const vectorField pcu = pcl + Sf;

                // // Coupling  contributions
                // pcoupleLower -= pcl;
                // pcoupleUpper -= pcu;
            }
        }

        // Insert component coeffs
        // cmptI is the 1st column for Ux, 2nd for Uy, 3rd for Uz
        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
        (
            divUCoeffs, jac, rowOffset, colOffset + cmptI, 1
        );
    }

    return 0;
}


int foamPetscSnesHelper::solve(const bool returnOnSnesError)
{
    if (!initialised_)
    {
        FatalErrorIn("void foamPetscSnesHelper::solve()")
            << "This function cannot be called because the foamPetscSnesHelper "
            << "object was not initialised during construction"
            << abort(FatalError);
    }

    // Solve the nonlinear system
    SNESSolve(snes_, NULL, x_);

    // Check convergence
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes_, &reason);

    if (reason < 0)
    {
        WarningIn
        (
            "void foamPetscSnesHelper::checkConvergence(SNES snes) const"
        )   << "PETSc SNES solver return error check disabled" << endl
            << "The SNES nonlinear solver did not converge." << nl
            << " PETSc SNES convergence error code: " << reason << nl
            << " PETSc SNES convergence reason: "
            << SNESConvergedReasons[reason] << endl;

        if (returnOnSnesError)
        {
            return reason;
        }
        else if (stopOnPetscError_)
        {
            FatalErrorIn
            (
                "void foamPetscSnesHelper::checkConvergence(SNES snes) const"
            )   << "Stopping because of the PETSc SNES error" << nl
                << "Set `stopOnPetscError` to `false` to continue on PETSc "
                << "SNES errors"
                << abort(FatalError);
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
