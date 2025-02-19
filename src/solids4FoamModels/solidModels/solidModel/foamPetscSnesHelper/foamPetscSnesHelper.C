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

    // Map the solution to an OpenFOAM field
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;
    const int blockSize = user->solMod_.blockSize();

    // Initialise the matrix if it has yet to be allocated; otherwise zero all
    // entries
    MatInfo info;
    MatGetInfo(B, MAT_LOCAL, &info);
    // PetscInt rstart, rend;
    // MatGetOwnershipRange(B, &rstart, &rend);
    if (info.nz_used)
    {
        // Zero the matrix but do not reallocate the space
        // The "-snes_lag_jacobian -2" PETSc option can be used to avoid
        // re-building the matrix
        CHKERRQ(MatZeroEntries(B));
    }
    else
    {
        Foam::Info<< "    Initialising the matrix..." << Foam::flush;

        // Set the block size
        CHKERRQ(MatSetBlockSize(B, blockSize));

        // Number of vector unknowns on this processor
        const int blockn = user->solMod_.globalCells().localSize();

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
        const Foam::fvMesh& mesh = user->solMod_.fmesh();

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
        //CHKERRQ(MatMPIAIJSetPreallocation(B, 0, d_nnz, 0, o_nnz));
        // Allocate parallel matrix with the same conservative stencil per node
        //CHKERRQ(MatMPIAIJSetPreallocation(B, d_nz, NULL, 0, NULL));
        CHKERRQ(MatMPIBAIJSetPreallocation(B, blockSize, 0, d_nnz, 0, o_nnz));

        // Raise an error if mallocs are required during matrix assembly
        MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(foamPetscSnesHelper, 0);


// * * * * * * * * * * * * * * * Private Function  * * * * * * * * * * * * * //

void foamPetscSnesHelper::makeNeiProcFields() const
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

    const fvMesh& mesh = mesh_;

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


const PtrList<labelList>& foamPetscSnesHelper::neiProcGlobalIDs() const
{
    if (neiProcGlobalIDs_.size() == 0)
    {
        makeNeiProcFields();
    }

    return neiProcGlobalIDs_;
}


const PtrList<scalarField>& foamPetscSnesHelper::neiProcVolumes() const
{
    if (neiProcVolumes_.size() == 0)
    {
        makeNeiProcFields();
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

    boolList useBoundaryFaceValues(mesh_.boundary().size(), false);
    forAll(mesh_.boundary(), patchI)
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
    const fvMesh& mesh,
    const label blockSize,
    const labelListList& fieldDefs,
    const Switch stopOnPetscError,
    const Switch initialise
)
:
    initialised_(initialise),
    mesh_(mesh),
    globalCells_(mesh.nCells()),
    blockSize_(blockSize),
    stopOnPetscError_(stopOnPetscError),
    snesPtr_(),
    xPtr_(),
    APtr_(),
    snesUserPtr_(),
    neiProcGlobalIDs_(),
    neiProcVolumes_()
{
    if (initialise)
    {
        // Expand the options file name
        optionsFile.expand();

        // Initialise PETSc with an options file
        PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

        // Create the PETSc SNES object
        snesPtr_.set(new SNES());
        SNES& snes = snesPtr_();
        SNESCreate(PETSC_COMM_WORLD, &snes);

        // Create user data context
        snesUserPtr_.set(new appCtxfoamPetscSnesHelper(*this));
        appCtxfoamPetscSnesHelper& user = snesUserPtr_();

        // Set the user context
        SNESSetApplicationContext(snes, &user);

        // Set the residual function
        SNESSetFunction(snes, NULL, formResidualFoamPetscSnesHelper, &user);

        // Global system size
        const int blockN = globalCells_.size();
        const int N = blockSize_*blockN;

        // Local (this processor) system size
        const int blockn = globalCells_.localSize();
        const int n = blockSize_*blockn;

        // Create the Jacobian matrix
        APtr_.set(new Mat());
        Mat& A = APtr_();
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, n, n, N, N);
        MatSetFromOptions(A);
        MatSetType(A, MATMPIAIJ);
        //MatSetType(A, MATMPIBAIJ);

        // Set the Jacobian function
        SNESSetJacobian(snes, A, A, formJacobianFoamPetscSnesHelper, &user);

        // Set solver options
        // Uses default options, can be overridden by command line options
        SNESSetFromOptions(snes);

        // Create the solution vector
        xPtr_.set(new Vec());
        Vec& x = xPtr_();
        VecCreate(PETSC_COMM_WORLD, &x);
        VecSetSizes(x, n, N);
        VecSetBlockSize(x, blockSize_);
        VecSetType(x, VECMPI);
        PetscObjectSetName((PetscObject) x, "Solution");
        VecSetFromOptions(x);
        VecZeroEntries(x);

        // Set up the fieldsplit index sets if fieldDefs is non-empty.
        // After setting the fieldsplit index sets, SNES/KSP will use them via
        // the PC.
        // We assume the ordering of the global vector is in blocks:
        // For each cell, the entries are:
        //     global_index = cell*blockSize_ + local_index,
        // where local_index runs from 0 to blockSize_-1.
        // The user provides fieldDefs such as, e.g.
        //     fieldDefs = { {0,1,2}, {3} }
        // to indicate that indices 0-2 in each block belong to field 0
        // (momentum) and index 3 belongs to field 1 (pressure).
        if (fieldDefs.size() > 0)
        {
            // Get the KSP from SNES, then retrieve the PC
            KSP ksp;
            SNESGetKSP(snes, &ksp);
            PC pc;
            KSPGetPC(ksp, &pc);

            // Obtain the local ownership range of the solution vector
            PetscInt rstart, rend;
            VecGetOwnershipRange(x, &rstart, &rend);
            // Determine which cells (blocks) belong locally
            const PetscInt localCellStart = rstart/blockSize_;
            const PetscInt localCellEnd = rend/blockSize_;

            // Loop over each field as defined in fieldDefs.
            for
            (
                label fieldIndex = 0;
                fieldIndex < fieldDefs.size();
                ++fieldIndex
            )
            {
                // Get the list of local indices for this field.
                const labelList& field = fieldDefs[fieldIndex];

                // Create an empty OpenFOAM List for the PETSc indices.
                DynamicList<PetscInt> indices(rend - rstart);

                // Loop over each locally owned cell
                for
                (
                    PetscInt cell = localCellStart; cell < localCellEnd; ++cell
                )
                {
                    // For each degree-of-freedom in this field...
                    for (label j = 0; j < field.size(); j++)
                    {
                        PetscInt idx = cell*blockSize_ + field[j];

                        // Only add the index if it is within our local ownership
                        if (idx >= rstart && idx < rend)
                        {
                            indices.append(idx);
                        }
                    }
                }

                // Create a PETSc index set from the indices.
                // Note: indices.begin() returns a pointer to the first element
                IS is;
                ISCreateGeneral
                (
                    PETSC_COMM_WORLD,
                    indices.size(),
                    indices.begin(),
                    PETSC_COPY_VALUES,
                    &is
                );

                // Build a field name, "field0", "field1", etc.
                std::string fieldName = "field" + std::to_string(fieldIndex);
                PCFieldSplitSetIS(pc, fieldName.c_str(), is);
                ISDestroy(&is);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

foamPetscSnesHelper::~foamPetscSnesHelper()
{
    if (initialised_)
    {
        SNESDestroy(&snesPtr_());
        snesPtr_.clear();

        VecDestroy(&xPtr_());
        xPtr_.clear();

        MatDestroy(&APtr_());
        APtr_.clear();

        snesUserPtr_.clear();

        PetscFinalize();
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
    const fvMesh& mesh = mesh_;
    const leastSquaresS4fVectors& lsv = lsVectors(p);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const scalarField& VI = mesh.V();

    const scalar sign = flipSign ? -1.0 : 1.0;

    // Initialise the block coefficient
    const label nCoeffCmpts = blockSize_*blockSize_;
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
            values[cmptI*blockSize_ + colOffset] =
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
            values[cmptI*blockSize_ + colOffset] =
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
            const labelList& neiGlobalFaceCells = neiProcGlobalIDs()[patchI];
            const scalarField& patchNeiVols = neiProcVolumes()[patchI];
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
                    values[cmptI*blockSize_ + colOffset] =
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
                    values[cmptI*blockSize_ + colOffset] =
                        patchNeiVols[patchFaceI]*patchNeiLs[patchFaceI][cmptI];
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
            isA<symmetryPolyPatch>(fp) || isA<symmetryPlanePolyPatch>(fp)
        )
        {
            // Treat symmetry planes consistently with internal faces
            // Use the mirrored face-cell values rather than the patch face
            // values
            // See https://doi.org/10.1080/10407790.2022.2105073
            const vectorField nHat(fp.nf());
            forAll(faceCells, patchFaceI)
            {
                // Explicit calculation
                // lsGrad[faceCells[patchFaceI]] +=
                //      patchOwnLs[patchFaceI]
                //     *(
                //         transform
                //         (
                //             I - 2.0*sqr(nHat[patchFaceI]),
                //             vsf[faceCells[patchFaceI]]
                //         )
                //       - vsf[faceCells[patchFaceI]]
                //     );

                // Subtract "patchOwnLs[patchFaceI] & 2*sqr(nHat[patchFaceI])"
                // from (faceCells[patchFaceI], faceCells[patchFaceI])

                // Local block row ID
                const label ownCellID = faceCells[patchFaceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal(ownCellID);

                // mat(ownCellID, ownCellID) -=
                //     VI[ownCellID]*patchOwnLs[patchFaceI]
                //   & 2*sqr(nHat[patchFaceI])
                const vector coeff
                (
                    -sign*VI[ownCellID]*patchOwnLs[patchFaceI]
                   & (2.0*sqr(nHat[patchFaceI]) )
                );
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize_ + colOffset] = coeff[cmptI];
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
                        values[cmptI*blockSize_ + colOffset] =
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


int foamPetscSnesHelper::solve(const bool returnOnSnesError)
{
    if (!initialised_)
    {
        FatalErrorIn("void foamPetscSnesHelper::solve()")
            << "This function cannot be called because the foamPetscSnesHelper "
            << "object was not initialised during construction"
            << abort(FatalError);
    }

    // Take a reference to the SNES solver object
    SNES& snes = snesPtr_();

    // Take a reference to the solution vector
    Vec& x = xPtr_();

    // Solve the nonlinear system
    SNESSolve(snes, NULL, x);

    // Check convergence
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes, &reason);

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
