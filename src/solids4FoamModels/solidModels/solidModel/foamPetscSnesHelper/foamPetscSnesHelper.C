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
    // Note
    // The "-snes_lag_jacobian -2" PETSc option can be used to avoid
    // re-building the matrix

    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Access the OpenFOAM data
    appCtxfoamPetscSnesHelper *user = (appCtxfoamPetscSnesHelper *)ctx;

    // Zero the matrix but do not reallocate the space
    // For a nested matrix, this will zero all sub-matrices
    CHKERRQ(MatZeroEntries(B));

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
    Mat& jac,
    const fvMesh& mesh,
    const label blockSize,
    const bool createMat
)
{
    // Count number of local blocks and local scalar equations
    const label blockn = mesh.nCells();
    const label n = blockn*blockSize;

    // Global system size: total number of scalar equation across all
    // processors
    const label N = returnReduce(n, sumOp<label>());

    // // Set the Jacobian matrix size
    if (createMat)
    {
        MatCreate(PETSC_COMM_WORLD, &jac);
        MatSetFromOptions(jac);
        MatSetSizes(jac, n, n, N, N);
        MatSetType(jac, MATMPIAIJ);
    }

    // Set the block size
    CHKERRQ(MatSetBlockSize(jac, blockSize));

    // Count the number of non-zeros in the matrix
    // Note: we assume a compact stencil, i.e. face only face neighbours

    // Number of on-processor non-zeros per row
    int* d_nnz = (int*)malloc(blockn*sizeof(int));

    // Number of off-processor non-zeros per row
    int* o_nnz = (int*)malloc(blockn*sizeof(int));

    // Initialise d_nnz to one and o_nnz to zero
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


Foam::label Foam::initialiseSolution
(
    Vec& x,
    const fvMesh& mesh,
    const label blockSize,
    const bool createVec
)
{
    if (createVec)
    {
        // Count number of local blocks and local scalar equations
        const label blockn = mesh.nCells();
        const label n = blockn*blockSize;

        // Global system size: total number of scalar equation across all
        // processors
        const label N = returnReduce(n, sumOp<label>());

        x = Vec();
        CHKERRQ(VecCreate(PETSC_COMM_WORLD, &x));
        CHKERRQ(VecSetSizes(x, n, N));
        CHKERRQ(VecSetType(x, VECMPI));
    }

    CHKERRQ(VecSetBlockSize(x, blockSize));
    CHKERRQ(PetscObjectSetName((PetscObject) x, "Solution"));
    CHKERRQ(VecSetFromOptions(x));
    CHKERRQ(VecZeroEntries(x));

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


label foamPetscSnesHelper::initialiseSnes()
{
    if (snes_)
    {
        FatalErrorInFunction
            << "Pointer already set" << abort(FatalError);
    }

    // Create the PETSc SNES object
    snes_ = SNES();
    CHKERRQ(SNESCreate(PETSC_COMM_WORLD, &snes_));

    // Create user data context
    snesUserPtr_.set(new appCtxfoamPetscSnesHelper(*this));
    appCtxfoamPetscSnesHelper& user = snesUserPtr_();

    // Set the user context
    CHKERRQ(SNESSetApplicationContext(snes_, &user));

    // Set the residual function
    CHKERRQ
    (
        SNESSetFunction(snes_, NULL, formResidualFoamPetscSnesHelper, &user)
    );

    // The derived class initialises A
    CHKERRQ(initialiseJacobian(A_));

    // Set the Jacobian function
    CHKERRQ
    (
        SNESSetJacobian(snes_, A_, A_, formJacobianFoamPetscSnesHelper, &user)
    );

    // Set solver options
    // Uses default options, can be overridden by command line options
    CHKERRQ(SNESSetFromOptions(snes_));

    // The derived class initialises the solution vector
    CHKERRQ(initialiseSolution(x_));

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


foamPetscSnesHelper::foamPetscSnesHelper
(
    fileName optionsFile,
    const label nLocalBlocks,
    const Switch stopOnPetscError,
    const Switch initialise
)
:
    initialised_(initialise),
    options_(nullptr),
    stopOnPetscError_(stopOnPetscError),
    snes_(nullptr),
    x_(nullptr),
    A_(nullptr),
    snesUserPtr_(),
    globalCellsPtr_
    (
        nLocalBlocks >= 0 ? new globalIndex(nLocalBlocks) : nullptr
    ),
    neiProcGlobalIDs_(),
    neiProcVolumes_()
{
    if (initialise)
    {
        // Expand the options file name
        optionsFile.expand();

        // Check the options file exists
        IFstream is(optionsFile);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot find the PETSc options file: " << optionsFile
                << abort(FatalError);
        }

        // Initialise PETSc with an options file, but we will need to reset it
        // before each call to snes.solve as PETSc actively uses only one
        // options database and there may be several objects using PETSc, e.g.
        // solid solver, fluid solver, mesh motion
        PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

        // Create and store the options database
        PetscOptionsCreate(&options_);
        PetscOptionsInsertFile
        (
            PETSC_COMM_WORLD, options_, optionsFile.c_str(), PETSC_TRUE
        );
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

foamPetscSnesHelper::~foamPetscSnesHelper()
{
    if (initialised_)
    {
        PetscOptionsDestroy(&options_);
        SNESDestroy(&snes_);
        VecDestroy(&x_);
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
            coeff = -sign*(wI[faceI])*UI[ownCellID]*SfI[faceI];
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
        else if
        (
            patch.type() != "empty"
         && !isA<symmetryFvPatchField<vector>>(pU)
         && !isA<symmetryPlaneFvPatchField<vector>>(pU)
         && !pU.fixesValue()
        )
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

    // Initialise the SNES object
    if (!snes_)
    {
        if (initialiseSnes() != 0)
        {
            FatalErrorInFunction
                << "initialiseSnes failed" << abort(FatalError);
        }
    }

    // Load the correct options database
    PetscOptionsPush(options_);
    SNESSetFromOptions(snes_);

    // Solve the nonlinear system
    SNESSolve(snes_, NULL, x_);

    // Un-load the options file
    PetscOptionsPop();

    // Check convergence
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes_, &reason);

    if (reason < 0)
    {
        WarningInFunction
            << "PETSc SNES solver return error check disabled" << endl
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
            FatalErrorInFunction
                << "Stopping because of the PETSc SNES error" << nl
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
