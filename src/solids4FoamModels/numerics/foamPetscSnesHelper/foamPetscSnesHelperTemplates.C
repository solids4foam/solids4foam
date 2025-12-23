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
#include "processorFvPatch.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class Type>
label foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
(
    fvMatrix<Type>& fvM,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns
) const
{
    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise block coeff
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    // Get the number of components per field element from the traits (e.g.,
    // 1 for scalar, 3 for vector, etc.)
    //const label vfBlockSize = pTraits<Type>::nComponents;

    // Take a reference to the mesh
    const fvMesh& mesh = fvM.psi().mesh();

    // Insert the diagonal
    if (fvM.hasDiag())
    {
        // Retrieve the matrix diagonal
        const Field<Type> diag(fvM.DD());

        // Insert the coeffs
        forAll(diag, blockRowI)
        {
            // Obtain a pointer to the underlying scalar data in
            // diag[blockRowI] (assumes contiguous storage)
            const scalar* curDiagPtr =
                reinterpret_cast<const scalar*>(&diag[blockRowI]);

            // Construct the diag block coefficient
            for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
            {
                for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                {
                    if (cmptI == cmptJ)
                    {
                        values
                        [
                            (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                        ] = curDiagPtr[cmptI];
                    }
                }
            }

            // Get the global row index
            const label globalBlockRowI =
                foamPetscSnesHelper::globalCells().toGlobal(blockRowI);

            // Insert block coefficient
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );
        }
    }

    // Insert the upper off-diagonal
    const bool hasUpper = fvM.hasUpper();
    const bool hasLower = fvM.hasLower();
    if (hasUpper)
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& upper = fvM.upper();

        forAll(own, faceI)
        {
            // Construct the upper block coefficient
            const scalar upperScalarCoeff = upper[faceI];
            for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
            {
                for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                {
                    if (cmptI == cmptJ)
                    {
                        values
                        [
                            (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                        ] = upperScalarCoeff;
                    }
                }
            }

            const label blockRowI = own[faceI];
            const label blockColI = nei[faceI];

            const label globalBlockRowI =
                foamPetscSnesHelper::globalCells().toGlobal(blockRowI);
            const label globalBlockColI =
                foamPetscSnesHelper::globalCells().toGlobal(blockColI);

            if (hasUpper)
            {
                // Insert upper block coeff
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );
            }

            if (!hasLower) // symmetric matrix: use the upper coeff
            {
                // Insert symmetric lower block coeff
                CHKERRQ
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockColI, 1, &globalBlockRowI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );
            }
        }
    }

    if (hasLower)
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& lower = fvM.lower();

        forAll(own, faceI)
        {
            // Construct the lower block coefficient
            const scalar lowerScalarCoeff = lower[faceI];
            for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
            {
                for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                {
                    if (cmptI == cmptJ)
                    {
                        values
                        [
                            (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                        ] = lowerScalarCoeff;
                    }
                }
            }

            const label blockRowI = own[faceI];
            const label blockColI = nei[faceI];

            const label globalBlockRowI =
                foamPetscSnesHelper::globalCells().toGlobal(blockRowI);
            const label globalBlockColI =
                foamPetscSnesHelper::globalCells().toGlobal(blockColI);

            // Insert lower coeff
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );
        }
    }

    // Get the global IDs of the neighbour cells across processor patches
    const PtrList<labelList>& neiProcGlobalIDs = this->neiProcGlobalIDs(mesh);

    // Insert the off-processor coefficients
    forAll(mesh.boundaryMesh(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];
        if (fp.type() == "processor")
        {
            const Field<Type>& intCoeffs = fvM.internalCoeffs()[patchI];
            const Field<Type>& neiCoeffs = fvM.boundaryCoeffs()[patchI];
            const labelUList& faceCells = mesh.boundary()[patchI].faceCells();
            const labelList& neiGlobalFaceCells = neiProcGlobalIDs[patchI];

            forAll(fp, faceI)
            {
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal
                    (
                        faceCells[faceI]
                    );

                // On-proc diagonal coefficient
                {
                    // Obtain a pointer to the underlying scalar data in
                    // intCoeffs[faceI] (assumes contiguous storage)
                    const scalar* curIntCoeffsPtr =
                        reinterpret_cast<const scalar*>
                        (
                            &intCoeffs[faceI]
                        );

                    for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                    {
                        for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                        {
                            if (cmptI == cmptJ)
                            {
                                values
                                [
                                    (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                                ] = curIntCoeffsPtr[cmptI];
                            }
                        }
                    }

                    CHKERRQ
                    (
                        MatSetValuesBlocked
                        (
                            jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                            values.cdata(),
                            ADD_VALUES
                        )
                    );
                }

                // Off-proc off-diagonal coefficient
                {
                    // Obtain a pointer to the underlying scalar data in
                    // neiCoeffs[faceI] (assumes contiguous storage)
                    const scalar* curNeiCoeffsPtr =
                        reinterpret_cast<const scalar*>
                        (
                            &neiCoeffs[faceI]
                        );

                    for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                    {
                        for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                        {
                            if (cmptI == cmptJ)
                            {
                                // Take care: we need to flip the sign
                                values
                                [
                                    (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                                ] = -curNeiCoeffsPtr[cmptI];
                            }
                        }
                    }

                    const label globalBlockColI = neiGlobalFaceCells[faceI];

                    CHKERRQ
                    (
                        MatSetValuesBlocked
                        (
                            jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                            values.cdata(),
                            ADD_VALUES
                        )
                    );
                }
            }
        }
        else if (fp.coupled()) // coupled but not a processor boundary
        {
            FatalErrorInFunction
                << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        // else non-coupled boundary contributions have already been added to
        // the diagonal
    }

    return 0;
}


template <class Type>
void foamPetscSnesHelper::ExtractFieldComponents
(
    const PetscScalar *x,
    Field<Type>& vf,
    const label offset,
    const label xBlockSize,
    const labelList& compIndices
) const
{
    // Obtain a pointer to the underlying scalar array in vf
    // We assume vf is stored as contiguous data, i.e.
    // (vx1, vy1, vz1, vx2, vy2, vz2, ..., vxN, vyN, vzN)
    scalar* vfPtr = reinterpret_cast<scalar*>(vf.begin());

    // Get the number of components per field element from the traits (e.g.,
    // 1 for scalar, 3 for vector, etc.)
    const label vfBlockSize = pTraits<Type>::nComponents;

    // Decide which components to extract: if compIndices is empty, default to
    // all components
    labelList compList;
    if (compIndices.size() == 0)
    {
        compList.setSize(vfBlockSize);
        for (label i = 0; i < vfBlockSize; ++i)
        {
            compList[i] = i;  // Default: copy all components in order
        }
    }
    else
    {
        compList = compIndices;
    }

    // The number of components to insert per block is determined by the size of
    // compList
    const label nCompToInsert = compList.size();

    if (location_ == solutionLocation::CELLS)
    {
        // Loop over each block in vf and copy the selected components in x into the
        // appropriate positions in vf.
        // The source index for cell 'cellI' is computed as:
        //   offset + (cellI*blockSize) + local component index.
        forAll(vf, cellI)
        {
            for (label compI = 0; compI < nCompToInsert; ++compI)
            {
                const label xIndex = offset + cellI*xBlockSize + compI;
                const label vfIndex  = cellI*vfBlockSize + compList[compI];
                vfPtr[vfIndex] = x[xIndex];
            }
        }
    }
    else if (location_ == solutionLocation::POINTS)
    {
        notImplemented
        (
            "Function not implemented for solutionLocation::POINTS.\n"
            "Instead, use the alternative function:\n"
            "template <class Type>\n"
            "void foamPetscSnesHelper::ExtractFieldComponents\n"
            "(\n"
            "    const Vec x,\n"
            "    Field<Type>& vf,\n"
            "    const label offset,\n"
            "    const labelList& compIndices\n"
            ") const\n"
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = " << location_ << exit(FatalError);
    }
}


template <class Type>
void foamPetscSnesHelper::ExtractFieldComponents
(
    const Vec x,
    Field<Type>& vf,
    const label offset,
    const labelList& compIndices
) const
{
    // Access the x data
    const PetscScalar *xx;
    VecGetArrayRead(x, &xx);

    // Lookup x block size
    label xBlockSize;
    VecGetBlockSize(x, &xBlockSize);

    // Insert vf into xx
    if (location_ == solutionLocation::CELLS)
    {
        ExtractFieldComponents(xx, vf, offset, xBlockSize, compIndices);
    }
    else if (location_ == solutionLocation::POINTS)
    {
        // Obtain a pointer to the underlying scalar array in vf
        // We assume vf is stored as contiguous data, i.e.
        // (vx1, vy1, vz1, vx2, vy2, vz2, ..., vxN, vyN, vzN)
        scalar* vfPtr = reinterpret_cast<scalar*>(vf.begin());

        // Get the number of components per field element from the traits (e.g.,
        // 1 for scalar, 3 for vector, etc.)
        const label vfBlockSize = pTraits<Type>::nComponents;

        // Decide which components to extract: if compIndices is empty, default to
        // all components
        labelList compList;
        if (compIndices.size() == 0)
        {
            compList.setSize(vfBlockSize);
            for (label i = 0; i < vfBlockSize; ++i)
            {
                compList[i] = i;  // Default: copy all components in order
            }
        }
        else
        {
            compList = compIndices;
        }

        // The number of components to extract per block is determined by the size of
        // compList
        const label nCompToExtract = compList.size();
        // Extract the local processor values
        // Loop over each point in vf, which is owned by this processor, and
        // copy the selected components into the appropriate positions in x
        // We will also count the number of non-local points
        label localPointI = 0;
        label nNotOwnedByThisProc = 0;
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();
        forAll(vf, pointI)
        {
            if (ownedByThisProc[pointI])
            {
                for (label compI = 0; compI < nCompToExtract; ++compI)
                {
                    const label xIndex =
                        offset + localPointI*xBlockSize + compI;
                    const label vfIndex  = pointI*vfBlockSize + compList[compI];
                    vfPtr[vfIndex] = xx[xIndex];
                }

                localPointI++;
            }
            else
            {
                nNotOwnedByThisProc++;
            }
        }

        // Extract values not owned by this processor; we need to sync these from
        // other processors
        if (Pstream::parRun())
        {
            // Convert to the number of scalar unknowns on this processor but
            // are not owned by this processor
            nNotOwnedByThisProc *= xBlockSize;

            // Collect the global solution scalar indices of the non-local point
            // values

            // Create a list of global indices on the on-proc points not owned by
            // this proc
            List<PetscInt> indices(nNotOwnedByThisProc);
            label index = 0;
            const labelList& localToGlobalPointMap =
                globalPoints().localToGlobalPointMap();
            forAll(ownedByThisProc, pointI)
            {
                if (!ownedByThisProc[pointI])
                {
                    // Calculate the global point index for this point
                    const label globalPointI = localToGlobalPointMap[pointI];
                    for (label compI = 0; compI < nCompToExtract; ++compI)
                    {
                        indices[index++] =
                            offset + globalPointI*xBlockSize + compI;
                    }
                }
            }

            // Note: there are several PETSc objects below that could be stored
            // The challenge is that this function could be called for different
            // OpenFOAM fields (e.g., displacement, pressure) and the resulting
            // scalar indexing sets will be different (although they are
            // obviously related in terms of block rows). For now, leave it, but
            // it is possible to make this more efficient.

            // Create the index set
            // Could be stored
            IS indexSet;
            ISCreateGeneral
            (
                PETSC_COMM_WORLD,
                nNotOwnedByThisProc,
                indices.data(),
                PETSC_COPY_VALUES,
                &indexSet
            );

            // Create a local vector for holding non-local values
            // Could be stored
            Vec xNotOwned;
            VecCreate(PETSC_COMM_WORLD, &xNotOwned);
            VecSetSizes(xNotOwned, nNotOwnedByThisProc, PETSC_DECIDE);
            VecSetType(xNotOwned, VECMPI);
            VecSetUp(xNotOwned);

            // Create a context, which is required for syncing data
            // Could be stored
            VecScatter ctx;
            VecScatterCreate(x, indexSet, xNotOwned, NULL, &ctx);
            VecScatterBegin(ctx, x, xNotOwned, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(ctx, x, xNotOwned, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterDestroy(&ctx);

            // Populate not-owned values array
            // Could be stored
            PetscScalar* xNotOwnedArr;
            VecGetArray(xNotOwned, &xNotOwnedArr);
            {
                label index = 0;
                forAll(vf, pointI)
                {
                    if (!ownedByThisProc[pointI])
                    {
                        for (label compI = 0; compI < nCompToExtract; ++compI)
                        {
                            const label vfIndex  =
                                pointI*vfBlockSize + compList[compI];
                            vfPtr[vfIndex] = xNotOwnedArr[index++];
                        }
                    }
                }
            }
            VecRestoreArray(xNotOwned, &xNotOwnedArr);

            // Destroy xNotOwned
            VecDestroy(&xNotOwned);

            // Destroy the index set
            ISDestroy(&indexSet);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = " << location_ << exit(FatalError);
    }

    // Restore the x vector
    VecRestoreArrayRead(x, &xx);
}


template <class Type>
void foamPetscSnesHelper::InsertFieldComponents
(
    const Field<Type>& vf,
    PetscScalar *x,
    const label offset,
    const label xBlockSize,
    const labelList& compIndices
) const
{
    // Get a pointer to the underlying scalar data of vf (assumes contiguous
    // storage)
    const scalar* vfPtr = reinterpret_cast<const scalar*>(vf.begin());

    // Determine the number of components per field element (e.g., 1 for
    // scalar, 3 for vector, etc.)
    const label vfBlockSize = pTraits<Type>::nComponents;

    // Decide which components to extract: if compIndices is empty, default to
    // all components
    labelList compList;
    if (compIndices.size() == 0)
    {
        compList.setSize(vfBlockSize);
        for (label i = 0; i < vfBlockSize; ++i)
        {
            compList[i] = i;  // Default: copy all components in order
        }
    }
    else
    {
        compList = compIndices;
    }

    // The number of components to insert per cell is determined by the size of
    // compList
    const label nCompToInsert = compList.size();

    if (location_ == solutionLocation::CELLS)
    {
        // Loop over each cell in vf and copy the selected components into the
        // appropriate positions in x.
        // The destination index for cell 'cellI' is computed as:
        //   offset + (cellI*xBlockSize) + local component index.
        forAll(vf, cellI)
        {
            for (label compI = 0; compI < nCompToInsert; ++compI)
            {
                const label xIndex = offset + cellI*xBlockSize + compI;
                const label vfIndex  = cellI*vfBlockSize + compList[compI];
                x[xIndex] = vfPtr[vfIndex];
            }
        }
    }
    else if (location_ == solutionLocation::POINTS)
    {
        // We assume here that the vf field has already been synced in parallel
        // i.e., that the points on processor boundaries are already updated and
        // consistent
        // We only insert the values owned by this processor
        label localPointI = 0;
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();
        forAll(vf, pointI)
        {
            if (ownedByThisProc[pointI])
            {
                for (label compI = 0; compI < nCompToInsert; ++compI)
                {
                    const label xIndex =
                        offset + localPointI*xBlockSize + compI;
                    const label vfIndex  = pointI*vfBlockSize + compList[compI];
                    x[xIndex] = vfPtr[vfIndex];
                }

                localPointI++;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = " << location_ << exit(FatalError);
    }
}


template <class Type>
void foamPetscSnesHelper::InsertFieldComponents
(
    const Field<Type>& vf,
    Vec x,
    const label offset,
    const labelList& compIndices
) const
{
    // Access the x data
    PetscScalar *xx;
    VecGetArray(x, &xx);

    // Lookup x block size
    label xBlockSize;
    VecGetBlockSize(x, &xBlockSize);

    // Insert vf into xx
    InsertFieldComponents(vf, xx, offset, xBlockSize, compIndices);

    // Restore the x vector
    VecRestoreArray(x, &xx);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
