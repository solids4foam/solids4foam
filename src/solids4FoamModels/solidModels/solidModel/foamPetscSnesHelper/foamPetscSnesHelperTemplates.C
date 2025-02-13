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
    const fvMatrix<Type>& fvM,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns
) const
{
    // Initialise block coeff
    const label nCoeffCmpts = blockSize_*blockSize_;
    PetscScalar values[nCoeffCmpts];
    std::memset(values, 0, sizeof(values));

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
                //reinterpret_cast<const scalar*>(diag[blockRowI].begin());
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
                            (cmptI + rowOffset)*blockSize_ + cmptJ + colOffset
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
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                    ADD_VALUES
                );
            );
        }
    }

    // Insert the off-diagonal
    if (fvM.hasUpper())
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
                            (cmptI + rowOffset)*blockSize_ + cmptJ + colOffset
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

            // Insert block coefficients
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI, values,
                    ADD_VALUES
                );
            );
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI, values,
                    ADD_VALUES
                );
            );
        }
    }

    // Collect the global cell indices from neighbours at processor boundaries
    // These are used to insert the off-processor coefficients
    // First, send the data
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == "processor")
        {
            // Take a copy of the faceCells (local IDs) and convert them to
            // global IDs
            labelList globalFaceCells(mesh.boundary()[patchI].faceCells());
            foamPetscSnesHelper::globalCells().inplaceToGlobal(globalFaceCells);

            // Send global IDs to the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh.boundary()[patchI]);
            procPatch.send
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );
        }
    }

    // Next, receive the data
    PtrList<labelList> neiProcGlobalIDs(mesh.boundaryMesh().size());
    forAll(mesh.boundaryMesh(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];
        if (fp.type() == "processor")
        {
            neiProcGlobalIDs.set(patchI, new labelList(fp.size()));
            labelList& globalFaceCells = neiProcGlobalIDs[patchI];

            // Receive global IDs from the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh.boundary()[patchI]);
            procPatch.receive
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );
        }
    }

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
                            // intCoeffs[faceI].begin()
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
                                    (cmptI + rowOffset)*blockSize_ + cmptJ + colOffset
                                ] = curIntCoeffsPtr[cmptI];
                            }
                        }
                    }

                    CHKERRQ
                    (
                        MatSetValuesBlocked
                        (
                            jac, 1, &globalBlockRowI, 1, &globalBlockRowI, values,
                            ADD_VALUES
                        );
                    );
                }

                // Off-proc off-diagonal coefficient
                {
                    // Obtain a pointer to the underlying scalar data in
                    // neiCoeffs[faceI] (assumes contiguous storage)
                    const scalar* curNeiCoeffsPtr =
                        reinterpret_cast<const scalar*>
                        (
                            // neiCoeffs[faceI].begin()
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
                                    (cmptI + rowOffset)*blockSize_ + cmptJ + colOffset
                                ] = -curNeiCoeffsPtr[cmptI];
                            }
                        }
                    }

                    const label globalBlockColI = neiGlobalFaceCells[faceI];

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
    const label nScalarEqns
) const
{
    // Obtain a pointer to the underlying scalar array in vf
    // We assume vf is stored as contiguous data, i.e.
    // (vx1, vy1, vz1, vx2, vy2, vz2, ..., vxN, vyN, vzN)
    scalar* vfPtr = reinterpret_cast<scalar*>(vf.begin());

    // Get the number of components per field element from the traits (e.g.,
    // 1 for scalar, 3 for vector, etc.)
    const label vfBlockSize = pTraits<Type>::nComponents;

    // Loop over each element in vf and copy the appropriate components from x
    forAll(vf, blockI)
    {
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            vfPtr[blockI*vfBlockSize + cmptI] = x[offset + blockI*blockSize_ + cmptI];
        }
    }
}


template <class Type>
void foamPetscSnesHelper::ExtractFieldComponents
(
    const Vec x,
    Field<Type>& vf,
    const label offset,
    const label nScalarEqns
) const
{
    // Access the x data
    const PetscScalar *xx;
    VecGetArrayRead(x, &xx);

    // Insert vf into xx
    ExtractFieldComponents(xx, vf, offset, nScalarEqns);

    // Restore the x vector
    VecRestoreArrayRead(x, &xx);

}


template <class Type>
void foamPetscSnesHelper::InsertFieldComponents
(
    const Field<Type>& vf,
    PetscScalar *x,
    const label offset,
    const label nScalarEqns
) const
{
    // Obtain a pointer to the underlying scalar data in vf (assumes contiguous
    // storage)
    const scalar* vfPtr = reinterpret_cast<const scalar*>(vf.begin());

    // Get the number of components per field element from the traits (e.g.,
    // 1 for scalar, 3 for vector, etc.)
    const label vfBlockSize = pTraits<Type>::nComponents;

    // Loop over each element in vf and copy its components into the corresponding position in x
    forAll(vf, blockI)
    {
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            x[offset + blockI*blockSize_ + cmptI] = vfPtr[blockI*vfBlockSize + cmptI];
        }
    }
}


template <class Type>
void foamPetscSnesHelper::InsertFieldComponents
(
    const Field<Type>& vf,
    Vec x,
    const label offset,
    const label nScalarEqns
) const
{
    // Access the x data
    PetscScalar *xx;
    VecGetArray(x, &xx);

    // Insert vf into xx
    InsertFieldComponents(vf, xx, offset, nScalarEqns);

    // Restore the x vector
    VecRestoreArray(x, &xx);

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
