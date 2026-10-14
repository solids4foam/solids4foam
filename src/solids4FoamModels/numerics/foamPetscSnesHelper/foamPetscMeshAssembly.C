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

#include "foamPetscMeshAssembly.H"
#include "processorFvPatch.H"
#include "symmetryFvPatchFields.H"
#include "leastSquaresVectors.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "fvm.H"
#include "petscErrorHandling.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "symmetryPlaneFvPatchFields.H"
#endif

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(foamPetscMeshAssembly, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void foamPetscMeshAssembly::checkMesh(const fvMesh& mesh) const
{
    if (&mesh != &mesh_)
    {
        FatalErrorInFunction
            << "The field is on mesh " << mesh.name() << " but this assembly "
            << "is for mesh " << mesh_.name() << abort(FatalError);
    }
}


void foamPetscMeshAssembly::makeNeiProcFields(const fvMesh& mesh) const
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
#ifdef OPENFOAM_COM
            foamPetscMeshAssembly::globalCells().inplaceToGlobal(globalFaceCells);
#else
            forAll(globalFaceCells, cI)
            {
                const label localCellID = globalFaceCells[cI];
                globalFaceCells[cI] =
                    foamPetscMeshAssembly::globalCells().toGlobal(localCellID);
            }
#endif

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


const PtrList<labelList>& foamPetscMeshAssembly::neiProcGlobalIDs
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


const PtrList<scalarField>& foamPetscMeshAssembly::neiProcVolumes
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


const leastSquaresS4fVectors& foamPetscMeshAssembly::lsVectors
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

#ifdef OPENFOAM_COM
    return leastSquaresS4fVectors::New(lsName, p.mesh(), useBoundaryFaceValues);
#else
    return leastSquaresS4fVectors::New(p.mesh(), useBoundaryFaceValues);
#endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamPetscMeshAssembly::foamPetscMeshAssembly
(
    const fvMesh& mesh,
    const solutionLocation location
)
:
    mesh_(mesh),
    location_(location),
    globalCellsPtr_
    (
        location == solutionLocation::CELLS
      ? new globalIndex(mesh.nCells())
      : nullptr
    ),
    globalPointsPtr_
    (
        location == solutionLocation::POINTS
      ? new globalPointIndices(mesh)
      : nullptr
    ),
    neiProcGlobalIDs_(),
    neiProcVolumes_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label foamPetscMeshAssembly::initialiseJacobian
(
    Mat& jac,
    const label blockSize,
    const bool createMat
) const
{
    const fvMesh& mesh = mesh_;

    // Number of local block equations
    label blockn = -1;

    // Number of local scalar equations, i.e. blockn*blockSize
    label n = -1;

    // Number of scalar equations across all processors
    label N = -1;

    if (location_ == solutionLocation::CELLS)
    {
        // Set the number local block equations
        blockn = mesh.nCells();

        // Set the number local scalar equations
        n = mesh.nCells()*blockSize;

        // Set the global system size
        N = returnReduce(n, sumOp<label>());
    }
    else if (location_ == solutionLocation::POINTS)
    {
        // Note: the size of x is, in general, not equal to the number of
        // points on this processor, as x only contains the points owned by
        // this processor. To access the values not-owned by the proc, we
        // need to request the values from the other processors

        // Take references for brevity and efficiency
        const labelList& localToGlobalPointMap =
            globalPoints().localToGlobalPointMap();
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();

        // Find size of global system
        // i.e. the highest global point index + 1
        const label blockN = gMax(localToGlobalPointMap) + 1;
        N = blockSize*blockN;

        // Find the start and end global point indices for this proc
        label blockStartID = N;
        label blockEndID = -1;
        forAll(ownedByThisProc, pI)
        {
            if (ownedByThisProc[pI])
            {
                blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
                blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
            }
        }

        // Find size of local system, i.e. the range of points owned by this
        // proc
        blockn = blockEndID - blockStartID + 1;
        n = blockSize*blockn;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = "
            << solutionLocationNames_[location_]
            << exit(FatalError);
    }

    // Set the Jacobian matrix size
    if (createMat)
    {
        MatCreate(PETSC_COMM_WORLD, &jac);
        MatSetFromOptions(jac);
        MatSetSizes(jac, n, n, N, N);

        MatSetType(jac, MATMPIAIJ);
    }

    // Set the block size
    AssertPETSc(MatSetBlockSize(jac, blockSize));

    // Count the number of non-zeros in the matrix
    // Note: we assume a compact stencil, i.e. face only face neighbours

    // Number of on-processor non-zeros per row
    // Initialise d_nnz to one
    labelList d_nnz(blockn, 1);

    // Number of off-processor non-zeros per row
    // Initialise o_nnz to zero
    labelList o_nnz(blockn, 0);

    // Count the neighbours
    if (location_ == solutionLocation::CELLS)
    {
        List<labelHashSet> gradCols(blockn);
        List<labelHashSet> rowCols(blockn);

        forAll(rowCols, cellI)
        {
            gradCols[cellI] = labelHashSet(16);
            rowCols[cellI] = labelHashSet(32);

            const label globalCellID = globalCells().toGlobal(cellI);
            gradCols[cellI].insert(globalCellID);
            rowCols[cellI].insert(globalCellID);
        }

        const Foam::labelUList& own = mesh.owner();
        const Foam::labelUList& nei = mesh.neighbour();
        forAll(own, faceI)
        {
            const Foam::label ownCellID = own[faceI];
            const Foam::label neiCellID = nei[faceI];

            const label globalOwnCellID = globalCells().toGlobal(ownCellID);
            const label globalNeiCellID = globalCells().toGlobal(neiCellID);

            gradCols[ownCellID].insert(globalNeiCellID);
            gradCols[neiCellID].insert(globalOwnCellID);

            rowCols[ownCellID].insert(globalNeiCellID);
            rowCols[neiCellID].insert(globalOwnCellID);
        }

        const PtrList<labelList>& neiProcGlobalIDs =
            this->neiProcGlobalIDs(mesh);

        forAll(mesh.boundary(), patchI)
        {
            if (mesh.boundary()[patchI].type() == "processor")
            {
                const Foam::labelUList& faceCells =
                    mesh.boundary()[patchI].faceCells();
                const labelList& neiGlobalFaceCells =
                    neiProcGlobalIDs[patchI];

                forAll(faceCells, fcI)
                {
                    const Foam::label cellID = faceCells[fcI];
                    gradCols[cellID].insert(neiGlobalFaceCells[fcI]);
                    rowCols[cellID].insert(neiGlobalFaceCells[fcI]);
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

        auto insertSet =
            [](labelHashSet& row, const labelHashSet& cols)
            {
                forAllConstIter(labelHashSet, cols, iter)
                {
                    row.insert(iter.key());
                }
            };

        // The pressure stabilisation may contain div(interpolate(grad(p))).
        // This reaches the gradient stencil of the neighbouring cell across
        // each internal face, not just the direct face-neighbour stencil.
        forAll(own, faceI)
        {
            const Foam::label ownCellID = own[faceI];
            const Foam::label neiCellID = nei[faceI];

            insertSet(rowCols[ownCellID], gradCols[neiCellID]);
            insertSet(rowCols[neiCellID], gradCols[ownCellID]);
        }

#if defined(OPENFOAM_ORG) || defined(FOAMEXTEND)
        const label myProcNo = Pstream::myProcNo();
        const label gStart = globalCells().offset(myProcNo);
        const label gEnd = gStart + globalCells().localSize(myProcNo);
#else
        const label gStart = globalCells().localStart();
        const label gEnd = globalCells().localEnd();
#endif

        d_nnz = 0;

        forAll(rowCols, cellI)
        {
            forAllConstIter(labelHashSet, rowCols[cellI], iter)
            {
                const label gCol = iter.key();

                if (gCol >= gStart && gCol < gEnd)
                {
                    ++d_nnz[cellI];
                }
                else
                {
                    ++o_nnz[cellI];
                }
            }
        }
    }
    else if (location_ == solutionLocation::POINTS)
    {
        // Take references
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();
        const labelList& stencilSizeOwned = globalPoints().stencilSizeOwned();
        const labelList& stencilSizeNotOwned =
            globalPoints().stencilSizeNotOwned();

        // Reset d_nnz to 0
        d_nnz = 0;

        // We can optionally track the max on-core non-zeros to zero
        //label d_nz = 0;

        // Count non-zeros
        label rowI = 0;
        forAll(stencilSizeOwned, blockRowI)
        {
            if (ownedByThisProc[blockRowI])
            {
                const label nCompOwned = stencilSizeOwned[blockRowI];
                const label nCompNotOwned = stencilSizeNotOwned[blockRowI];

                d_nnz[rowI] += nCompOwned;
                o_nnz[rowI++] += nCompNotOwned;

                //d_nz = max(d_nz, nCompOwned);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown solution location = " << location_
            << exit(FatalError);
    }

    // Allocate parallel matrix
    // Note: MatXAIJSetPreallocation is used rather than a type-specific
    // routine, e.g. MatMPIBAIJSetPreallocation, because the type-specific
    // routines dispatch through PetscTryMethod and so silently do nothing when
    // the matrix is of another type, leaving the matrix un-preallocated:
    // PETSc-3.19 and newer then quietly call MatSetUp, whereas older versions
    // give a "wrong state" error in debug builds and segfault in optimised
    // builds. MatXAIJSetPreallocation takes the same block-row counts and
    // dispatches to whichever type the matrix actually has, so the matrix is
    // preallocated whether it is AIJ, BAIJ or SBAIJ, e.g. if the user
    // overrides the type with -mat_type
    AssertPETSc
    (
        MatXAIJSetPreallocation
        (
            jac, blockSize, d_nnz.data(), o_nnz.data(), NULL, NULL
        )
    );

    // Raise an error if mallocs are required during matrix assembly
    MatSetOption(jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    return 0;
}


label foamPetscMeshAssembly::initialiseSolution
(
    Vec& x,
    const label blockSize,
    const bool createVec
) const
{
    const fvMesh& mesh = mesh_;

    if (createVec)
    {
        // Number of local scalar equations
        label n = -1;

        // Number of scalar equations across all processors
        label N = -1;

        if (location_ == solutionLocation::CELLS)
        {
            // Set the number local scalar equations
            n = mesh.nCells()*blockSize;

            // Set the global system size
            N = returnReduce(n, sumOp<label>());
        }
        else if (location_ == solutionLocation::POINTS)
        {
            // Note: the size of x is, in general, not equal to the number of
            // points on this processor, as x only contains the points owned by
            // this processor. To access the values not-owned by the proc, we
            // need to request the values from the other processors

            // Take references for brevity and efficiency
            const labelList& localToGlobalPointMap =
                globalPoints().localToGlobalPointMap();
            const boolList& ownedByThisProc = globalPoints().ownedByThisProc();

            // Find size of global system
            // i.e. the highest global point index + 1
            const label blockN = gMax(localToGlobalPointMap) + 1;
            N = blockSize*blockN;

            // Find the start and end global point indices for this proc
            label blockStartID = N;
            label blockEndID = -1;
            forAll(ownedByThisProc, pI)
            {
                if (ownedByThisProc[pI])
                {
                    blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
                    blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
                }
            }

            //const label startID = blockSize*blockStartID;
            //const label endID = blockSize*(blockEndID + 1) - 1;

            // Find size of local system, i.e. the range of points owned by this
            // proc
            const label blockn = blockEndID - blockStartID + 1;
            n = blockSize*blockn;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown solution location = " << location_
                << exit(FatalError);
        }

        x = Vec();
        AssertPETSc(VecCreate(PETSC_COMM_WORLD, &x));
        AssertPETSc(VecSetSizes(x, n, N));
        AssertPETSc(VecSetType(x, VECMPI));
    }

    AssertPETSc(VecSetBlockSize(x, blockSize));
    AssertPETSc(PetscObjectSetName((PetscObject) x, "Solution"));
    AssertPETSc(VecSetFromOptions(x));
    AssertPETSc(VecZeroEntries(x));

    return 0;
}


label foamPetscMeshAssembly::InsertFvmGradIntoPETScMatrix
(
    const volScalarField& p,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const scalar scale
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

    // const scalar sign = flipSign ? -1.0 : 1.0;
    const scalar sign = scale;

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise the block coefficient
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    forAll(own, faceI)
    {
        // Local block row ID
        const label ownCellID = own[faceI];

        // Local block column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            foamPetscMeshAssembly::globalCells().toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            foamPetscMeshAssembly::globalCells().toGlobal(neiCellID);

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
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                values.cdata(),
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
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockRowI, 1, &globalBlockColI, values.cdata(),
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
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockColI, 1, &globalBlockColI, values.cdata(),
                ADD_VALUES
            )
        );

        // lsGrad[neiCellID] -= neiLs[faceI]*(vsf[neiCellID] - vsf[ownCellID]);
        // mat(neIFaceI, ownCellID) -= VI[neiCellID]*neiLs[faceI];
        // Flip the sign of the values
        for (label i = 0; i < nCoeffCmpts; ++i)
        {
            values[i] = -values[i];
        }
        AssertPETSc
        (
            MatSetValuesBlocked
            (
                jac, 1, &globalBlockColI, 1, &globalBlockRowI, values.cdata(),
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

            forAll(fp, patchFaceI)
            {
                // Local block row ID
                const label ownCellID = faceCells[patchFaceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscMeshAssembly::globalCells().toGlobal(ownCellID);

                // On-proc diagonal coefficient
                // mat(ownCellID, ownCellID) -= VI[ownCellID]*ownLs[faceI];
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize + colOffset] =
                        -sign*VI[ownCellID]*patchOwnLs[patchFaceI][cmptI];
                }
                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );

                // Neighbour global cell ID
                const label globalBlockColI = neiGlobalFaceCells[patchFaceI];

                // Off-proc off-diagonal coefficient
                // mat(ownCellID, neiCellID) += VI[ownCellID]*ownLs[faceI];
                for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                {
                    values[cmptI*blockSize + colOffset] =
                        sign*VI[ownCellID]
                       *patchOwnLs[patchFaceI][cmptI];
                }
                AssertPETSc
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
        else if (fp.coupled()) // coupled but not processor
        {
            FatalErrorInFunction
                << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        else if
        (
            isA<symmetryPolyPatch>(fp.patch())
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(fp.patch())
#endif
        )
        {
            // The delta in scalar p across the symmetry is zero by definition
            // so the symmetry plane does not contribution coefficients
        }
        else
        {
            if (useBoundaryFaceValues[patchI])
            {
                const fvPatchScalarField& pp = p.boundaryField()[patchI];
                const scalarField intCoeffs(pp.valueInternalCoeffs(fp.weights()));

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
                        foamPetscMeshAssembly::globalCells().toGlobal(ownCellID);

                    // mat(ownCellID, ownCellID) -=
                    //     VI[ownCellID]*patchOwnLs[patchFaceI];
                    for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
                    {
                        values[cmptI*blockSize + colOffset] =
                            sign*VI[ownCellID]
                           *(intCoeffs[patchFaceI] - 1.0)
                           *patchOwnLs[patchFaceI][cmptI];
                    }
                    AssertPETSc
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
        }
    }

    return 0;
}


label foamPetscMeshAssembly::InsertFvcDivGradInterpolateIntoPETScMatrix
(
    const volScalarField& p,
    const surfaceScalarField& gamma,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const scalar scale
) const
{
    const fvMesh& mesh = p.mesh();
    const leastSquaresS4fVectors& lsv = lsVectors(p);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();
    const boolList& useBoundaryFaceValues = lsv.useBoundaryFaceValues();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    List<DynamicList<label>> gradCols(mesh.nCells());
    List<DynamicList<vector>> gradCoeffs(mesh.nCells());

    auto appendGradCoeff =
        [&](const label cellI, const label globalColI, const vector& coeff)
        {
            gradCols[cellI].append(globalColI);
            gradCoeffs[cellI].append(coeff);
        };

    forAll(own, faceI)
    {
        const label ownCellID = own[faceI];
        const label neiCellID = nei[faceI];

        const label globalOwnCellID = globalCells().toGlobal(ownCellID);
        const label globalNeiCellID = globalCells().toGlobal(neiCellID);

        appendGradCoeff(ownCellID, globalOwnCellID, -ownLs[faceI]);
        appendGradCoeff(ownCellID, globalNeiCellID, ownLs[faceI]);

        appendGradCoeff(neiCellID, globalNeiCellID, -neiLs[faceI]);
        appendGradCoeff(neiCellID, globalOwnCellID, neiLs[faceI]);
    }

    const PtrList<labelList>& neiProcGlobalIDs = this->neiProcGlobalIDs(mesh);

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];
        const labelUList& faceCells = fp.faceCells();
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchI];

        if (fp.type() == "processor")
        {
            const labelList& neiGlobalFaceCells = neiProcGlobalIDs[patchI];

            forAll(fp, patchFaceI)
            {
                const label ownCellID = faceCells[patchFaceI];
                const label globalOwnCellID = globalCells().toGlobal(ownCellID);

                appendGradCoeff
                (
                    ownCellID, globalOwnCellID, -patchOwnLs[patchFaceI]
                );
                appendGradCoeff
                (
                    ownCellID,
                    neiGlobalFaceCells[patchFaceI],
                    patchOwnLs[patchFaceI]
                );
            }
        }
        else if (fp.coupled())
        {
            // Coupled non-processor patches are not handled elsewhere in the
            // PETSc helper either.
            FatalErrorInFunction
                << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        else if
        (
            isA<symmetryPolyPatch>(fp.patch())
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(fp.patch())
#endif
        )
        {
            // The scalar jump across a symmetry plane is zero.
        }
        else if (useBoundaryFaceValues[patchI])
        {
            const fvPatchScalarField& pp = p.boundaryField()[patchI];
            const scalarField intCoeffs(pp.valueInternalCoeffs(fp.weights()));

            forAll(faceCells, patchFaceI)
            {
                const label ownCellID = faceCells[patchFaceI];
                const label globalOwnCellID = globalCells().toGlobal(ownCellID);

                appendGradCoeff
                (
                    ownCellID,
                    globalOwnCellID,
                    (intCoeffs[patchFaceI] - 1.0)*patchOwnLs[patchFaceI]
                );
            }
        }
    }

    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);
    const label valueOffset = rowOffset*blockSize + colOffset;

    auto insertGradDivStencil =
        [&]
        (
            const label rowCellID,
            const vector& faceCoeff,
            const DynamicList<label>& cols,
            const DynamicList<vector>& coeffs
        )
        {
            const label globalBlockRowI = globalCells().toGlobal(rowCellID);

            forAll(cols, coeffI)
            {
                const label globalBlockColI = cols[coeffI];

                // Insert unconditionally to establish the full sparsity
                // pattern on the first assembly.  PETSc compresses the
                // storage after MatAssemblyEnd, so skipping zero entries
                // would lose their positions; on a moving mesh the
                // coefficients may later become nonzero, triggering a
                // MAT_NEW_NONZERO_ALLOCATION_ERR.
                values[valueOffset] = faceCoeff & coeffs[coeffI];
                AssertPETSc
                (
                    MatSetValuesBlocked
                    (
                        jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                        values.cdata(),
                        ADD_VALUES
                    )
                );
                values[valueOffset] = 0.0;
            }
        };

    const scalarField& gammaI = gamma;
    const vectorField& SfI = mesh.Sf();
    const scalarField& wI = mesh.weights();

    forAll(own, faceI)
    {
        const label ownCellID = own[faceI];
        const label neiCellID = nei[faceI];

        const vector ownFaceCoeff =
            scale*gammaI[faceI]*wI[faceI]*SfI[faceI];
        const vector neiFaceCoeff =
            scale*gammaI[faceI]*(1.0 - wI[faceI])*SfI[faceI];

        insertGradDivStencil
        (
            ownCellID, ownFaceCoeff,
            gradCols[ownCellID], gradCoeffs[ownCellID]
        );
        insertGradDivStencil
        (
            ownCellID, neiFaceCoeff,
            gradCols[neiCellID], gradCoeffs[neiCellID]
        );

        insertGradDivStencil
        (
            neiCellID, -ownFaceCoeff,
            gradCols[ownCellID], gradCoeffs[ownCellID]
        );
        insertGradDivStencil
        (
            neiCellID, -neiFaceCoeff,
            gradCols[neiCellID], gradCoeffs[neiCellID]
        );
    }

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];
        const labelUList& faceCells = fp.faceCells();
        const vectorField& pSf = fp.Sf();
        const fvsPatchScalarField& pGamma = gamma.boundaryField()[patchI];

        if (fp.type() == "empty")
        {
            continue;
        }
        else if (fp.coupled())
        {
            const fvsPatchScalarField& pw =
                mesh.weights().boundaryField()[patchI];

            forAll(fp, patchFaceI)
            {
                const label ownCellID = faceCells[patchFaceI];
                const vector faceCoeff =
                    scale*pGamma[patchFaceI]*pw[patchFaceI]*pSf[patchFaceI];

                insertGradDivStencil
                (
                    ownCellID, faceCoeff,
                    gradCols[ownCellID], gradCoeffs[ownCellID]
                );
            }
        }
        else
        {
            forAll(fp, patchFaceI)
            {
                const label ownCellID = faceCells[patchFaceI];
                const vector faceCoeff =
                    scale*pGamma[patchFaceI]*pSf[patchFaceI];

                insertGradDivStencil
                (
                    ownCellID, faceCoeff,
                    gradCols[ownCellID], gradCoeffs[ownCellID]
                );
            }
        }
    }

    return 0;
}


label foamPetscMeshAssembly::InsertFvmDivPhiUIntoPETScMatrix
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const scalar scale
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
    // const scalar sign = flipSign ? -1 : 1;
    const scalar sign = scale;
    tensor coeff = tensor::zero;

    // The first-order upwind discretisation for cell P is
    //     div(phi,U)_P = \sum phi*Uf
    // where
    //     Uf = Up     if phi > 0
    //     Uf = Un     otherwise
    // and
    //     phi = Sf & [w*Up + (1 - w)*Un]
    //
    // The Newton linearisation of div(phi,U) for cell P is
    //     \sum (d/dUp)(phi*Uf) = \sum phi*(dUf/dUp) + Uf*(dphi/dUp)
    //
    // Using the following
    //     dUp/dUp = I
    //     dUn/dUn = I
    //     dphi/dUp = w*Sf
    //     dphi/dUn = (1 - w)*Sf
    // we can calculate the block coefficient for a face:
    // if phi > 0
    //     Uf = Up
    //     (d/dUp)(phi*Up) = phi*I + w*Up*Sf        => diagonal
    //     (d/dUn)(phi*Up) = (1 - w)*Up*Sf          => off-diagonal
    // else
    //     Uf = Un
    //     (d/dUp)(phi*Un) = w*Un*Sf                => diagonal
    //     (d/dUn)(phi*Un) = phi*I + (1 - w)*Un*Sf  => off-diagonal
    //
    // The w*Sf*Up term is missing from fvm::div(phi, U) as it requires a
    // coupled solver. The caller inserts the segregated fvm::div(phi,U)
    // matrix; this routine adds only the missing tensor terms.
    // From the neighbours perspective, the coefficient is (1 - w)*Sf*Un

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Prepare coeff array
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    forAll(phiI, faceI)
    {
        // Local row ID
        const label ownCellID = own[faceI];

        // Local column ID
        const label neiCellID = nei[faceI];

        // Global block row ID
        const label globalBlockRowI =
            foamPetscMeshAssembly::globalCells().toGlobal(ownCellID);

        // Global block column ID
        const label globalBlockColI =
            foamPetscMeshAssembly::globalCells().toGlobal(neiCellID);

        // Take references for brevity, clarity and efficiency
        const scalar w = wI[faceI];
        const vector& Sf = SfI[faceI];
        const vector& Up = UI[ownCellID];
        const vector& Un = UI[neiCellID];

        if (phiI[faceI] > 0)
        {
            // Add w*Up*Sf to owner eqn
            coeff = sign*w*Up*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add (1 - w)*Up*Sf as nei to contribution to own eqn
            coeff = sign*(1.0 - w)*Up*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Up*Sf to the neighbour diagonal
            // Flip the sign as Sf is opposite
            coeff = -sign*(1.0 - w)*Up*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -w*Up*Sf as the own contribution to the nei eqn
            // Flip the sign as Sf is opposite
            coeff = -sign*w*Up*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );
        }
        else
        {
            // Add w*Un*Sf to owner diagonal
            coeff = sign*w*Un*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add (1 - w)*Un*Sf as nei to contribution to own eqn
            // coeff = sign*w*Un*Sf;
            coeff = sign*(1.0 - w)*Un*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockRowI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -(1 - w)*Un*Sf to neighbour diagonal
            // Flip the sign as Sf is opposite
            coeff = -sign*(1.0 - w)*Un*Sf;
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
            AssertPETSc
            (
                MatSetValuesBlocked
                (
                    jac, 1, &globalBlockColI, 1, &globalBlockColI,
                    values.cdata(),
                    ADD_VALUES
                )
            );

            // Add -w*Un*Sf as own to contribution to nei eqn
            // Flip the sign as Sf is opposite
            coeff = -sign*w*Un*Sf;
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
            AssertPETSc
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
            // TODO: add the extra Newton linearisation tensor terms
            // (w*Sf*U and (1-w)*Sf*U) for processor faces. Omitting
            // these only affects the Jacobian accuracy at processor
            // boundaries; the segregated fvm::div part is already
            // handled correctly by InsertFvMatrixIntoPETScMatrix.
        }
        else if
        (
            patch.type() != "empty"
         && !isA<symmetryFvPatchField<vector>>(pU)
#ifdef OPENFOAM_NOT_EXTEND
         && !isA<symmetryPlaneFvPatchField<vector>>(pU)
#endif
         && !pU.fixesValue()
        )
        {
            forAll(fc, faceI)
            {
                // Local row ID
                const label ownCellID = fc[faceI];

                // Global block row ID
                const label globalBlockRowI =
                    foamPetscMeshAssembly::globalCells().toGlobal(ownCellID);

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
                AssertPETSc
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
    }

    return 0;
}


label foamPetscMeshAssembly::InsertFvmDivUIntoPETScMatrix
(
    const volScalarField& p,
    const volVectorField& U,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const scalar scale
) const
{
    // Take references for efficiency and brevity
    const fvMesh& mesh = p.mesh();
    const scalar sign = scale;

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

        // Interior face coefficients (scaled by sign so the whole row
        // is scaled consistently). Note: these are the negation of
        // OpenFOAM's `fvm::div(phi, U)` coefficient convention, so the
        // assembled matrix for sign=+1 represents -V*div(U) per the
        // operator action at interior cells but with a self-consistent
        // boundary treatment that pairs with the gradient helper. In
        // practice the empirically self-consistent sign for the
        // displacement-pressure Newton system on this discretisation
        // requires the caller to pass scale = -alpha when assembling
        // J_pD = -alpha*V*div(D); see the call site in
        // linGeomTotalDispSolid::formJacobian for details.
        lower = sign*w*Sf.component(cmptI);
        upper = lower - sign*Sf.component(cmptI);

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

            if (patch.coupled())
            {
                if (patch.type() == "processor")
                {
                    // Convention: InsertFvMatrixIntoPETScMatrix inserts
                    //   +internalCoeffs to the diagonal
                    //   -boundaryCoeffs to the off-diagonal
                    //
                    // Interior face convention (lower = +sign*w*Sf,
                    // upper = sign*(w-1)*Sf, negSumDiag adds
                    // -sign*w*Sf to the owner diagonal per face).
                    // Reproduce the same stencil across processor
                    // boundaries by setting:
                    //   internalCoeffs[f] = -sign*w*Sf
                    //     (matches negSumDiag's contribution).
                    //   boundaryCoeffs[f] = sign*(1-w)*Sf
                    //     so -boundaryCoeffs = sign*(w-1)*Sf = upper[f].
                    scalarField& patchIntCoeffs =
                        divUCoeffs.internalCoeffs()[patchI];
                    scalarField& patchBndCoeffs =
                        divUCoeffs.boundaryCoeffs()[patchI];

                    forAll(pU, faceI)
                    {
                        patchIntCoeffs[faceI] =
                            -sign*pw[faceI]*Sf[faceI][cmptI];
                        patchBndCoeffs[faceI] =
                            sign*(1.0 - pw[faceI])*Sf[faceI][cmptI];
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Coupled boundary type " << patch.type()
                        << " not yet supported in "
                        << "InsertFvmDivUIntoPETScMatrix"
                        << abort(FatalError);
                }
            }
            else
            {
                // Non-coupled boundary contribution.
                // With the fvm::div-style interior coefficients above,
                // negSumDiag already yields the correct V*div diagonal
                // at Dirichlet boundary cells (where the BC's
                // valueInternalCoeffs is zero). For BCs where the face
                // value depends on the interior cell value
                // (e.g. zeroGradient, where valueInternalCoeffs = 1),
                // an additional `-sign*intCoeffs*Sf` correction per
                // face is required to recover the true V*div diagonal.
                const vectorField internalCoeffs
                (
                    pU.valueInternalCoeffs(pw)
                );

                forAll(pU, faceI)
                {
                    diag[fc[faceI]] -=
                        sign
                       *internalCoeffs[faceI][cmptI]
                       *Sf[faceI][cmptI];
                }
            }
        }

        // Insert component coeffs
        // cmptI is the 1st column for Ux, 2nd for Uy, 3rd for Uz
        foamPetscMeshAssembly::InsertFvMatrixIntoPETScMatrix<scalar>
        (
            divUCoeffs, jac, rowOffset, colOffset + cmptI, 1
        );
    }

    return 0;
}


label foamPetscMeshAssembly::InsertFvmGradPGaussIntoPETScMatrix
(
    const volScalarField& p,
    Mat jac,
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns,
    const scalar scale
) const
{
    // Gauss-style discrete gradient of p, inserted as the J_Dp block.
    // Interior face contribution to V*grad(p)|_own_i :
    //   = Sf_i * (w * p[own] + (1-w) * p[nei])
    // The scalar matrix is built identically to InsertFvmDivUIntoPETScMatrix
    // (same upper/lower/negSumDiag coefficients per face). The difference
    // is the destination block position: this helper inserts at
    //   (rowOffset + cmptI, colOffset)
    // i.e. cmptI varies the displacement-equation ROW (D_x, D_y, D_z),
    // while InsertFvmDivU varies the displacement-field COLUMN. Together
    // they produce the structurally consistent J_pD/J_Dp pair for the
    // Gauss discretisation of (-V*div D) in the pressure equation and
    // (-V*grad p) in the momentum equation.
    const fvMesh& mesh = p.mesh();
    const scalar sign = scale;

    for (label cmptI = 0; cmptI < 3; ++cmptI)
    {
        if (nScalarEqns == 2 && cmptI == 2)
        {
            break;
        }

        fvScalarMatrix gradPCoeffs(p, dimArea*dimPressure);

        const vectorField& Sf = mesh.Sf();
        const surfaceScalarField& weights = mesh.weights();
        const scalarField& w = weights;

        scalarField& upper = gradPCoeffs.upper();
        scalarField& lower = gradPCoeffs.lower();

        // Interior coefficients (mirror InsertFvmDivU)
        lower = sign*w*Sf.component(cmptI);
        upper = lower - sign*Sf.component(cmptI);

        gradPCoeffs.negSumDiag();

        // Boundary contributions: use p's BC (not U's)
        scalarField& diag = gradPCoeffs.diag();
        forAll(mesh.boundary(), patchI)
        {
            const fvPatchScalarField& pP = p.boundaryField()[patchI];
            const fvPatch& patch = pP.patch();
            const vectorField& Sfp = patch.Sf();
            const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
            const labelUList& fc = patch.faceCells();

            if (patch.coupled())
            {
                if (patch.type() == "processor")
                {
                    scalarField& patchIntCoeffs =
                        gradPCoeffs.internalCoeffs()[patchI];
                    scalarField& patchBndCoeffs =
                        gradPCoeffs.boundaryCoeffs()[patchI];

                    forAll(pP, faceI)
                    {
                        patchIntCoeffs[faceI] =
                            -sign*pw[faceI]*Sfp[faceI][cmptI];
                        patchBndCoeffs[faceI] =
                            sign*(1.0 - pw[faceI])*Sfp[faceI][cmptI];
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Coupled boundary type " << patch.type()
                        << " not yet supported in "
                        << "InsertFvmGradPGaussIntoPETScMatrix"
                        << abort(FatalError);
                }
            }
            else
            {
                // Non-coupled boundary: use p's BC valueInternalCoeffs.
                // For zeroGradient p, valueInternalCoeffs = 1 → adds
                //   -sign*Sf[cmpt] to diag per boundary face.
                // For fixedValue p, valueInternalCoeffs = 0 → no
                //   additional correction (and the diag stays as set
                //   by negSumDiag, which is the correct V*grad value
                //   for that BC because the face value of p doesn't
                //   depend on the interior).
                const scalarField internalCoeffs
                (
                    pP.valueInternalCoeffs(pw)
                );

                forAll(pP, faceI)
                {
                    diag[fc[faceI]] -=
                        sign*internalCoeffs[faceI]*Sfp[faceI][cmptI];
                }
            }
        }

        // Insert at (rowOffset + cmptI, colOffset)
        foamPetscMeshAssembly::InsertFvMatrixIntoPETScMatrix<scalar>
        (
            gradPCoeffs, jac, rowOffset + cmptI, colOffset, 1
        );
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif // #ifdef USE_PETSC
