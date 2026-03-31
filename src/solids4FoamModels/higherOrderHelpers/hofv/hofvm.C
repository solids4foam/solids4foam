/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#if defined(USE_PETSC) && !defined(FOAMEXTEND)

#include <algorithm>

#include "hofvm.H"
#include "fvc.H"
#include "compatibilityFunctions.H"
#include "surfaceFields.H"
#include "petscErrorHandling.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "fixedGradientFvPatchFields.H"

#ifdef OPENFOAM_COM
    #include "symmetryPlanePolyPatch.H"
#endif

namespace Foam
{

namespace
{

inline label nDisplacementEqns(const fvMesh& mesh)
{
    return mesh.nGeometricD() == 2 ? 2 : 3;
}


inline void setTensorBlockValues
(
    List<PetscScalar>& values,
    const tensor& coeff,
    const label blockSize,
    const label nScalarEqns,
    const label rowOffset,
    const label colOffset
)
{
    std::fill(values.begin(), values.end(), 0.0);

    for (label i = 0; i < nScalarEqns; ++i)
    {
        for (label j = 0; j < nScalarEqns; ++j)
        {
            values[(i + rowOffset)*blockSize + j + colOffset] = coeff[i*3 + j];
        }
    }
}


inline void addTensorCoeff
(
    Mat matrix,
    List<PetscScalar>& values,
    const tensor& coeff,
    const PetscInt globalBlockRowI,
    const PetscInt globalBlockColI,
    const label blockSize,
    const label nScalarEqns,
    const label rowOffset,
    const label colOffset
)
{
    setTensorBlockValues
    (
        values,
        coeff,
        blockSize,
        nScalarEqns,
        rowOffset,
        colOffset
    );

    AssertPETSc
    (
        MatSetValuesBlocked
        (
            matrix,
            1,
            &globalBlockRowI,
            1,
            &globalBlockColI,
            values.cdata(),
            ADD_VALUES
        )
    );
}

} // End anonymous namespace


Foam::tensor Foam::hofvm::laplacianCoeff
(
    const scalar& gammaMagSf,
    const scalar& quadWeight,
    const vector& gradInterpCoeff,
    const vector& faceNormal
)
{
    return I*gammaMagSf*quadWeight*(gradInterpCoeff & faceNormal);
}


Foam::tensor Foam::hofvm::laplacianTransposeCoeff
(
    const scalar& gammaMagSf,
    const scalar& quadWeight,
    const vector& gradInterpCoeff,
    const vector& faceNormal
)
{
    const vector& c = gradInterpCoeff;
    const vector& n = faceNormal;

    return gammaMagSf*quadWeight
       * tensor
         (
             c.x()*n.x(), c.x()*n.y(), c.x()*n.z(),
             c.y()*n.x(), c.y()*n.y(), c.y()*n.z(),
             c.z()*n.x(), c.z()*n.y(), c.z()*n.z()
         );
}


Foam::tensor Foam::hofvm::laplacianTraceCoeff
(
    const scalar& gammaMagSf,
    const scalar& quadWeight,
    const vector& gradInterpCoeff,
    const vector& faceNormal
)
{
    const vector& c = gradInterpCoeff;
    const vector& n = faceNormal;

    return gammaMagSf*quadWeight
       * tensor
         (
             c.x()*n.x(), c.y()*n.x(), c.z()*n.x(),
             c.x()*n.y(), c.y()*n.y(), c.z()*n.y(),
             c.x()*n.z(), c.y()*n.z(), c.z()*n.z()
         );
}


Foam::tmp<Foam::fvVectorMatrix> Foam::hofvm::d2dt2(const volVectorField& D)
{
    // Just temporarly solution
    return fvmDdtVectorCompat(D);
}


namespace hofvm
{

static label hofvmLaplacianPETSc
(
    Mat matrix,
    const foamPetscSnesHelper& petscSnesHelper,
    const movingLeastSquares& mls,
    const volVectorField& D,
    const volScalarField& diffusivity,
    tensor (*calcCoeff)
    (
        const scalar& gammaMagSf,
        const scalar& quadWeight,
        const vector& gradInterpCoeff,
        const vector& faceNormal
    ),
    const label rowOffset,
    const label colOffset
)
{
    const fvMesh& mesh = D.mesh();

    const movingLeastSquaresStencil& stencil = mls.stencilData();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const scalarField& magSfI = mesh.magSf().internalField();
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const label nScalarEqns = nDisplacementEqns(mesh);

    // Diffusion coefficient linearly interpolated to face centres.
    // If Gamma is not constant, next step is to interpolate diffusivity to
    // quad points using  hofvc::interpolate
    const surfaceScalarField gamma(fvc::interpolate(diffusivity));
    const scalarField& gammaI = gamma.internalField();

    // Face quadrature points weights, stencils and gradient interpolation
    // coefficients
    const CompactListList<scalar>& faceQuadWeights =
        mls.quadrature().faceQuadWeights();
    const CompactListList<label>& stencils = stencil.facesStencil();
    const List<CompactListList<vector>>& gradCoeffs = mls.faceGradCoeffsData();

    // Get the blockSize
    label blockSize = -1;
    AssertPETSc(MatGetBlockSize(matrix, &blockSize));

    // Initialise the block coefficient
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    // Loop over internal faces
    forAll(owner, faceI)
    {
        const vector& faceNormal = n[faceI];
        const scalar gammaMagSf = magSfI[faceI]*gammaI[faceI];
        const label ownCellID = owner[faceI];
        const label neiCellID = neighbour[faceI];
        const PetscInt globalOwnRow =
            petscSnesHelper.globalCells().toGlobal(ownCellID);
        const PetscInt globalNeiRow =
            petscSnesHelper.globalCells().toGlobal(neiCellID);
        const labelUList stencil = stencils[faceI];

        forAll(faceQuadWeights[faceI], qpI)
        {
            const scalar quadPointW = faceQuadWeights[faceI][qpI];

            forAll(stencil, cI)
            {
                const PetscInt globalCellID = stencil[cI];
                const tensor coeff =
                    calcCoeff
                    (
                        gammaMagSf,
                        quadPointW,
                        gradCoeffs[faceI][qpI][cI],
                        faceNormal
                    );

                addTensorCoeff
                (
                    matrix,
                    values,
                    coeff,
                    globalOwnRow,
                    globalCellID,
                    blockSize,
                    nScalarEqns,
                    rowOffset,
                    colOffset
                );

                addTensorCoeff
                (
                    matrix,
                    values,
                    -coeff,
                    globalNeiRow,
                    globalCellID,
                    blockSize,
                    nScalarEqns,
                    rowOffset,
                    colOffset
                );
            }
        }
    }

    // Loop over boundary faces
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const fvPatchVectorField& patchField = D.boundaryField()[patchI];

        if (isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        if (isA<processorPolyPatch>(pp))
        {
            const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
            const scalarField& pGamma = gamma.boundaryField()[patchI];
            const vectorField patchNormal(mesh.boundary()[patchI].nf());
            const label start = pp.start();

            forAll(pp, faceI)
            {
                const label faceID = faceI + start;
                const label ownCellID = mesh.faceOwner()[faceID];
                const PetscInt globalOwnRow =
                    petscSnesHelper.globalCells().toGlobal(ownCellID);
                const vector& faceNormal = patchNormal[faceI];
                const scalar gammaMagSf = pMagSf[faceI]*pGamma[faceI];
                const labelUList stencil = stencils[faceID];

                forAll(faceQuadWeights[faceID], qpI)
                {
                    const scalar quadPointW = faceQuadWeights[faceID][qpI];

                    forAll(stencil, cI)
                    {
                        const PetscInt globalCellID = stencil[cI];
                        const tensor coeff =
                            calcCoeff
                            (
                                gammaMagSf,
                                quadPointW,
                                gradCoeffs[faceID][qpI][cI],
                                faceNormal
                            );

                        addTensorCoeff
                        (
                            matrix,
                            values,
                            coeff,
                            globalOwnRow,
                            globalCellID,
                            blockSize,
                            nScalarEqns,
                            rowOffset,
                            colOffset
                        );
                    }
                }
            }

            continue;
        }

        if (pp.coupled())
        {
            // Cyclic not supprted but can be added if needed
            FatalErrorInFunction
                << "Coupled boundary patches are not implemented for the "
                << "high-order Jacobian."
                << abort(FatalError);
        }

        const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
        const scalarField& pGamma = gamma.boundaryField()[patchI];
        const vectorField patchNormal(mesh.boundary()[patchI].nf());
        const label start = pp.start();

        if
        (
            isA<symmetryPolyPatch>(pp)
#ifdef OPENFOAM_COM
         || isA<symmetryPlanePolyPatch>(pp)
#endif
        )
        {
            forAll(pp, faceI)
            {
                const label faceID = faceI + start;
                const label ownCellID = mesh.faceOwner()[faceID];
                const PetscInt globalOwnRow =
                    petscSnesHelper.globalCells().toGlobal(ownCellID);

                // Householder reflection matrix
                const vector& faceNormal = patchNormal[faceI];
                const tensor R = I - 2.0*sqr(faceNormal);

                const scalar gammaMagSf = pMagSf[faceI]*pGamma[faceI];
                const labelUList stencil = stencils[faceID];
                const label stencilSize = stencil.size();

                forAll(faceQuadWeights[faceID], qpI)
                {
                    const scalar quadPointW = faceQuadWeights[faceID][qpI];

                    for (label cI = 0; cI < stencilSize; ++cI)
                    {
                        const PetscInt globalCellID = stencil[cI];

                        const tensor coeff =
                            calcCoeff
                            (
                                gammaMagSf,
                                quadPointW,
                                gradCoeffs[faceID][qpI][cI],
                                faceNormal
                            );

                        const tensor mirrorCoeff =
                            calcCoeff
                            (
                                gammaMagSf,
                                quadPointW,
                                gradCoeffs[faceID][qpI][cI + stencilSize],
                                faceNormal
                            ) & R;

                        addTensorCoeff
                        (
                            matrix,
                            values,
                            coeff,
                            globalOwnRow,
                            globalCellID,
                            blockSize,
                            nScalarEqns,
                            rowOffset,
                            colOffset
                        );

                        addTensorCoeff
                        (
                            matrix,
                            values,
                            mirrorCoeff,
                            globalOwnRow,
                            globalCellID,
                            blockSize,
                            nScalarEqns,
                            rowOffset,
                            colOffset
                        );
                    }
                }
            }

            continue;
        }

        if (patchField.fixesValue())
        {
            forAll(pp, faceI)
            {
                const label faceID = faceI + start;
                const label ownCellID = mesh.faceOwner()[faceID];
                const PetscInt globalOwnRow =
                    petscSnesHelper.globalCells().toGlobal(ownCellID);
                const vector& faceNormal = patchNormal[faceI];
                const scalar gammaMagSf = pMagSf[faceI]*pGamma[faceI];
                const labelUList stencil = stencils[faceID];

                forAll(faceQuadWeights[faceID], qpI)
                {
                    const scalar quadPointW = faceQuadWeights[faceID][qpI];

                    forAll(stencil, cI)
                    {
                        const PetscInt globalCellID = stencil[cI];
                        const tensor coeff =
                            calcCoeff
                            (
                                gammaMagSf,
                                quadPointW,
                                gradCoeffs[faceID][qpI][cI],
                                faceNormal
                            );

                        addTensorCoeff
                        (
                            matrix,
                            values,
                            coeff,
                            globalOwnRow,
                            globalCellID,
                            blockSize,
                            nScalarEqns,
                            rowOffset,
                            colOffset
                        );
                    }
                }
            }

            continue;
        }

        if (isA<fixedGradientFvPatchVectorField>(patchField))
        {
            continue;
        }

        FatalErrorInFunction
            << "Unsupported patch type '" << pp.type()
            << "' for high-order Jacobian assembly on patch '"
            << pp.name() << "'."
            << abort(FatalError);
    }

    return 0;
}

} // End namespace hofvm


label Foam::hofvm::initialiseJacobian
(
    Mat& jac,
    const foamPetscSnesHelper& petscSnesHelper,
    const movingLeastSquares& mls,
    const volVectorField& D,
    const label blockSize,
    const bool createMat
)
{
    const fvMesh& mesh = D.mesh();

    if (jac)
    {
        FatalErrorInFunction
            << "Jacobian matrix already initialised"
            << abort(FatalError);
    }

    const label blockn = mesh.nCells();
    const PetscInt bs = static_cast<PetscInt>(blockSize);
    const PetscInt n = static_cast<PetscInt>(blockn)*bs;
    const PetscInt N =
        static_cast<PetscInt>(returnReduce(label(n), sumOp<label>()));

    if (createMat)
    {
        AssertPETSc(MatCreate(PETSC_COMM_WORLD, &jac));
        AssertPETSc(MatSetSizes(jac, n, n, N, N));
        AssertPETSc(MatSetType(jac, MATMPIBAIJ));
        AssertPETSc(MatSetBlockSize(jac, blockSize));
        AssertPETSc(MatSetFromOptions(jac));
    }

    const CompactListList<label>& stencils = mls.stencilData().facesStencil();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    List<labelHashSet> rowCols(blockn);

    forAll(rowCols, cellI)
    {
        rowCols[cellI] = labelHashSet(64);
        rowCols[cellI].insert
        (
            petscSnesHelper.globalCells().toGlobal(cellI)
        );
    }

    forAll(owner, faceI)
    {
        const label ownCellID = owner[faceI];
        const label neiCellID = neighbour[faceI];
        const labelUList faceStencil = stencils[faceI];

        forAll(faceStencil, cI)
        {
            const label globalCellID = faceStencil[cI];
            rowCols[ownCellID].insert(globalCellID);
            rowCols[neiCellID].insert(globalCellID);
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const fvPatchVectorField& patchField = D.boundaryField()[patchI];

        if (isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        if (isA<processorPolyPatch>(pp))
        {
            const label start = pp.start();

            forAll(pp, faceI)
            {
                const label faceID = faceI + start;
                const label ownCellID = mesh.faceOwner()[faceID];
                const labelUList faceStencil = stencils[faceID];

                forAll(faceStencil, cI)
                {
                    rowCols[ownCellID].insert(faceStencil[cI]);
                }
            }

            continue;
        }

        // Cyclic not supported
        if (pp.coupled())
        {
            FatalErrorInFunction
                << "Coupled boundary patches are not implemented for the "
                << "high-order Jacobian."
                << abort(FatalError);
        }

        if
        (
            patchField.fixesValue()
         || isA<symmetryPolyPatch>(pp)
#ifdef OPENFOAM_COM
         || isA<symmetryPlanePolyPatch>(pp)
#endif
        )
        {
            const label start = pp.start();

            forAll(pp, faceI)
            {
                const label faceID = faceI + start;
                const label ownCellID = mesh.faceOwner()[faceID];
                const labelUList faceStencil = stencils[faceID];

                forAll(faceStencil, cI)
                {
                    rowCols[ownCellID].insert(faceStencil[cI]);
                }
            }
        }
    }

#ifdef OPENFOAM_ORG
    const label myProcNo = Pstream::myProcNo();
    const label gStart = petscSnesHelper.globalCells().offset(myProcNo);
    const label gEnd =
        gStart + petscSnesHelper.globalCells().localSize(myProcNo);
#else
    const label gStart = petscSnesHelper.globalCells().localStart();
    const label gEnd = petscSnesHelper.globalCells().localEnd();
#endif
    List<PetscInt> d_nnz(blockn, 0);
    List<PetscInt> o_nnz(blockn, 0);

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

    AssertPETSc
    (
        MatMPIBAIJSetPreallocation
        (
            jac,
            blockSize,
            0,
            d_nnz.data(),
            0,
            o_nnz.data()
        )
    );

    AssertPETSc(MatSetOption(jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    AssertPETSc(MatSetOption(jac, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
    AssertPETSc(MatSetOption(jac, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    return 0;
}


void Foam::hofvm::laplacianIntoPETScMatrix
(
    Mat jac,
    const foamPetscSnesHelper& petscSnesHelper,
    const movingLeastSquares& mls,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const label rowOffset,
    const label colOffset
)
{
    hofvm::hofvmLaplacianPETSc
    (
        jac,
        petscSnesHelper,
        mls,
        D,
        diffusivity,
        hofvm::laplacianCoeff,
        rowOffset,
        colOffset
    );
}


void Foam::hofvm::laplacianTransposeIntoPETScMatrix
(
    Mat jac,
    const foamPetscSnesHelper& petscSnesHelper,
    const movingLeastSquares& mls,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const label rowOffset,
    const label colOffset
)
{
    hofvm::hofvmLaplacianPETSc
    (
        jac,
        petscSnesHelper,
        mls,
        D,
        diffusivity,
        hofvm::laplacianTransposeCoeff,
        rowOffset,
        colOffset
    );
}


void Foam::hofvm::laplacianTraceIntoPETScMatrix
(
    Mat jac,
    const foamPetscSnesHelper& petscSnesHelper,
    const movingLeastSquares& mls,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const label rowOffset,
    const label colOffset
)
{
    hofvm::hofvmLaplacianPETSc
    (
        jac,
        petscSnesHelper,
        mls,
        D,
        diffusivity,
        hofvm::laplacianTraceCoeff,
        rowOffset,
        colOffset
    );
}

} // End namespace Foam

#endif

// ************************************************************************* //
