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

#ifdef OPENFOAM_COM

#include "vertexCentredLinGeomSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrix.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "sparseMatrixTools.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#include "pointFieldFunctions.H"
#ifdef USE_PETSC
    #include <petscksp.h>
    #include <petscsnes.h>
#endif

// * * * * * * * * * * * * * * External Functions  * * * * * * * * * * * * * //

#ifdef USE_PETSC

// User data "context" for PETSc functions
typedef struct appCtx
{
    // Reference to the solid model object
    Foam::solidModels::vertexCentredLinGeomSolid& solMod_;

    // Constructor
    appCtx
    (
        Foam::solidModels::vertexCentredLinGeomSolid& solMod
    )
    :
        solMod_(solMod)
    {}
} appCtx;


PetscErrorCode formResidual
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Vec f,        // residual
    void *ctx     // user context
)
{
    const PetscScalar *xx;
    PetscScalar       *ff;
    appCtx *user = (appCtx *)ctx;

    // Access x and f data
    // VecRestoreArrayRead() and VecRestoreArray() must be called when access is
    // no longer needed
    CHKERRQ(VecGetArrayRead(x,&xx));
    CHKERRQ(VecGetArray(f,&ff));

    // Map the solution to an OpenFOAM field
    const bool twoD = user->solMod_.twoD();
    Foam::pointVectorField pointD("petsc_pointD", user->solMod_.pointD());
    Foam::vectorField& pointDI = pointD;
    const Foam::boolList& ownedByThisProc =
        user->solMod_.gPointIndices().ownedByThisProc();
    {
        int index = 0;
        forAll(pointDI, i)
        {
            if (ownedByThisProc[i])
            {
                pointDI[i].x() = xx[index++];
                pointDI[i].y() = xx[index++];

                if (!twoD)
                {
                    pointDI[i].z() = xx[index++];
                }
            }
        }
    }

    if (Foam::Pstream::parRun())
    {
        // I need to get values not owned by this proc
        Foam::FatalError
            << "Fix solution retrieval" << Foam::abort(Foam::FatalError);
    }

    // Don't call correctBCs as we want to allow perturbation of known DOFs
    //pointD.correctBoundaryConditions();

    // Compute the residual
    const Foam::vectorField res(user->solMod_.residualMomentum(pointD));

    // Map the data to f
    const Foam::labelList& localToGlobalPointMap =
        user->solMod_.gPointIndices().localToGlobalPointMap();
    const int blockSize = twoD ? 2 : 3;
    forAll(res, localBlockRowI)
    {
        const Foam::vector& resI = res[localBlockRowI];
        const int blockRowI = localToGlobalPointMap[localBlockRowI];

        ff[blockRowI*blockSize] = resI.x();
        ff[blockRowI*blockSize + 1] = resI.y();

        if (!twoD)
        {
            ff[blockRowI*blockSize + 2] = resI.z();
        }
    }

    if (Foam::Pstream::parRun())
    {
        // What happens for points not owned by this proc?
        // I guess the source should be summed across procs
        Foam::FatalError
            << "Form residual: do we need to use owned by?"
            << abort(Foam::FatalError);
    }

    // Restore vectors
    CHKERRQ(VecRestoreArrayRead(x,&xx));
    CHKERRQ(VecRestoreArray(f,&ff));

    return 0;
}


PetscErrorCode formJacobian
(
    SNES snes,    // snes object
    Vec x,        // current solution
    Mat A,        // Jacobian
    Mat B,        // Jaconian precondioner (can be A)
    void *ctx     // user context
)
{
    // Get pointer to solution data
    const PetscScalar *xx;
    CHKERRQ(VecGetArrayRead(x, &xx));

    // Map the solution to an OpenFOAM field
    appCtx *user = (appCtx *)ctx;
    const bool twoD = user->solMod_.twoD();
    Foam::pointVectorField pointD("petsc_pointD", user->solMod_.pointD());
    Foam::vectorField& pointDI = pointD;
    const Foam::boolList& ownedByThisProc =
        user->solMod_.gPointIndices().ownedByThisProc();
    {
        int index = 0;
        forAll(pointDI, i)
        {
            if (ownedByThisProc[i])
            {
                pointDI[i].x() = xx[index++];
                pointDI[i].y() = xx[index++];

                if (!twoD)
                {
                    pointDI[i].z() = xx[index++];
                }
            }
        }
    }

    if (Foam::Pstream::parRun())
    {
        // I need to get values not owned by this proc
        Foam::FatalError
            << "Fix solution retrieval" << abort(Foam::FatalError);
    }

    // This may not be needed
    pointD.correctBoundaryConditions();

    // Restore solution vector
    CHKERRQ(VecRestoreArrayRead(x, &xx));


    // Compute Jacobian in OpenFOAM format
    Foam::sparseMatrix matrix;
    matrix += user->solMod_.JacobianMomentum(pointD)();


    // Set matrix coefficients, if any, to zero
    MatInfo info;
    MatGetInfo(A, MAT_LOCAL, &info);
    if (info.nz_used)
    {
        CHKERRQ(MatZeroEntries(A));
    }


    // Insert OpenFOAM matrix into PETSc matrix
    // Note: we use global indices when inserting coefficients
    const Foam::sparseMatrixData& data = matrix.data();
    const Foam::labelList& localToGlobalPointMap =
        user->solMod_.gPointIndices().localToGlobalPointMap();
    PetscScalar values2d[4];
    PetscScalar values3d[9];
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
        const Foam::tensor& coeff = iter();
        const int blockRowI = localToGlobalPointMap[iter.key()[0]];
        const int blockColI = localToGlobalPointMap[iter.key()[1]];

        if (twoD)
        {
            // Prepare values
            values2d[0] = coeff.xx();
            values2d[1] = coeff.xy();
            values2d[2] = coeff.yx();
            values2d[3] = coeff.yy();

            // Insert tensor coefficient
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    A, 1, &blockRowI, 1, &blockColI, values2d, ADD_VALUES
                )
            );
        }
        else // 3-D
        {
            // Prepare values
            values3d[0] = coeff.xx();
            values3d[1] = coeff.xy();
            values3d[2] = coeff.xz();
            values3d[3] = coeff.yx();
            values3d[4] = coeff.yy();
            values3d[5] = coeff.yz();
            values3d[6] = coeff.zx();
            values3d[7] = coeff.zy();
            values3d[8] = coeff.zz();

            // Insert tensor coefficient
            CHKERRQ
            (
                MatSetValuesBlocked
                (
                    A, 1, &blockRowI, 1, &blockColI, values3d, ADD_VALUES
                )
            );
        }
    }
    if (Foam::Pstream::parRun())
    {
        Foam::FatalError
            << "formJac: check if owned by is needed" << abort(Foam::FatalError);
    }

    // Complete matrix assembly
    CHKERRQ(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    if (A != B)
    {
        CHKERRQ(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
        CHKERRQ(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
    }

    return 0;
}
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredLinGeomSolid, 0);
addToRunTimeSelectionTable(solidModel, vertexCentredLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void vertexCentredLinGeomSolid::updatePointDivSigma
(
    const pointVectorField& pointD,
    surfaceTensorField& dualGradDf,
    surfaceSymmTensorField& dualSigmaf,
    pointVectorField& pointDivSigma
)
{
    if (debug)
    {
        Info<< "void explicitVertexCentredLinGeomSolid::"
            << "updatePointDivSigma(...): start" << endl;
    }

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.0));

    // Calculate gradD at dual faces
    dualGradDf = vfvc::fGrad
    (
        pointD,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta,
        debug
    );

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Calculate the tractions at the dual faces
    surfaceVectorField dualTraction
    (
        (dualMesh().Sf()/dualMesh().magSf()) & dualSigmaf
    );

    // Enforce extract tractions on traction boundaries
    enforceTractionBoundaries
    (
        pointD, dualTraction, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
        }
    }

    // Calculate divergence of stress (force per unit volume) for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh().magSf());

    // Calculate absolute divergence of stress (force)
    // We do this to allow syncing of forces at points on processor boundaries
    const vectorField dualDivSigmaAbs(dualDivSigma*dualMesh().V());

    // Map dual cell field to primary mesh point field
    // We temporarily use the pointDivSigma field to hold absolute forces
    // but convert them back to force per unit volume below
    vectorField& pointDivSigmaI = pointDivSigma;
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualDivSigmaAbs, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigmaI[pointID] = dualDivSigmaAbs[dualCellI];
    }

    // Sum absolute forces in parallel
    pointConstraints::syncUntransformedData
    (
        mesh(), pointDivSigma, plusEqOp<vector>()
    );

    // Convert force to force per unit volume
    // Perform calculation per point to avoid dimension checks
    const scalarField& pointGlobalVolI = pointGlobalVol_;
    forAll(pointDivSigmaI, pointI)
    {
        pointDivSigmaI[pointI] /= pointGlobalVolI[pointI];
    }

    if (debug)
    {
        Info<< "void explicitVertexCentredLinGeomSolid::"
            << " updatePointDivSigma(...): end" << endl;
    }
}


void vertexCentredLinGeomSolid::setFixedDofs
(
    const pointVectorField& pointD,
    boolList& fixedDofs,
    pointField& fixedDofValues,
    symmTensorField& fixedDofDirections
) const
{
    // Flag all fixed DOFs

    forAll(pointD.boundaryField(), patchI)
    {
        if
        (
            // isA<uniformFixedValuePointPatchVectorField>
            isA<fixedValuePointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            // const uniformFixedValuePointPatchVectorField& dispPatch =
            //     refCast<const uniformFixedValuePointPatchVectorField>
            // const fixedValuePointPatchVectorField& dispPatch =
            //     refCast<const fixedValuePointPatchVectorField>
            //     (
            //         pointD.boundaryField()[patchI]
            //     );

            // const vector& disp = dispPatch.uniformValue();

            const labelList& meshPoints =
                pointD.mesh().mesh().boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                const vector& disp = pointD[pointID];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the new one
                    if
                    (
                        mag
                        (
                            fixedDofDirections[pointID]
                          & (fixedDofValues[pointID] - disp)
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredLinGeomSolid::setFixedDofs(...)"
                        )   << "Inconsistent displacements prescribed at point "
                            << "= " << pointD.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // Set all directions as fixed, just in case it was
                    // previously marked as a symmetry point
                    fixedDofDirections[pointID] = symmTensor(I);
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = disp;
                    fixedDofDirections[pointID] = symmTensor(I);
                }
            }
        }
        else if
        (
            isA<symmetryPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointD.mesh().boundary()[patchI].meshPoints();
            const vectorField& pointNormals =
                pointD.mesh().boundary()[patchI].pointNormals();

            scalarField normalDisp(meshPoints.size(), 0.0);
            if
            (
                isA<fixedDisplacementZeroShearPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                )
            )
            {
                normalDisp =
                (
                    pointNormals
                  & pointD.boundaryField()[patchI].patchInternalField()
                );

                if (debug)
                {
                    Info<< "normalDisp = " << normalDisp << endl;
                }
            }

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the current condition
                    if
                    (
                        mag
                        (
                            (pointNormals[pI] & fixedDofValues[pointID])
                          - normalDisp[pI]
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredLinGeomSolid::setFixedDofs(...)"
                        )   << "Inconsistent displacements prescribed at point "
                            << "= " << pointD.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // If the point is not fully fixed then make sure the normal
                    // direction is fixed
                    if (mag(fixedDofDirections[pointID] - symmTensor(I)) > 0)
                    {
                        // If the directions are orthogonal we can add them
                        const symmTensor curDir = sqr(pointNormals[pI]);
                        if (mag(fixedDofDirections[pointID] & curDir) > 0)
                        {
                            FatalError
                                << "Point " << pointID << " is fixed in two "
                                << "directions: this is only implemented for "
                                << "Cartesian axis directions" << abort(FatalError);
                        }

                        fixedDofDirections[pointID] += curDir;
                    }
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = normalDisp[pI]*pointNormals[pI];
                    fixedDofDirections[pointID] = sqr(pointNormals[pI]);
                }
            }
        }
    }
}


void vertexCentredLinGeomSolid::enforceTractionBoundaries
(
    const pointVectorField& pointD,
    surfaceVectorField& dualTraction,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    const pointMesh& pMesh = pointD.mesh();
    const fvMesh& dualMesh = dualTraction.mesh();

    forAll(pointD.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const solidTractionPointPatchVectorField& tracPatch =
                refCast<const solidTractionPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                );

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // Primary mesh point normals
            const vectorField& n =
                pMesh.boundary()[patchI].pointNormals();

            // Primary mesh point tractions
            const vectorField totalTraction
            (
                tracPatch.traction() - n*tracPatch.pressure()
            );

            // Create dual mesh faces traction field
            vectorField dualFaceTraction
            (
                dualMesh.boundaryMesh()[patchI].size(), vector::zero
            );

            // Multiple points map to each dual face so we will count them
            // and then divide the dualFaceTraction by this field so that it is
            // the average of all the points that map to it
            scalarField nPointsPerDualFace(dualFaceTraction.size(), 0.0);

            // Map from primary mesh point field to second mesh face field using
            // the pointToDualFaces map
            forAll(totalTraction, pI)
            {
                const label pointID = meshPoints[pI];
                const labelList& curDualFaces = pointToDualFaces[pointID];

                forAll(curDualFaces, dfI)
                {
                    const label dualFaceID = curDualFaces[dfI];

                    if (!dualMesh.isInternalFace(dualFaceID))
                    {
                        // Check which patch this dual face belongs to
                        const label dualPatchID =
                            dualMesh.boundaryMesh().whichPatch(dualFaceID);

                        if (dualPatchID == patchI)
                        {
                            // Find local face index
                            const label localDualFaceID =
                                dualFaceID
                              - dualMesh.boundaryMesh()[dualPatchID].start();

                            // Set dual face traction
                            dualFaceTraction[localDualFaceID] +=
                                totalTraction[pI];

                            // Update the count for this face
                            nPointsPerDualFace[localDualFaceID]++;
                        }
                    }
                }
            }

            if (gMin(nPointsPerDualFace) < 1)
            {
                FatalErrorIn
                (
                    "void vertexCentredLinGeomSolid::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;

            // Overwrite the dual patch face traction
            dualTraction.boundaryFieldRef()[patchI] = dualFaceTraction;
        }
        else if
        (
            isA<symmetryPointPatchVectorField>(pointD.boundaryField()[patchI])
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            // Set the dual patch face shear traction to zero
            const vectorField n(dualMesh.boundary()[patchI].nf());
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
        }
    }

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
        }
    }
}

bool vertexCentredLinGeomSolid::vertexCentredLinGeomSolid::converged
(
    const label iCorr,
    scalar& initResidual,
    const scalar res,
    const label nInterations,
    const pointVectorField& pointD,
    const vectorField& pointDcorr
) const
{
    // Calculate the residual as the root mean square of the correction
    const scalar residualAbs = gSum(magSqr(pointDcorr));

    // Store initial residual
    if (iCorr == 0)
    {
        initResidual = residualAbs;

        // If the initial residual is small then convergence has been achieved
        if (initResidual < SMALL)
        {
            Info<< "    Initial residual is less than " << SMALL << nl
                << "    Converged" << endl;
            return true;
        }
        Info<< "    Initial residual = " << initResidual << endl;
    }

    // Define a normalised residual wrt the initial residual
    const scalar residualNorm = residualAbs/initResidual;

    // Calculate the maximum displacement
    const scalar maxMagD = gMax(mag(pointD.primitiveField()));

    // Check for convergence
    bool converged = false;
    if (residualNorm < solutionTol())
    {
        Info<< "    Converged" << endl;
        converged = true;
    }
    else if (residualAbs < VSMALL)
    {
        Info<< "    Converged: absolute residual is less than " << VSMALL
            << endl;
        converged = true;
    }

    if (iCorr == 0)
    {
        Info<< "    Corr, res, relRes, resAbs, iters, maxMagD" << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged || iCorr >= nCorr() - 1)
    {
        Info<< "    " << iCorr
            << ", " << res
            << ", " << residualNorm
            << ", " << residualAbs
            << ", " << nInterations
            << ", " << maxMagD << endl;

        if (iCorr >= nCorr())
        {
            Warning
                << "Max iterations reached within the momentum loop"
                << endl;
            converged = true;
        }
    }

    return converged;
}


scalar vertexCentredLinGeomSolid::calculateLineSearchSlope
(
    const scalar eta,
    const vectorField& pointDcorr,
    pointVectorField& pointD,
    surfaceTensorField& dualGradDf,
    surfaceSymmTensorField& dualSigmaf,
    const scalar zeta
)
{
    // Store pointD as we will reset it after changing it
    pointD.storePrevIter();

    // Update pointD
    pointD.primitiveFieldRef() += eta*pointDcorr;
    pointD.correctBoundaryConditions();

    // Calculate gradD at dual faces
    dualGradDf = vfvc::fGrad
    (
        pointD,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta
    );

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Update the source vector
    pointD.correctBoundaryConditions();
    const vectorField source(-residualMomentum(pointD));

    // Reset pointD
    pointD = pointD.prevIter();

    // Return the slope
    return gSum(pointDcorr & source);
}


scalar vertexCentredLinGeomSolid::calculateLineSearchFactor
(
    const scalar rTol, // Slope reduction tolerance
    const int maxIter, // Maximum number of line search iterations
    const vectorField& pointDcorr, // Point displacement correction
    const vectorField& source, // Linear system source
    const scalar zeta // Discretisation parameter
)
{
    // Calculate initial slope
    const scalar s0 = gSum(pointDcorr & source);

    // Set initial line search parameter
    scalar eta = 1.0;
    int lineSearchIter = 0;

    // Perform backtracking to find suitable eta
    do
    {
        lineSearchIter++;

        // Calculate slope at eta
        const scalar s1 = calculateLineSearchSlope
        (
            eta, pointDcorr, pointD(), dualGradDf_, dualSigmaf_, zeta
        );

        // Calculate ratio of s1 to s0
        const scalar r = s1/s0;

        if (mag(r) < rTol)
        {
            break;
        }
        else
        {
            // Interpolate/extrapolate to find new eta
            // Limit it to be less than 10
            //eta = min(-1/(r - 1), 10);

            if (r < 0)
            {
                // Simple back tracking
                eta *= 0.5;
            }
            else
            {
                // Extrapolate
                eta = min(-1/(r - 1), 10);
            }
        }

        if (lineSearchIter == maxIter)
        {
            Warning
                << "Max line search iterations reached!" << endl;
        }
    }
    while (lineSearchIter < maxIter);

    // Update pointD and re-calculate source, then calculate s
    if (mag(eta - 1) > SMALL)
    {
        Info<< "        line search parameter = " << eta
            << ", iter = " << lineSearchIter << endl;
    }

    return eta;
}


void vertexCentredLinGeomSolid::makeDualImpKf() const
{
    if (dualImpKfPtr_.valid())
    {
        FatalErrorIn("void vertexCentredLinGeomSolid::makeDualImpKf() const")
            << "Pointer already set!" << abort(FatalError);
    }

    dualImpKfPtr_.set
    (
        new surfaceScalarField(dualMechanicalPtr_().impKf())
    );
}


const surfaceScalarField& vertexCentredLinGeomSolid::dualImpKf() const
{
    if (dualImpKfPtr_.empty())
    {
        makeDualImpKf();
    }

    return dualImpKfPtr_();
}


void vertexCentredLinGeomSolid::predict()
{
    Info<< "Predicting pointD" << endl;

    if (true)
    {
        // Assuming constant velocity
        pointD() = pointD().oldTime() + pointU_*runTime().deltaT();
    }
    else
    {
        // Assuming constant acceleration
        pointD() =
            pointD().oldTime()
          + pointU_*runTime().deltaT()
          + 0.5*pointA_*pow(runTime().deltaT(), 2);
    }
}


bool vertexCentredLinGeomSolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Evolving solid solver" << endl;

    // Update boundary conditions
    pointD().correctBoundaryConditions();

    // Create user data context
    appCtx user(*this);

    // Lookup the PETSc options file
    fileName optionsFile(solidModelDict().lookup("optionsFile"));
    optionsFile.expand();

    // Initialise PETSc with an options file
    PetscInitialize(NULL, NULL, optionsFile.c_str(), NULL);

    // Create a SNES solver context
    SNES snes;
    SNESCreate(PETSC_COMM_WORLD, &snes);

    // Set the user context
    SNESSetApplicationContext(snes, &user);

    // Set the residual function
    SNESSetFunction(snes, NULL, formResidual, &user);

    // Initialise the Jacobian matrix

    const int blockSize = twoD_ ? 2 : 3;
    const Foam::labelList& localToGlobalPointMap =
        globalPointIndices_.localToGlobalPointMap();

    // Find size of global system, i.e. the highest global point index + 1
    const int blockN = Foam::gMax(localToGlobalPointMap) + 1;
    const int N = blockSize*blockN;

    // Find the start and end global point indices for this proc
    int blockStartID = N;
    int blockEndID = -1;
    const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    forAll(ownedByThisProc, pI)
    {
        if (ownedByThisProc[pI])
        {
            blockStartID = Foam::min(blockStartID, localToGlobalPointMap[pI]);
            blockEndID = Foam::max(blockEndID, localToGlobalPointMap[pI]);
        }
    }
    //const int startID = blockSize*blockStartID;
    //const int endID = blockSize*(blockEndID + 1) - 1;

    // Find size of local system, i.e. the range of points owned by this proc
    const int blockn = blockEndID - blockStartID + 1;
    const int n = blockSize*blockn;

    // Create the matrix
    Mat A;
    CHKERRQ(MatCreate(PETSC_COMM_WORLD, &A));

    // Set matrix characteristics
    CHKERRQ(MatSetSizes(A, n, n, N, N));
    //CHKERRQ(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
    CHKERRQ(MatSetFromOptions(A));
    CHKERRQ(MatSetType(A, MATMPIAIJ));
    //CHKERRQ(MatSetType(A, MATMPIBAIJ));
    CHKERRQ(MatSetBlockSize(A, blockSize));

    // Count the number of non-zeros
    // TODO: fix for parallel
    int* d_nnz = (int*)malloc(n*sizeof(int));
    int* o_nnz = (int*)malloc(n*sizeof(int));
    // label d_nnz[n];
    // label o_nnz[n];
    int d_nz = 0;
    Foam::sparseMatrixTools::setNonZerosPerRow
    (
        d_nnz,
        o_nnz,
        d_nz,
        n,
        blockSize,
        ownedByThisProc,
        globalPointIndices_.stencilSizeOwned(),
        globalPointIndices_.stencilSizeNotOwned()
    );

    // Allocate parallel matrix
    CHKERRQ(MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz));
    // TODO: change d_nnz/o_nnz to block sizes!
    //CHKERRQ(MatMPIBAIJSetPreallocation(A, blockSize, 0, d_nnz, 0, o_nnz));

    // TO BE FIXED: some mallocs are still needed in parallel!
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    CHKERRQ(MatSetUp(A));

    // Do not call the matrix assembly as we have not inserted any values
    //CHKERRQ(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    //CHKERRQ(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    // Set the Jacobian function
    SNESSetJacobian(snes, A, A, formJacobian, &user);

    // Set solver options
    // Uses default options, can be overridden by command line options
    SNESSetFromOptions(snes);

    // // Find the start and end global point indices for this proc
    // forAll(ownedByThisProc, pI)
    // {
    //     if (ownedByThisProc[pI])
    //     {
    //         blockStartID = min(blockStartID, localToGlobalPointMap[pI]);
    //         blockEndID = max(blockEndID, localToGlobalPointMap[pI]);
    //     }
    // }
    //const label startID = blockSize*blockStartID;
    //const label endID = blockSize*(blockEndID + 1) - 1;

    // Create the solution vector
    Vec x;
    CHKERRQ(VecCreate(PETSC_COMM_WORLD, &x));
    CHKERRQ(VecSetSizes(x, n, N));
    CHKERRQ(VecSetBlockSize(x, blockSize));
    CHKERRQ(VecSetType(x, VECMPI));
    CHKERRQ(PetscObjectSetName((PetscObject) x, "Solution"));
    CHKERRQ(VecSetFromOptions(x));

    // Set initial guess
    CHKERRQ(VecSet(x, 0.0));

    // Solve the nonlinear system
    CHKERRQ(SNESSolve(snes, NULL, x));

    // Check if SNES converged
    SNESConvergedReason reason;
    CHKERRQ(SNESGetConvergedReason(snes, &reason));
    if (reason < 0)
    {
        FatalErrorIn("vertexCentredLinGeomSolid::evolveSnes()")
            << "The SNES nonlinear solver did not converge."
            << " PETSc SNES error code = " << reason << abort(FatalError);
    }

    // Retrieve the solution
    const PetscScalar *xx;
    VecGetArrayRead(x,&xx);
    Foam::vectorField& pointDI = pointD();
    {
        int index = 0;
        forAll(pointDI, i)
        {
            if (ownedByThisProc[i])
            {
                pointDI[i].x() = xx[index++];
                pointDI[i].y() = xx[index++];

                if (!twoD_)
                {
                    pointDI[i].z() = xx[index++];
                }
            }
        }
    }
    if (Pstream::parRun())
    {
        // I need to get values not owned by this proc
        FatalError
            << "Fix solution retrieval" << abort(FatalError);
    }
    pointD().correctBoundaryConditions();

    // Destroy PETSc objects
    VecDestroy(&x);
    SNESDestroy(&snes);


    // Update point accelerations
    // Note: for NewmarkBeta, this needs to come before the pointU update
    pointA_.primitiveFieldRef() =
        vfvc::ddt
        (
            mesh().ddtScheme("ddt(pointU)"),
            mesh().d2dt2Scheme("d2dt2(pointD)"),
            pointU_
        );

    // Update point velocities
    pointU_.primitiveFieldRef() =
        vfvc::ddt
        (
            mesh().ddtScheme("ddt(pointD)"),
            mesh().d2dt2Scheme("d2dt2(pointD)"),
            pointD()
        );

    // Calculate gradD at dual faces
    dualGradDf_ = vfvc::fGrad
    (
        pointD(),
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta_,
        debug
    );

    // Update the increment of displacement
    pointDD() = pointD() - pointD().oldTime();

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1)
    {
        mechanical().mapGradToSubMeshes(gradD());
    }

    // Update dual face stress field
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate
    mechanical().correct(sigma());

    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());
#else
    FatalErrorIn("vertexCentredLinGeomSolid::evolveSnes()")
        << "PETSc not available. Please set the PETSC_DIR environment "
        << "variable and re-compile solids4foam" << abort(FatalError);
#endif

    return true;
}


bool vertexCentredLinGeomSolid::evolveImplicitCoupled()
{
    Info<< "Evolving solid solver" << endl;

    // Update boundary conditions
    pointD().correctBoundaryConditions();

    // Initialise matrix
    sparseMatrix matrix(sum(globalPointIndices_.stencilSize()));

    // Store material tangent field for dual mesh faces
    Field<scalarSquareMatrix> materialTangent
    (
        dualMechanicalPtr_().materialTangentFaceField()
    );

    // Solution field: point displacement correction
    vectorField pointDcorr(pointD().internalField().size(), vector::zero);

    // Newton-Raphson loop over momentum equation
    int iCorr = 0;
    scalar initResidual = 0.0;
    SolverPerformance<vector> solverPerf;
    do
    {
        // The source is the negative of the momentum residual
        pointD().correctBoundaryConditions();
        vectorField source(-residualMomentum(pointD()));

        if (fullNewton_ || iCorr == 0)
        {
            // Assemble the matrix once per outer iteration
            matrix.clear();

            // Update material tangent
            materialTangent = dualMechanicalPtr_().materialTangentFaceField();

            // Add div(sigma) coefficients
            vfvm::divSigma
            (
                matrix,
                mesh(),
                dualMesh(),
                dualMeshMap().dualFaceToCell(),
                dualMeshMap().dualCellToPoint(),
                materialTangent,
                zeta_
            );

            // Add d2dt2 coefficients
            vfvm::d2dt2
            (
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                runTime().deltaTValue(),
                pointD().name(),
                matrix,
                pointRho_.internalField(),
                pointVol_.internalField(),
                int(bool(debug))
            );
        }

        if (debug > 1)
        {
            // Print the matrix
            matrix.print();
        }

        // Enforce fixed DOF on the linear system
        sparseMatrixTools::enforceFixedDof
        (
            matrix,
            source,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofValues_,
            fixedDofScale_
        );

        if (debug > 1)
        {
            // Print the matrix
            matrix.print();
        }

        // Solve linear system for displacement correction
        if (debug)
        {
            Info<< "bool vertexCentredLinGeomSolid::evolve(): "
                << " solving linear system: start" << endl;
        }

        if (Switch(solidModelDict().lookup("usePETSc")))
        {
#ifdef USE_PETSC
            fileName optionsFile(solidModelDict().lookup("optionsFile"));
            solverPerf = sparseMatrixTools::solveLinearSystemPETSc
            (
                matrix,
                source,
                pointDcorr,
                twoD_,
                optionsFile,
                mesh().points(),
                globalPointIndices_.ownedByThisProc(),
                globalPointIndices_.localToGlobalPointMap(),
                globalPointIndices_.stencilSizeOwned(),
                globalPointIndices_.stencilSizeNotOwned(),
                solidModelDict().lookupOrDefault<bool>("debugPETSc", false)
            );
#else
            FatalErrorIn("vertexCentredLinGeomSolid::evolve()")
                << "PETSc not available. Please set the PETSC_DIR environment "
                << "variable and re-compile solids4foam" << abort(FatalError);
#endif
        }
        else
        {
            // Use Eigen SparseLU direct solver
            sparseMatrixTools::solveLinearSystemEigen
            (
                matrix, source, pointDcorr, twoD_, false, debug
            );
        }

        if (debug)
        {
            Info<< "bool vertexCentredLinGeomSolid::evolve(): "
                << " solving linear system: end" << endl;
        }

        // Update point displacement field
        if (Switch(solidModelDict().lookup("lineSearch")))
        {
            // Lookup target tolerance for slope reduction
            const scalar rTol
            (
                solidModelDict().lookupOrDefault<scalar>("lineSearchRTol", 0.8)
            );

            // Lookup the maximum number of line search iterations
            const int maxIter
            (
                solidModelDict().lookupOrDefault<scalar>
                (
                    "lineSearchMaxIter", 10
                )
            );

            // Calculate line search factor
            const scalar eta
            (
                calculateLineSearchFactor
                (
                    rTol, maxIter, pointDcorr, source, zeta_
                )
            );

            // Update displacement field
            pointD().primitiveFieldRef() += eta*pointDcorr;
        }
        else if (mesh().relaxField(pointD().name()))
        {
            // Relaxing the correction can help convergence

            const scalar rf
            (
                mesh().fieldRelaxationFactor(pointD().name())
            );

            pointD().primitiveFieldRef() += rf*pointDcorr;
        }
        else
        {
            pointD().primitiveFieldRef() += pointDcorr;
        }
        pointD().correctBoundaryConditions();

        // Update point accelerations
        // Note: for NewmarkBeta, this needs to come before the pointU update
        pointA_.primitiveFieldRef() =
            vfvc::ddt
            (
                mesh().ddtScheme("ddt(pointU)"),
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                pointU_
            );

        // Update point velocities
        pointU_.primitiveFieldRef() =
            vfvc::ddt
            (
                mesh().ddtScheme("ddt(pointD)"),
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                pointD()
            );
    }
    while
    (
        !converged
        (
            iCorr,
            initResidual,
            mag(solverPerf.finalResidual()),
            cmptMax(solverPerf.nIterations()),
            pointD(),
            pointDcorr
        ) && ++iCorr
    );

    // Calculate gradD at dual faces
    dualGradDf_ = vfvc::fGrad
    (
        pointD(),
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta_,
        debug
    );

    // Update the increment of displacement
    pointDD() = pointD() - pointD().oldTime();

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1)
    {
        mechanical().mapGradToSubMeshes(gradD());
    }

    // Update dual face stress field
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate
    mechanical().correct(sigma());

    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());

    return true;
}


bool vertexCentredLinGeomSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver" << endl;

    // Predict pointD
    if (predictor_ && newTimeStep())
    {
        predict();
    }

    // Initialise matrix
    sparseScalarMatrix matrixNoBCs(sum(globalPointIndices_.stencilSize()));

    // Lookup flag to indicate compact or large Laplacian stencil
    const Switch compactImplicitStencil
    (
        solidModelDict().lookupOrDefault<Switch>("compactImplicitStencil", true)
    );
    Info<< "compactImplicitStencil: " << compactImplicitStencil << endl;

    // Create scalar Laplacian discretisation matrix without boundary conditions
    vfvm::laplacian
    (
        matrixNoBCs,
        compactImplicitStencil,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        dualImpKf().primitiveField(),
        debug
    );

    // // Global point index lists
    // const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    // const labelList& localToGlobalPointMap =
    //     globalPointIndices_.localToGlobalPointMap();

    // Solution field: point displacement componnet correction
    scalarField pointDcorr(pointD().internalField().size(), 0.0);

    // Vector field version of pointDcorr
    vectorField pointDcorrVec(pointD().internalField().size(), vector::zero);

    // Initialise the source
    vectorField source(mesh().nPoints(), vector::zero);

    // Outer loop over momentum equation
    int iCorr = 0;
    scalar initResidual = 0.0;
    SolverPerformance<vector> solverPerf;
    do
    {
        // Store previous iteration of pointD for residual calculation
        pointD().storePrevIter();

        // Update the source vector
        pointD().correctBoundaryConditions();
        source = -residualMomentum(pointD());

        // Loop over solution directions (e.g., x, y, z)
        forAll(mesh().solutionD(), dirI)
        {
            if (mesh().solutionD()[dirI] < 0)
            {
                // Empty direction
                continue;
            }

            vector dir = vector::zero;
            if (dirI == 0)
            {
                dir = vector(1, 0, 0);
            }
            else if (dirI == 1)
            {
                dir = vector(0, 1, 0);
            }
            else
            {
                dir = vector(0, 0, 1);
            }

            // Take a copy of matrixNoBCs and enforce the boundary conditions
            // for this direction
            sparseScalarMatrix matrixDirI(matrixNoBCs);

            // Take a copy of the source for this direction
            scalarField sourceDirI(source.component(dirI));

            // Enforce fixed DOF on the linear system for this direction
            // We are not current using fixedDofDirections_
            sparseMatrixTools::enforceFixedDof
            (
                matrixDirI,
                sourceDirI,
                fixedDofs_,
                dir & (dir & fixedDofDirections_),
                fixedDofValues_.component(dirI),
                fixedDofScale_,
                debug
            );

            // Solve linear system for displacement component correction

            // For now, use Eigen as linear solver
            // We can add PETSc or other aproaches later
            {
                // Lookup exportToMatlab flag
                const Switch writeMatlabMatrix
                (
                    solidModelDict().lookup("writeMatlabMatrix")
                );

                // Use Eigen SparseLU direct solver
                sparseMatrixTools::solveLinearSystemEigen
                (
                    matrixDirI, sourceDirI, pointDcorr, writeMatlabMatrix, debug
                );
            }

            pointD().primitiveFieldRef().replace
            (
                dirI, pointD().primitiveField().component(dirI) + pointDcorr
            );

            pointD().correctBoundaryConditions();
        }

        // Update point accelerations
        // Note: for NewmarkBeta, this needs to come before the pointU update
        pointA_.primitiveFieldRef() =
            vfvc::ddt
            (
                mesh().ddtScheme("ddt(pointU)"),
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                pointU_
            );

        // Update point velocities
        pointU_.primitiveFieldRef() =
            vfvc::ddt
            (
                mesh().ddtScheme("ddt(pointD)"),
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                pointD()
            );

        if (twoD_)
        {
            twoDCorrector_.correctPoints(pointD());

            // Remove displacement in the empty directions
            forAll(mesh().geometricD(), dirI)
            {
                if (mesh().geometricD()[dirI] < 0)
                {
                    pointD().primitiveFieldRef().replace(dirI, 0.0);
                }
            }
        }

        // Relax pointD
        if (mesh().relaxField(pointD().name()))
        {
            pointD().relax(mesh().fieldRelaxationFactor(pointD().name()));
        }

        // Update correction vector field
        pointDcorrVec = pointD() - pointD().prevIter();
    }
    while
    (
        !converged
        (
            iCorr,
            initResidual,
            mag(solverPerf.finalResidual()),
            cmptMax(solverPerf.nIterations()),
            pointD(),
            pointDcorrVec
        )
     && ++iCorr < nCorr()
    );

    // Calculate gradD at dual faces
    dualGradDf_ = vfvc::fGrad
    (
        pointD(),
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta_,
        debug
    );

    // Update the increment of displacement
    pointDD() = pointD() - pointD().oldTime();

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1)
    {
        mechanical().mapGradToSubMeshes(gradD());
    }

    // Update dual face stress field
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate
    mechanical().correct(sigma());

    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());

    return true;
}


bool vertexCentredLinGeomSolid::evolveExplicit()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the solid momentum equation for pointD" << nl
            << "Simulation Time, Clock Time, Max Stress" << endl;
    }

    physicsModel::printInfo() = bool
    (
        time().timeIndex() % infoFrequency() == 0
     || mag(time().value() - time().endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< time().value() << " " << time().elapsedClockTime()
            << " " << max(mag(dualSigmaf_)).value() << endl;

        physicsModel::printInfo() = false;
    }

    // Central difference scheme

    // Take a reference to the current and previous time-step
    const dimensionedScalar& deltaT = time().deltaT();
    //const dimensionedScalar& deltaT0 = time().deltaT0();

    // Compute the velocity
    // Note: this is the velocity at the middle of the time-step
    //pointU_ = pointU_.oldTime() + 0.5*(deltaT + deltaT0)*pointA_.oldTime();
    pointU_ = pointU_.oldTime() + deltaT*pointA_.oldTime();

    // Compute displacement
    pointD() = pointD().oldTime() + deltaT*pointU_;

    // Enforce boundary conditions on the displacement field
    pointD().correctBoundaryConditions();

    if (twoD_)
    {
        twoDCorrector_.correctPoints(pointD());

        // Remove displacement in the empty directions
        forAll(mesh().geometricD(), dirI)
        {
            if (mesh().geometricD()[dirI] < 0)
            {
                pointD().primitiveFieldRef().replace(dirI, 0.0);
            }
        }
    }

    // Update the divergence of stress based on the latest pointD field
    updatePointDivSigma(pointD(), dualGradDf_, dualSigmaf_, pointDivSigma_);

    // Compute acceleration
    pointA_ = pointDivSigma_/pointRho_ - dampingCoeff()*pointU_ + g();

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredLinGeomSolid::vertexCentredLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    dualMechanicalPtr_
    (
        new dualMechanicalModel
        (
            dualMesh(),
            nonLinGeom(),
            incremental(),
            mechanical(),
            dualMeshMap().dualFaceToCell()
        )
    ),
    dualImpKfPtr_(),
    zeta_(solidModelDict().lookupOrDefault<scalar>("zeta", 0.0)),
    fullNewton_
    (
        (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
      ? Switch(solidModelDict().lookup("fullNewton"))
      : Switch(false)
    ),
    steadyState_(false),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    twoDCorrector_(mesh()),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    fixedDofs_(mesh().nPoints(), false),
    fixedDofValues_(fixedDofs_.size(), vector::zero),
    fixedDofDirections_(fixedDofs_.size(), symmTensor::zero),
    fixedDofScale_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "fixedDofScale",
            (
                average(mechanical().impK())
               *Foam::sqrt(gAverage(mesh().magSf()))
            ).value()
        )
    ),
    pointU_
    (
        IOobject
        (
            "pointU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pointA_
    (
        IOobject
        (
            "pointA",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity/dimTime, vector::zero)
    ),
    pointRho_
    (
        IOobject
        (
            "point(rho)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    ),
    pointVol_
    (
        IOobject
        (
            "pointVolumes",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    pointGlobalVol_
    (
        IOobject
        (
            "pointGlobalVolumes",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    pointDivSigma_
    (
        IOobject
        (
            "pointDivSigma",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimForce/dimVolume, vector::zero)
    ),
    dualGradDf_
    (
        IOobject
        (
            "grad(D)f",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualSigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        "calculated"
    ),
    globalPointIndices_(mesh()),
    pointVolInterp_(pMesh(), mesh())
{
    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // pointD field must be defined
    pointDisRequired();

    // Set fixed degree of freedom list
    setFixedDofs(pointD(), fixedDofs_, fixedDofValues_, fixedDofDirections_);

    // Set point density field
    mechanical().volToPoint().interpolate(rho(), pointRho_);

    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
    scalarField& pointVolI = pointVol_.primitiveFieldRef();
    scalarField& pointGlobalVolI = pointGlobalVol_;
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
        pointGlobalVolI[pointID] = dualCellVol[dualCellI];
    }

    // Sum the shared point volumes to create the point global volumes
    pointConstraints::syncUntransformedData
    (
        mesh(), pointGlobalVol_, plusEqOp<scalar>()
    );

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointU_.oldTime().storeOldTime();
    pointA_.storeOldTime();

    // Write fixed degree of freedom equation scale
    Info<< "fixedDofScale: " << fixedDofScale_ << endl;

    // Disable the writing of the unused fields
    D().writeOpt() = IOobject::NO_WRITE;
    D().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    DD().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    U().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

vertexCentredLinGeomSolid::~vertexCentredLinGeomSolid()
{
#ifdef USE_PETSC
    if
    (
        solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED
     && Switch(solidModelDict().lookup("usePETSc"))
    )
    {
        PetscFinalize();
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> vertexCentredLinGeomSolid::residualMomentum
(
    const pointVectorField& pointD
)
{
    if (debug)
    {
        Info<< "void vertexCentredLinGeomSolid::residualMomentum(...): start"
            << endl;
    }

    // Prepare result
    tmp<vectorField> tresidual(new vectorField(pointD.size(), vector::zero));
    vectorField& residual = tresidual.ref();

    // The residual vector is defined as
    // F = div(sigma) + rho*g - rho*d2dt2(D)

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Point density field
    const scalarField& pointRhoI = pointRho_.internalField();

    // Update the displacement gradient at the dual faces
    // const surfaceTensorField dualGradDf
    // NOTE: dualMechanicalPtr_().correct currently looks up dualGradDf from
    // the object registry by name so we must use this field
    dualGradDf_ =
    // (
        vfvc::fGrad
        (
            pointD,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta_,
            debug
        );
    // );

    // Calculate the stress at the dual faces
    surfaceSymmTensorField dualSigmaf
    (
        IOobject
        (
            "sigmaf",
            runTime().timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        "calculated"
    );
    dualMechanicalPtr_().correct(dualSigmaf);

    // Calculate the tractions on the dual faces
    surfaceVectorField dualTraction
    (
        (dualMesh().Sf()/dualMesh().magSf()) & dualSigmaf
    );

    // Enforce extract tractions on traction boundaries
    enforceTractionBoundaries
    (
        pointD, dualTraction, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Calculate divergence of stress for the dual cells
    // const vectorField dualDivSigma = fvc::div(dualMesh().Sf() & dualSigmaf_);
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh().magSf());

    // Map dual cell field to primary mesh point field
    // vectorField& pointDivSigma = pointDivSigma_;
    // pointDivSigma = vector::zero;
    vectorField pointDivSigma(pointD.size(), vector::zero);
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigma[pointID] = dualDivSigma[dualCellI];
    }

    // Add surface forces to source
    residual += pointDivSigma*pointVolI;

    // Add gravity body forces
    residual += pointRhoI*g().value()*pointVolI;

    // Add transient term
    residual -= vfvc::d2dt2
    (
        mesh().d2dt2Scheme("d2dt2(pointD)"),
        pointD,
        pointU_,
        pointA_,
        pointRho_,
        pointVol_,
        int(bool(debug))
    );

    // Enforce fixed DOFs
    // This is only important if this residual function is used to approximate the Jacobian using finite differences
    const vectorField& pointDI = pointD;
    forAll(residual, pointI)
    {
        if (fixedDofs_[pointI])
        {
            // Free direction
            const symmTensor& fixedDir = fixedDofDirections_[pointI];
            const symmTensor freeDir = I - fixedDir;

            residual[pointI] =
                (fixedDir & (-pointDI[pointI])) + (freeDir & residual[pointI]);
        }
    }


    if (debug)
    {
        Info<< "void vertexCentredLinGeomSolid::residualMomentum(...): end"
            << endl;
    }

    return tresidual;
}


tmp<sparseMatrix> vertexCentredLinGeomSolid::JacobianMomentum
(
    const pointVectorField& pointD
)
{
    // Initialise matrix
    tmp<sparseMatrix> tmatrix
    (
        new sparseMatrix(sum(globalPointIndices_.stencilSize()))
    );
    sparseMatrix& matrix = tmatrix.ref();

    // Update gradD at dual faces
    dualGradDf_ =
        vfvc::fGrad
        (
            pointD,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta_,
            debug
        );

    // Calculate the stress at the dual faces
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Update material tangent
    const Field<scalarSquareMatrix> materialTangent
    (
        dualMechanicalPtr_().materialTangentFaceField()
    );

    // Add div(sigma) coefficients
    vfvm::divSigma
    (
        matrix,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        materialTangent,
        zeta_
    );

    // Add d2dt2 coefficients
    vfvm::d2dt2
    (
        mesh().d2dt2Scheme("d2dt2(pointD)"),
        runTime().deltaTValue(),
        pointD.name(),
        matrix,
        pointRho_.internalField(),
        pointVol_.internalField(),
        int(bool(debug))
    );

    // Enforce fixed DOF on the linear system
    vectorField source(pointD.size()); // needed but not used
    sparseMatrixTools::enforceFixedDof
    (
        matrix,
        source, // not used
        fixedDofs_,
        fixedDofDirections_,
        fixedDofValues_,
        fixedDofScale_
    );

    return tmatrix;
}


void vertexCentredLinGeomSolid::setDeltaT(Time& runTime)
{
    if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        // Max wave speed in the domain
        const scalar waveSpeed = max
        (
            Foam::sqrt(mechanical().impK()/mechanical().rho())
        ).value();

        // deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
        // In the current discretisation, information can move two cells per
        // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
        // calculating the required stable time-step
        // i.e. deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
        // For safety, we should use a time-step smaller than this e.g. Abaqus uses
        // stableTimeStep/sqrt(2): we will default to this value
        const scalar requiredDeltaT =
            1.0/
            gMax
            (
                DimensionedField<scalar, Foam::surfaceMesh>
                (
                    dualMesh().surfaceInterpolation::
                    deltaCoeffs().internalField()
                   *waveSpeed
                )
            );

        // Lookup the desired Courant number
        const scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

        const scalar newDeltaT = maxCo*requiredDeltaT;

        // Update print info
        physicsModel::printInfo() = bool
        (
            runTime.timeIndex() % infoFrequency() == 0
         || mag(runTime.value() - runTime.endTime().value()) < SMALL
        );

        physicsModel::printInfo() = false;

        if (time().timeIndex() == 1)
        {
            Info<< nl << "Setting deltaT = " << newDeltaT
                << ", maxCo = " << maxCo << endl;
        }

        runTime.setDeltaT(newDeltaT);
    }
}


bool vertexCentredLinGeomSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    {
        return evolveImplicitCoupled();
    }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            << solutionAlgorithmNames_.names() << endl;
    }

    // Keep compiler happy
    return true;
}


void vertexCentredLinGeomSolid::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    // Get point field on patch
    const vectorField traction
    (
        globalPatches()[interfaceI].globalPointToPatch
        (
            globalPatches()[interfaceI].interpolator().faceToPointInterpolate
            (
                faceZoneTraction
            )()
        )
    );

    // Lookup point patch field
    pointPatchVectorField& ptPatch = pointD().boundaryFieldRef()[patchID];

    if (isA<solidTractionPointPatchVectorField>(ptPatch))
    {
        solidTractionPointPatchVectorField& patchD =
            refCast<solidTractionPointPatchVectorField>(ptPatch);

        patchD.traction() = traction;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::vertexCentredLinGeomSolid::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << ptPatch.type()
            << " for point patch " << ptPatch.patch().name()
            << " should instead be type "
            << solidTractionPointPatchVectorField::typeName
            << abort(FatalError);
    }
}


void vertexCentredLinGeomSolid::writeFields(const Time& runTime)
{
    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1)
    {
        mechanical().mapGradToSubMeshes(gradD());
    }

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate
    mechanical().correct(sigma());

    // Calculate gradD at the primary points using least squares: this should
    // be second-order accurate (... I think).
    const pointTensorField pGradD(vfvc::pGrad(pointD(), mesh()));

    // Calculate strain at the primary points based on pGradD
    // Note: the symm operator is not defined for pointTensorFields so we will
    // do it manually
    // const pointSymmTensorField pEpsilon("pEpsilon", symm(pGradD));
    pointSymmTensorField pEpsilon
    (
        IOobject
        (
            "pEpsilon",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedSymmTensor("0", dimless, symmTensor::zero)
    );

    pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
    pEpsilon.write();

    // Equivalent strain at the points
    pointScalarField pEpsilonEq
    (
        IOobject
        (
            "pEpsilonEq",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimless, 0.0)
    );

    pEpsilonEq.primitiveFieldRef() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));

    pEpsilonEq.write();

    Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << endl;

    // Stress at the points
    pointSymmTensorField pSigma
    (
        IOobject
        (
            "pSigma",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    // Calculate the stress at the points
    mechanical().correct(pSigma, pGradD);

    // Write point sigma
    pSigma.write();

    // Equivalent stress at the points
    pointScalarField pSigmaEq
    (
        IOobject
        (
            "pSigmaEq",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //sqrt((3.0/2.0)*magSqr(dev(pSigma)))
        pMesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    );

    pSigmaEq.primitiveFieldRef() =
        sqrt((3.0/2.0)*magSqr(dev(pSigma.primitiveField())));

    pSigmaEq.write();

    solidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // OPENFOAM_COM

// ************************************************************************* //
