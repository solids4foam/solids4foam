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

#include "vertexCentredNonLinGeomTotalLagSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "sparseMatrixTools.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#ifdef USE_PETSC
    #include <petscksp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, vertexCentredNonLinGeomTotalLagSolid, dictionary
);

const Enum<vertexCentredNonLinGeomTotalLagSolid::solutionAlgorithm>
vertexCentredNonLinGeomTotalLagSolid::solutionAlgorithmNames_
({
    {
        vertexCentredNonLinGeomTotalLagSolid::solutionAlgorithm::IMPLICIT_COUPLED,
        "implicitCoupled"
    },
    {
        vertexCentredNonLinGeomTotalLagSolid::solutionAlgorithm::IMPLICIT_SEGREGATED,
        "implicitSegregated"
    },
    {
        vertexCentredNonLinGeomTotalLagSolid::solutionAlgorithm::EXPLICIT,
        "explicit"
    },
});


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void vertexCentredNonLinGeomTotalLagSolid::updateSource
(
    vectorField& source,
    const labelList& dualCellToPoint
)
{
    if (debug)
    {
        Info<< "void updateSource(...): start" << endl;
    }

    // Reset to zero
    source = vector::zero;

    // The source vector is -F, where:
    // F = div(JF^-T * sigma) + rho*g - rho*d2dt2(D)

    // Point density field
    const scalarField& pointRhoI = pointRho_.internalField();

    // Dual face area vectors at deformed configuration
    const surfaceVectorField deformedDualSf
    (
        dualJf_*dualFinvf_.T() & dualMesh().Sf()
    );

    // Magnitude of deformedDualSf
    const surfaceScalarField deformedDualMagSf(mag(deformedDualSf));

    // Dual face unit normals at deformed configuration
    const surfaceVectorField deformedDualN(deformedDualSf/deformedDualMagSf);

    // Calculate the Cauchy tractions on the dual faces
    surfaceVectorField dualTraction(deformedDualN & dualSigmaf_);

    // Enforce extract tractions on traction boundaries
    enforceTractionBoundaries
    (
        pointD(),
        dualTraction,
        deformedDualN,
        mesh(),
        dualMeshMap().pointToDualFaces()
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

    // Calculate divergence of stress for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*deformedDualMagSf);

    // // Calculate absolute divergence of stress (force)
    // // We do this to allow syncing of forces at points on processor boundaries
    // const vectorField dualDivSigmaAbs(dualDivSigma*dualMesh().V());

    // // Map dual cell field to primary mesh point field
    // // We temporarily use the pointDivSigma field to hold absolute forces
    // // but convert them back to force per unit volume below
    vectorField& pointDivSigmaI = pointDivSigma_;
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
	pointDivSigmaI[pointID] = dualDivSigma[dualCellI];
    }

    // forAll(dualDivSigmaAbs, dualCellI)
    // {
    //     const label pointID = dualCellToPoint[dualCellI];
    //     pointDivSigmaI[pointID] = dualDivSigmaAbs[dualCellI];
    // }

    // // Sum absolute forces in parallel
    // pointConstraints::syncUntransformedData
    // (
    //     mesh(), pointDivSigma_, plusEqOp<vector>()
    // );

    // // Convert force to force per unit volume
    // // Perform calculation per point to avoid dimension checks
    // // Point volume field
    // const scalarField& pointVolI = pointVol_.internalField();
    // forAll(pointDivSigmaI, pointI)
    // {
    //     pointDivSigmaI[pointI] /= pointVolI[pointI];
    // }

    const scalarField& pointVolI = pointVol_.internalField();

    // Add surface forces to source
    source -= pointDivSigmaI*pointVolI;

    // Add gravity body forces
    source -= pointRhoI*g().value()*pointVolI;

    // Add transient term
    source += vfvc::d2dt2
    (
        mesh().d2dt2Scheme("d2dt2(pointD)"),
        pointD(),
        pointU_,
        pointA_,
        pointRho_,
        pointVol_,
        int(bool(debug))
    );

    if (debug)
    {
        Info<< "void vertexCentredNonLinGeomTotalLagSolid::"
            << "updateSource(...): end"
            << endl;
    }
}


void vertexCentredNonLinGeomTotalLagSolid::updatePointDivSigma
(
    const pointVectorField& pointD,
    surfaceTensorField& dualGradDf,
    surfaceTensorField& dualFf,
    surfaceTensorField& dualFinvf,
    surfaceScalarField& dualJf,
    surfaceSymmTensorField& dualSigmaf,
    pointVectorField& pointDivSigma
)
{
    if (debug)
    {
        Info<< "void updatePointDivSigma(...): start" << endl;
    }

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.1));

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

    // Update F
    dualFf = I + dualGradDf.T();

    // Update Finv
    dualFinvf = inv(dualFf);

    // Update J
    dualJf = det(dualFf);

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Dual face area vectors at deformed configuration
    const surfaceVectorField deformedDualSf
    (
        dualJf*dualFinvf.T() & dualMesh().Sf()
    );

    // Magnitude of deformedDualSf
    const surfaceScalarField deformedDualMagSf(mag(deformedDualSf));

    // Dual face unit normals at deformed configuration
    const surfaceVectorField deformedDualN(deformedDualSf/deformedDualMagSf);

    // Calculate the Cauchy tractions on the dual faces
    surfaceVectorField dualTraction(deformedDualN & dualSigmaf);

    // Enforce extract tractions on traction boundaries
    enforceTractionBoundaries
    (
        pointD,
        dualTraction,
        deformedDualN,
        mesh(),
        dualMeshMap().pointToDualFaces()
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
    const vectorField dualDivSigma = fvc::div(dualTraction*deformedDualMagSf);

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
    // forAll(dualDivSigma, dualCellI)
    // {
    //     const label pointID = dualCellToPoint[dualCellI];
    //     pointDivSigma[pointID] = dualDivSigma[dualCellI];
    // }

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


void vertexCentredNonLinGeomTotalLagSolid::setFixedDofs
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
                            "void vertexCentredNonLinGeomTotalLagSolid::setFixedDofs(...)"
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
                            "void vertexCentredNonLinGeomTotalLagSolid::setFixedDofs(...)"
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


void vertexCentredNonLinGeomTotalLagSolid::enforceTractionBoundaries
(
    const pointVectorField& pointD,
    surfaceVectorField& dualTraction,
    const surfaceVectorField& dualDeformedNormals,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
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

            // Primary mesh point tractions
            const vectorField& pointTraction = tracPatch.traction();
            const scalarField& pointPressure = tracPatch.pressure();

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
            forAll(pointTraction, pI)
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

                            // Dual face deformed unit normal
                            const vector& n =
                                dualDeformedNormals.boundaryField()
                                [
                                    patchI
                                ][localDualFaceID];

                            // Set dual face traction
                            // Use the deformed unit normal for this face for
                            // the pressure
                            dualFaceTraction[localDualFaceID] +=
                                pointTraction[pI] - n*pointPressure[pI];

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
                    "void vertexCentredNonLinGeomTotalLagSolid::"
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
            // It is assumed that the deformedN is the same as the initial
            // reference normal
            const vectorField n(dualMesh.boundary()[patchI].nf());
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
        }
    }
}

bool vertexCentredNonLinGeomTotalLagSolid::
vertexCentredNonLinGeomTotalLagSolid::converged
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
            Info<< "    Initial residual is less than 1e-15"
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
    else if (residualAbs < SMALL)
    {
        Info<< "    Converged: absolute residual is less than " << SMALL
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


scalar vertexCentredNonLinGeomTotalLagSolid::calculateLineSearchSlope
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

    // Update F
    dualFf_ = I + dualGradDf.T();

    // Update Finv
    dualFinvf_ = inv(dualFf_);

    // Update J
    dualJf_ = det(dualFf_);

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Update the source vector
    vectorField source(mesh().nPoints(), vector::zero);
    pointD.correctBoundaryConditions();
    updateSource(source, dualMeshMap().dualCellToPoint());

    // Reset pointD
    pointD = pointD.prevIter();

    // Return the slope
    return gSum(pointDcorr & source);
}


scalar vertexCentredNonLinGeomTotalLagSolid::calculateLineSearchFactor
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

Foam::tmp<Foam::Field<Foam::RectangularMatrix<Foam::scalar>>>
vertexCentredNonLinGeomTotalLagSolid::geometricStiffnessField
(
    const surfaceVectorField& SfUndef, //Undeformed surface area vector field
    const surfaceTensorField& gradDRef //Reference gradD
) const
{
    // Prepare tmp field
    // Todo: switch to 3 x 6 matrix instead of 3 x 9
    tmp<Field<RectangularMatrix<scalar>>> tresult
    (
        new Field<RectangularMatrix<scalar>>
        (
            dualMesh().nFaces(), RectangularMatrix<scalar>(3, 9, 0.0)
        )
    );
    Field<RectangularMatrix<scalar>>& result = tresult.ref();

    // For small strain the geometric stiffness is zero
    if (!geometricStiffness_)
    {
        return tresult;
    }

    // Calculate surface vector as per the total Lagrangian formulation
    // gamma = JF^-T*Sf0

    // Calculate unperturbed F
    const surfaceTensorField FRef(I + gradDRef.T());

    // Calculate unperturbed invF
    const surfaceTensorField invFRef(inv(FRef));

    // Calculate unperturbed J
    const surfaceScalarField JRef(det(FRef));

    // Calculate unperturbed deformed Sf
    const surfaceVectorField SfRef((JRef*invFRef.T()) & SfUndef);

    // Create field to be used for perturbations
    surfaceTensorField gradDPerturb("gradDPerturb", gradDRef);

    // Small number used for perturbations
    const scalar eps(solidModelDict().lookupOrDefault<scalar>("tangentEps", 1e-10));

    // For each component of gradD, sequentially apply a perturbation and
    // then calculate the resulting sigma
    for (label cmptI = 0; cmptI < tensor::nComponents; cmptI++)
    {

        // Reset gradDPerturb and multiply by 1.0 to avoid it being removed
        // from the object registry
        gradDPerturb = 1.0*gradDRef;

        // Perturb this component of gradD and calculate SfPerturb
        gradDPerturb.replace(cmptI, gradDRef.component(cmptI) + eps);

        //Info << gradDPerturb[0] - gradDRef[0] << endl;

        const surfaceTensorField FPerturb(I + gradDPerturb.T());
        const surfaceTensorField invFPerturb(inv(FPerturb));
        const surfaceScalarField JPerturb(det(FPerturb));
        const surfaceVectorField SfPerturb((JPerturb*invFPerturb.T()) & SfUndef);

        // Calculate each component
        const surfaceVectorField tangCmpt((SfPerturb - SfRef)/eps);
        const vectorField& tangCmptI = tangCmpt.internalField();

        // Insert each component
        forAll(tangCmptI, faceI)
        {
            if (cmptI == tensor::XX)
            {
                result[faceI](0,0) = tangCmptI[faceI][0];
                result[faceI](1,0) = tangCmptI[faceI][1];
                result[faceI](2,0) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::XY)
            {
                result[faceI](0,1) = tangCmptI[faceI][0];
                result[faceI](1,1) = tangCmptI[faceI][1];
                result[faceI](2,1) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::XZ)
            {
                result[faceI](0,2) = tangCmptI[faceI][0];
                result[faceI](1,2) = tangCmptI[faceI][1];
                result[faceI](2,2) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::YX)
            {
                result[faceI](0,3) = tangCmptI[faceI][0];
                result[faceI](1,3) = tangCmptI[faceI][1];
                result[faceI](2,3) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::YY)
            {
                result[faceI](0,4) = tangCmptI[faceI][0];
                result[faceI](1,4) = tangCmptI[faceI][1];
                result[faceI](2,4) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::YZ)
            {
                result[faceI](0,5) = tangCmptI[faceI][0];
                result[faceI](1,5) = tangCmptI[faceI][1];
                result[faceI](2,5) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::ZX)
            {
                result[faceI](0,6) = tangCmptI[faceI][0];
                result[faceI](1,6) = tangCmptI[faceI][1];
                result[faceI](2,6) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::ZY)
            {
                result[faceI](0,7) = tangCmptI[faceI][0];
                result[faceI](1,7) = tangCmptI[faceI][1];
                result[faceI](2,7) = tangCmptI[faceI][2];
            }
            else if (cmptI == tensor::ZZ)
            {
                result[faceI](0,8) = tangCmptI[faceI][0];
                result[faceI](1,8) = tangCmptI[faceI][1];
                result[faceI](2,8) = tangCmptI[faceI][2];
            }
        }

        forAll(tangCmpt.boundaryField(), patchI)
        {
            const vectorField& tangCmptP =
                tangCmpt.boundaryField()[patchI];
            const label start = mesh().boundaryMesh()[patchI].start();

            forAll(tangCmptP, fI)
            {
                const label faceID = start + fI;

                if (cmptI == tensor::XX)
                {
                    result[faceID](0,0) = tangCmptI[fI][0];
                    result[faceID](1,0) = tangCmptI[fI][1];
                    result[faceID](2,0) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::XY)
                {
                    result[faceID](0,1) = tangCmptI[fI][0];
                    result[faceID](1,1) = tangCmptI[fI][1];
                    result[faceID](2,1) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::XZ)
                {
                    result[faceID](0,2) = tangCmptI[fI][0];
                    result[faceID](1,2) = tangCmptI[fI][1];
                    result[faceID](2,2) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::YX)
                {
                    result[faceID](0,3) = tangCmptI[fI][0];
                    result[faceID](1,3) = tangCmptI[fI][1];
                    result[faceID](2,3) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::YY)
                {
                    result[faceID](0,4) = tangCmptI[fI][0];
                    result[faceID](1,4) = tangCmptI[fI][1];
                    result[faceID](2,4) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::YZ)
                {
                    result[faceID](0,5) = tangCmptI[fI][0];
                    result[faceID](1,5) = tangCmptI[fI][1];
                    result[faceID](2,5) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::ZX)
                {
                    result[faceID](0,6) = tangCmptI[fI][0];
                    result[faceID](1,6) = tangCmptI[fI][1];
                    result[faceID](2,6) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::ZY)
                {
                    result[faceID](0,7) = tangCmptI[fI][0];
                    result[faceID](1,7) = tangCmptI[fI][1];
                    result[faceID](2,7) = tangCmptI[fI][2];
                }
                else if (cmptI == tensor::ZZ)
                {
                    result[faceID](0,8) = tangCmptI[fI][0];
                    result[faceID](1,8) = tangCmptI[fI][1];
                    result[faceID](2,8) = tangCmptI[fI][2];
                }
            }
        }
    }

    return tresult;
}


void vertexCentredNonLinGeomTotalLagSolid::makeDualImpKf() const
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


const surfaceScalarField&
vertexCentredNonLinGeomTotalLagSolid::dualImpKf() const
{
    if (dualImpKfPtr_.empty())
    {
        makeDualImpKf();
    }

    return dualImpKfPtr_();
}


surfaceScalarField&
vertexCentredNonLinGeomTotalLagSolid::dualImpKf()
{
    if (dualImpKfPtr_.empty())
    {
        makeDualImpKf();
    }

    return dualImpKfPtr_();
}


void vertexCentredNonLinGeomTotalLagSolid::predict()
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


bool vertexCentredNonLinGeomTotalLagSolid::evolveImplicitCoupled()
{
    Info<< "Evolving solid solver" << endl;

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.1));
    if (debug)
    {
        Info<< "zeta: " << zeta << endl;
    }

    // Initialise matrix
    sparseMatrix matrix(sum(globalPointIndices_.stencilSize()));

    // Calculate F
    dualFf_ = I + dualGradDf_.T();

    // Calculate Finv
    dualFinvf_ = inv(dualFf_);

    // Calculate J
    dualJf_ = det(dualFf_);

    // Store the initial material tangent field for dual mesh faces
    Field<scalarSquareMatrix> materialTangent
    (
        dualMechanicalPtr_().materialTangentFaceField()
    );

    // Store the dual undeformed area vectors
    const surfaceVectorField& dualSf = dualMesh().Sf();

    // Calculate the initial geometric stiffness field for dual mesh faces
    Field<RectangularMatrix<scalar>> geometricStiffness
    (
        geometricStiffnessField(dualSf, dualGradDf_)
    );

    // Calculate stress field at dual faces
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Calculate stress for primary cells
    mechanical().correct(sigma());

    // Global point index lists
    const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    const labelList& localToGlobalPointMap =
        globalPointIndices_.localToGlobalPointMap();

    if (!fullNewton_)
    {
        // Assemble matrix once per time-step
        Info<< "    Assembling the matrix" << endl;

        // Add div(sigma) coefficients
        vfvm::divSigma
        (
            matrix,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            materialTangent,
            geometricStiffness,
            dualSigmaf_,
            dualGradDf_,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofScale_,
            zeta,
            debug
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

    // Solution field: point displacement correction
    vectorField pointDcorr(pointD().internalField().size(), vector::zero);

    // Initialise the source
    vectorField source(mesh().nPoints(), vector::zero);

    // Newton-Raphson loop over momentum equation
    int iCorr = 0;
    scalar initResidual = 0.0;
    SolverPerformance<vector> solverPerf;
    do
    {
        // Calculate gradD at dual faces
        dualGradDf_ = vfvc::fGrad
        (
            pointD(),
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta,
            debug
        );

        // Update F
        dualFf_ = I + dualGradDf_.T();

        // Update Finv
        dualFinvf_ = inv(dualFf_);

        // Update J
        dualJf_ = det(dualFf_);

        // Calculate stress at dual faces
        dualMechanicalPtr_().correct(dualSigmaf_);

        // Update the source vector
        pointD().correctBoundaryConditions();
        updateSource(source, dualMeshMap().dualCellToPoint());

        if (fullNewton_)
        {
            // Assemble the matrix once per outer iteration
            matrix.clear();

            // Update the material tangent
            materialTangent = dualMechanicalPtr_().materialTangentFaceField();

            // Update the geometric stiffness
            geometricStiffness = geometricStiffnessField(dualSf, dualGradDf_);

            // Add div(sigma) coefficients
            vfvm::divSigma
            (
                matrix,
                mesh(),
                dualMesh(),
                dualMeshMap().dualFaceToCell(),
                dualMeshMap().dualCellToPoint(),
                materialTangent,
                geometricStiffness,
                dualSigmaf_,
                dualGradDf_,
                fixedDofs_,
                fixedDofDirections_,
                fixedDofScale_,
                zeta,
                debug
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

        if (debug)
        {
            Info<< "bool evolve(): solving linear system: start" << endl;
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
                ownedByThisProc,
                localToGlobalPointMap,
                globalPointIndices_.stencilSizeOwned(),
                globalPointIndices_.stencilSizeNotOwned(),
                solidModelDict().lookupOrDefault<bool>("debugPETSc", false)
            );
#else
            FatalErrorIn
            (
                "vertexCentredNonLinGeomTotalLagSolid::evolve()"
            )   << "PETSc not available. Please set the PETSC_DIR environment "
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
            Info<< "bool vertexCentredNonLinGeomTotalLagSolid::evolve(): "
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
                    rTol, maxIter, pointDcorr, source, zeta
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
        zeta,
        debug
    );

    // Update F
    dualFf_ = I + dualGradDf_.T();

    // Update Finv
    dualFinvf_ = inv(dualFf_);

    // Update J
    dualJf_ = det(dualFf_);

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
    // In contrast, pSigma will be 2nd order
    mechanical().correct(sigma());

    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());

    return true;
}


bool vertexCentredNonLinGeomTotalLagSolid::evolveImplicitSegregated()
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
    if (debug)
    {
        Info<< "compactImplicitStencil: " << compactImplicitStencil << endl;
    }

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

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.1));
    if (debug)
    {
        Info<< "zeta: " << zeta << endl;
    }

#ifdef USE_PETSC
    // Global point index lists
    const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    const labelList& localToGlobalPointMap =
        globalPointIndices_.localToGlobalPointMap();
#endif

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
    SolverPerformance<scalar> solverScalarPerf;
    do
    {
        // Store previous iteration of pointD for residual calculation
        pointD().storePrevIter();

        // Calculate gradD at dual faces
        dualGradDf_ = vfvc::fGrad
        (
            pointD(),
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta,
            debug
        );

        // Update F
        dualFf_ = I + dualGradDf_.T();

        // Update Finv
        dualFinvf_ = inv(dualFf_);

        // Update J
        dualJf_ = det(dualFf_);

        // Calculate stress at dual faces
        dualMechanicalPtr_().correct(dualSigmaf_);

        // Update the source vector
        pointD().correctBoundaryConditions();
        updateSource(source, dualMeshMap().dualCellToPoint());

        // Update Laplacian diffusivity
        if (fullNewton_)
        {
            // Assemble the matrix once per outer iteration
            matrixNoBCs.clear();

            // Calculate the full material tangent
            Field<scalarSquareMatrix> materialTangent
            (
                dualMechanicalPtr_().materialTangentFaceField()
            );

            // Update impK
            scalarField& impKfI = dualImpKf().primitiveFieldRef();
            const label XX = symmTensor::XX;
            const label YY = symmTensor::YY;
            const label ZZ = symmTensor::ZZ;
            forAll(impKfI, dualFaceI)
            {
                // Todo: try XX coeffs for X direction, etc.
                // Todo: lookup update frequency
                impKfI[dualFaceI] =
                    (
                        materialTangent[dualFaceI](XX, XX)
                      + materialTangent[dualFaceI](YY, YY)
                      + materialTangent[dualFaceI](ZZ, ZZ)
                    )/3.0;
                    // min
                    // (
                    //     materialTangent[dualFaceI](XX, XX),
                    //     min
                    //     (
                    //         materialTangent[dualFaceI](YY, YY),
                    //         materialTangent[dualFaceI](ZZ, ZZ)
                    //     )
                    // );
            }

            // Update the matrix
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
        }

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
            // We can add PETSc or other aproaches (OpenFOAM solvers?) later
            if (Switch(solidModelDict().lookup("usePETSc")))
            {
#ifdef USE_PETSC
                fileName optionsFile(solidModelDict().lookup("optionsFile"));
                solverScalarPerf = sparseMatrixTools::solveLinearSystemPETSc
                (
                    matrixDirI,
                    sourceDirI,
                    pointDcorr,
                    twoD_,
                    optionsFile,
                    mesh().points(),
                    ownedByThisProc,
                    localToGlobalPointMap,
                    globalPointIndices_.stencilSizeOwned(),
                    globalPointIndices_.stencilSizeNotOwned(),
                    solidModelDict().lookupOrDefault<bool>("debugPETSc", false)
                );

                solverPerf.finalResidual() =
                    solverScalarPerf.nIterations()*vector::one;
                solverPerf.nIterations() =
                    solverScalarPerf.nIterations()*vector::one;
#else
                FatalErrorIn("evolveImplicitSegregated()")
                    << "PETSc not available. Please set the PETSC_DIR environment "
                    << "variable and re-compile solids4foam" << abort(FatalError);
#endif
            }
            else
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
        zeta,
        debug
    );

    // Update F
    dualFf_ = I + dualGradDf_.T();

    // Update Finv
    dualFinvf_ = inv(dualFf_);

    // Update J
    dualJf_ = det(dualFf_);

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
    // In contrast, pSigma will be 2nd order
    mechanical().correct(sigma());

    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());

    return true;
}


bool vertexCentredNonLinGeomTotalLagSolid::evolveExplicit()
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
    updatePointDivSigma
    (
        pointD(),
        dualGradDf_,
        dualFf_,
        dualFinvf_,
        dualJf_,
        dualSigmaf_,
        pointDivSigma_
    );

    // Compute acceleration
    pointA_ = pointDivSigma_/pointRho_ - dampingCoeff()*pointU_ + g();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredNonLinGeomTotalLagSolid::
vertexCentredNonLinGeomTotalLagSolid
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
    solutionAlgorithm_
    (
        solutionAlgorithmNames_.get("solutionAlgorithm", solidModelDict())
    ),
    dualImpKfPtr_(),
    fullNewton_
    (
        (
            solutionAlgorithm_ == solutionAlgorithm::IMPLICIT_COUPLED
         || solutionAlgorithm_ == solutionAlgorithm::IMPLICIT_SEGREGATED
        )
      ? Switch(solidModelDict().lookup("fullNewton"))
      : Switch(false)
    ),
    geometricStiffness_
    (
        (solutionAlgorithm_ == solutionAlgorithm::IMPLICIT_COUPLED)
      ? Switch(solidModelDict().lookup("geometricStiffness"))
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
    pointP_
    (
        IOobject
        (
            "pointP",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimPressure, 0.0)
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
            IOobject::READ_IF_PRESENT,
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
    dualFf_
    (
        IOobject
        (
            "Ff",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualFinvf_
    (
        IOobject
        (
            "Finvf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualJf_
    (
        IOobject
        (
            "Jf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedScalar("zero", dimless, 1.0),
        "calculated"
    ),
    dualPf_
    (
        IOobject
        (
            "pf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        "calculated"
    ),
        volP_
        (
                IOobject
                (
                    "volP",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
                ),
                mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
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

    // Set the pointVol and pointGlobalVol fields
    // Map dualMesh cell volumes to the primary mesh points
    scalarField& pointVolI = pointVol_;
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

vertexCentredNonLinGeomTotalLagSolid::~vertexCentredNonLinGeomTotalLagSolid()
{
#ifdef USE_PETSC
    if (Switch(solidModelDict().lookup("usePETSc")))
    {
        PetscFinalize();
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vertexCentredNonLinGeomTotalLagSolid::setDeltaT(Time& runTime)
{
    if (solutionAlgorithm_ == solutionAlgorithm::EXPLICIT)
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
                    dualMesh().surfaceInterpolation::deltaCoeffs().internalField()
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


bool vertexCentredNonLinGeomTotalLagSolid::evolve()
{
    if (solutionAlgorithm_ == solutionAlgorithm::IMPLICIT_COUPLED)
    {
        return evolveImplicitCoupled();
    }
    else if (solutionAlgorithm_ == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    else if (solutionAlgorithm_ == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
    else
    {
        FatalErrorIn("bool evolve()")
            << "Unrecognised solution algorithm. Available options are "
            << solutionAlgorithmNames_.names() << endl;
    }

    // Keep compiler happy
    return true;
}


void vertexCentredNonLinGeomTotalLagSolid::setTraction
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
            "void Foam::vertexCentredNonLinGeomTotalLagSolid::setTraction\n"
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

void vertexCentredNonLinGeomTotalLagSolid::writeFields(const Time& runTime)
{
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

    solidModel::writeFields(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // OPENFOAM_COM

// ************************************************************************* //
