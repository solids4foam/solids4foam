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

#include "vertexCentredNonLinGeomTotalLagPressureDisplacementSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrixExtended.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPointExtended.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "sparseMatrixExtendedTools.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#include "neoHookeanElasticMisesPlastic.H"
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

defineTypeNameAndDebug(vertexCentredNonLinGeomTotalLagPressureDisplacementSolid, 0);
addToRunTimeSelectionTable(solidModel, vertexCentredNonLinGeomTotalLagPressureDisplacementSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::updateSource
(
    Field<RectangularMatrix<scalar>>& source,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const scalar& pressureSmoothing,
    const scalar& zeta,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::updateSource(...): start"
            << endl;
    }

    // Reset to zero
    source = RectangularMatrix<scalar>(4,1,0);

    // The source vector for the momentum equation is -F, where:
    // F = div(JF^-T * sigma) + rho*g - rho*d2dt2(D)

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

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
        pointD(), dualTraction, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAM_NOT_EXTEND
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
#else
            dualTraction.boundaryField()[patchI] = vector::zero;
#endif
        }
    }

    // Calculate divergence of stress for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*deformedDualMagSf);

    // Map dual cell fields to primary mesh point fields
    vectorField pointDivSigmaI(mesh().nPoints(), vector::zero);
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigmaI[pointID] = dualDivSigma[dualCellI];       	 
    }
    
    // Insert momentum coefficients into the source

    forAll (source, pointI)
    {
        for (int i = 0; i < 3; i++)
        {
            // Add surface forces
            source[pointI](i,0) -= pointDivSigmaI[pointI].component(i)*pointVolI[pointI];

            // Add gravity body forces
            source[pointI](i,0) -= pointRhoI[pointI]*g().value().component(i)*pointVolI[pointI];
        }        
    }   

    // The source vector for the pressure equation is -F, where:
    // F = p - gamma*laplacian(p) - pBar(D)

    // Calculate laplacian(p) field
    const pointScalarField laplacianP
    (
        vfvc::laplacian
        (
            pointP_,
            mesh(),
            dualMesh(),
            dualFaceToCell,
            dualCellToPoint,
            zeta,
            debug
        )
    );

    //Calculate the bulk modulus (TO FIX LATER)
    const scalar E = 70e9;
    const scalar nu = 0.3;
    const scalar K( (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*(E/(2.0*(1.0 + nu))) );
    
    // Calculate the pBar field (TO FIX LATER)
    scalarField pBar(-0.5*K*pow(pointJ_, 2.0));    
    pBar += 0.5*K;
    pBar /= pointJ_;

    // Insert pressure coefficients into the source
    forAll (source, pointI)
    {
        // Add pressure term
        source[pointI](3,0) -= pointP_[pointI]*pointVolI[pointI];

        // Add gamma*laplacian(p) term
        source[pointI](3,0) += pressureSmoothing*laplacianP[pointI]*pointVolI[pointI];

        // Add pBar term
        source[pointI](3,0) += pBar[pointI]*pointVolI[pointI];
    }

    if (debug)
    {
        Info<< "void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::updateSource(...): end"
            << endl;
    }
}


void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::setFixedDofs
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
                            "void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::setFixedDofs(...)"
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
                            "void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::setFixedDofs(...)"
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


void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::enforceTractionBoundaries
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

            // Todo: use deformedN

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
                    "void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;

            // Overwrite the dual patch face traction
#ifdef OPENFOAM_NOT_EXTEND
            dualTraction.boundaryFieldRef()[patchI] = dualFaceTraction;
#else
            dualTraction.boundaryField()[patchI] = dualFaceTraction;
#endif
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

            // Todo: use deformedN

            const vectorField n(dualMesh.boundary()[patchI].nf());
#ifdef OPENFOAM_NOT_EXTEND
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#else
            dualTraction.boundaryField()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#endif
        }
    }
}


bool vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::converged
(
    const label iCorr,
    scalar& initResidualD,
    scalar& initResidualP,
    const label nInterations,
    const pointVectorField& pointD,
    const pointScalarField& pointP,
    const Field<RectangularMatrix<scalar>>& pointDPcorr
) const
{
    scalar residualDAbs = 0;
    scalar residualPAbs = 0;
    scalar nPoints = 0;
    // Calculate the residuals as the root mean square of the correction
    // Displacement residual
    forAll(pointDPcorr, pointI)
    {
        // Displacement residual
        residualDAbs += sqr(pointDPcorr[pointI](0,0)) + sqr(pointDPcorr[pointI](1,0)) +
            sqr(pointDPcorr[pointI](2,0));

        // Pressure residual
        residualPAbs += sqr(pointDPcorr[pointI](3,0));
        
        nPoints ++;
    }
    
    residualDAbs /= sqrt(residualDAbs/nPoints);
    residualPAbs /= sqrt(residualPAbs/nPoints);
    
    // Store initial residual
    if (iCorr == 0)
    {
        initResidualD = residualDAbs;
        initResidualP = residualPAbs;

        // If the initial residual is small then convergence has been achieved
        if (initResidualD < SMALL && initResidualP < SMALL)
        {
            Info<< "    Both displacement and pressure residuals are less than 1e-15"
                << "    Converged" << endl;
            return true;
        }
        Info<< "    Initial displacement residual = " << initResidualD << endl;
        Info<< "    Initial pressure residual = " << initResidualP << endl;
    }

    // Define a normalised residual wrt the initial residual
    const scalar residualDNorm = residualDAbs/initResidualD;
    const scalar residualPNorm = residualPAbs/initResidualP;

    // Calculate the maximum displacement
#ifdef OPENFOAM_NOT_EXTEND
    const scalar maxMagD = gMax(mag(pointD.primitiveField()));
#else
    const scalar maxMagD = gMax(mag(pointD.internalField()));
#endif

    // Calculate the maximum pressure
#ifdef OPENFOAM_NOT_EXTEND
    const scalar maxMagP = gMax(mag(pointP.primitiveField()));
#else
    const scalar maxMagP = gMax(mag(pointP.internalField()));
#endif

    // Print information for the displacement
    Info<< "    Displacement residuals: " << endl
        << "    Iter = " << iCorr
        << ", relRef = " << residualDNorm
        << ", resAbs = " << residualDAbs
        << ", nIters = " << nInterations
        << ", maxD = " << maxMagD << endl;

    // Print information for the pressure
    Info<< "    Pressure residuals: " << endl
        << "    Iter = " << iCorr
        << ", relRef = " << residualPNorm
        << ", resAbs = " << residualPAbs
        << ", nIters = " << nInterations
        << ", maxP = " << maxMagP << endl;

    // Displacement tolerance
    const scalar DTol = solidModelDict().lookupOrDefault<scalar>("solutionDTolerance", 1e-11);

    // Pressure tolerance
    const scalar PTol = solidModelDict().lookupOrDefault<scalar>("solutionPTolerance", 1e-6);

    // Check for convergence
    if (residualDNorm < DTol && residualPNorm < PTol)
    {
        Info<< "    Converged" << endl;
        return true;
    }
    else if (iCorr >= nCorr() - 1)
    {
        if (nCorr() > 1)
        {
            Warning
                << "Max iterations reached within the momentum Newton-Raphson "
                "loop" << endl;
        }

        return true;
    }

    // Convergence has not been reached
    return false;
}


//scalar vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::calculateLineSearchSlope
//(
//    const scalar eta,
//    const vectorField& pointDcorr,
//    pointVectorField& pointD,
//    surfaceTensorField& dualGradDf,
//    surfaceSymmTensorField& dualSigmaf,
//    const scalar zeta
//)
//{
//    // Store pointD as we will reset it after changing it
//    pointD.storePrevIter();

//    // Update pointD
//#ifdef OPENFOAM_NOT_EXTEND
//    pointD.primitiveFieldRef() += eta*pointDcorr;
//#else
//    pointD.internalField() += eta*pointDcorr;
//#endif
//    pointD.correctBoundaryConditions();

//    // Calculate gradD at dual faces
//    dualGradDf = vfvc::fGrad
//    (
//        pointD,
//        mesh(),
//        dualMesh(),
//        dualMeshMap().dualFaceToCell(),
//        dualMeshMap().dualCellToPoint(),
//        zeta
//    );

//    // Update F
//    dualFf_ = I + dualGradDf.T();

//    // Update Finv
//    dualFinvf_ = inv(dualFf_);

//    // Update J
//    dualJf_ = det(dualFf_);

//    // Calculate stress at dual faces
//    dualMechanicalPtr_().correct(dualSigmaf);

//    // Update the source vector
//    vectorField source(mesh().nPoints(), vector::zero);
//    pointD.correctBoundaryConditions();
//    updateSource(source, dualMeshMap().dualCellToPoint());

//    // Reset pointD
//    pointD = pointD.prevIter();

//    // Return the slope
//    return gSum(pointDcorr & source);
//}


//scalar vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::calculateLineSearchFactor
//(
//    const scalar rTol, // Slope reduction tolerance
//    const int maxIter, // Maximum number of line search iterations
//    const vectorField& pointDcorr, // Point displacement correction
//    const vectorField& source, // Linear system source
//    const scalar zeta // Discretisation parameter
//)
//{
//    // Calculate initial slope
//    const scalar s0 = gSum(pointDcorr & source);

//    // Set initial line search parameter
//    scalar eta = 1.0;
//    int lineSearchIter = 0;

//    // Perform backtracking to find suitable eta
//    do
//    {
//        lineSearchIter++;

//        // Calculate slope at eta
//        const scalar s1 = calculateLineSearchSlope
//        (
//            eta, pointDcorr, pointD(), dualGradDf_, dualSigmaf_, zeta
//        );

//        // Calculate ratio of s1 to s0
//        const scalar r = s1/s0;

//        if (mag(r) < rTol)
//        {
//            break;
//        }
//        else
//        {
//            // Interpolate/extrapolate to find new eta
//            // Limit it to be less than 10
//            //eta = min(-1/(r - 1), 10);

//            if (r < 0)
//            {
//                // Simple back tracking
//                eta *= 0.5;
//            }
//            else
//            {
//                // Extrapolate
//                eta = min(-1/(r - 1), 10);
//            }
//        }

//        if (lineSearchIter == maxIter)
//        {
//            Warning
//                << "Max line search iterations reached!" << endl;
//        }
//    }
//    while (lineSearchIter < maxIter);

//    // Update pointD and re-calculate source, then calculate s
//    if (mag(eta - 1) > SMALL)
//    {
//        Info<< "        line search parameter = " << eta
//            << ", iter = " << lineSearchIter << endl;
//    }

//    return eta;
//}

Foam::tmp<Foam::Field<Foam::RectangularMatrix<Foam::scalar>>>
vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::geometricStiffnessField
(
        const surfaceVectorField SfUndef,
        const surfaceTensorField gradDRef
) const
{
    // Prepare tmp field
    tmp<Field<Foam::RectangularMatrix<Foam::scalar>>> tresult
    (
        new Field<Foam::RectangularMatrix<Foam::scalar>>(dualMesh().nFaces(), Foam::RectangularMatrix<scalar>(3,9,0))
    );
#ifdef OPENFOAM_NOT_EXTEND
    Field<Foam::RectangularMatrix<Foam::scalar>>& result = tresult.ref();
#else
    Field<Foam::RectangularMatrix<Foam::scalar>>& result = tresult();
#endif

    //For small strain the geometric stiffness is zero
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


Foam::tmp<tensorField>
vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::pBarSensitivityField
(
    const pointTensorField pGradDRef
) const
{
    // Prepare tmp field
    tmp<tensorField> tresult
    (
        new tensorField(mesh().nPoints(), Foam::tensor::zero)
    );
#ifdef OPENFOAM_NOT_EXTEND
    tensorField& result = tresult.ref();
#else
    tensorField& result = tresult();
#endif

    //Calculate the bulk modulus
    const scalar E = 70e9;
    const scalar nu = 0.3;
    const scalar K( (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*(E/(2.0*(1.0 + nu))) );

    //d(pBar)/d(gradD) = -0.5*K*d(J^2)/d(gradD)
    
    // Calculate unperturbed F
    const tensorField FRef(I + pGradDRef.T());

    // Calculate unperturbed J
    const scalarField JRef(det(FRef));

    // Calculate unperturbed J
    const scalarField pBarRef(-0.5*K*(pow(JRef, 2.0) - 1.0)/JRef);

    // Create field to be used for perturbations
    pointTensorField pGradDPerturb("pGradDPerturb", pGradDRef);

    // Small number used for perturbations
    const scalar eps(solidModelDict().lookupOrDefault<scalar>("tangentEps", 1e-10));

    // For each component of gradD, sequentially apply a perturbation and
    // then calculate the resulting sigma
    for (label cmptI = 0; cmptI < tensor::nComponents; cmptI++)
    {
        // Reset gradDPerturb and multiply by 1.0 to avoid it being removed
        // from the object registry
        pGradDPerturb = 1.0*pGradDRef;

        // Perturb this component of pGradD field
        forAll(pGradDPerturb, pointI)
        {
            pGradDPerturb[pointI].component(cmptI) = pGradDRef[pointI].component(cmptI) + eps;
        }

        const tensorField FPerturb(I + pGradDPerturb.T());
        const scalarField JPerturb(det(FPerturb));
        const scalarField pBarPerturb(-0.5*K*(pow(JPerturb, 2.0) - 1.0)/JPerturb);

        // Calculate each component
        const scalarField tangCmpt((pBarPerturb - pBarRef)/eps);
        
//        Info << "pBarPerturb: " << pBarPerturb << endl;
//        Info << "pBarRef: " << pBarRef << endl;

        // Insert components
        forAll(tangCmpt, pointI)
        {
            if (cmptI == tensor::XX)
            {
                result[pointI][tensor::XX] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::XY)
            {
                result[pointI][tensor::XY] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::XZ)
            {
                result[pointI][tensor::XZ] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::YX)
            {
                result[pointI][tensor::YX] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::YY)
            {
                result[pointI][tensor::YY] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::YZ)
            {
                result[pointI][tensor::YZ] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::ZX)
            {
                result[pointI][tensor::ZX] = tangCmpt[pointI];
            }
            else if (cmptI == tensor::ZY)
            {
                result[pointI][tensor::ZY] = tangCmpt[pointI];
            }
            else // if (cmptI == tensor::ZZ)
            {
                result[pointI][tensor::ZZ] = tangCmpt[pointI];
            }
        }
    }
    
    //Info << result[0] << endl;
    //Info << pGradDRef << endl;
    //pointScalarField pointJ = det(I + pGradDRef.T());
    //Info << pointJ << endl;
//    tensorField pBar = -0.5*K*( 
//    	(2*(I + tr(pGradDRef)*I - pGradDRef.T() + cof(pGradDRef))) -
//    	(I + tr(pGradDRef)*I - pGradDRef.T() + cof(pGradDRef))*(1 - 1/sqr(pointJ))
//    	);
    
    //(sqr(pointJ) - 1)*(I + tr(pGradDRef)*I - pGradDRef.T() + cof(pGradDRef)));
    //Info << pBar[0] << endl;
    	
    return tresult;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::vertexCentredNonLinGeomTotalLagPressureDisplacementSolid
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
    fullNewton_(solidModelDict().lookup("fullNewton")),
    geometricStiffness_(solidModelDict().lookup("geometricStiffness")),
    steadyState_(false),
    compactStencil_(false),
    pressureSmoothing_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "pressureSmoothing", 0.5
        )
    ),
//    K_(mechanical().bulkModulus()),
//    Kf_(dualMechanicalPtr_().bulkModulus()),
    twoD_(sparseMatrixExtendedTools::checkTwoD(mesh())),
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
    pointJ_
    (
        IOobject
        (
            "pointJ",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        "calculated"
    ),
    pointGradD_
    (
        IOobject
        (
            "pGrad(D)",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
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
        dimensionedScalar("zero", dimless, 0.0),
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
			IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        "calculated"
    ),
    globalPointIndices_(mesh())
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
#ifdef OPENFOAM_NOT_EXTEND
    scalarField& pointVolI = pointVol_;
    // scalarField& pointVolI = pointVol_.primitiveFieldRef();
#else
    scalarField& pointVolI = pointVol_.internalField();
#endif
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
    }

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointP_.oldTime().storeOldTime();
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

vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::~vertexCentredNonLinGeomTotalLagPressureDisplacementSolid()
{
#ifdef USE_PETSC
    if (Switch(solidModelDict().lookup("usePETSc")))
    {
        PetscFinalize();
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;
    
    ////// Prepare fields at the beginning of each time step //////

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.2));
    if (debug)
    {
        Info<< "zeta: " << zeta << endl;
    }
    
    //Calculate the bulk modulus
    const scalar E = 70e9;
    const scalar nu = 0.3;
    const scalar K( (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*(E/(2.0*(1.0 + nu))) );

    // Initialise matrix where each coefficient is a 4x4 tensor
    sparseMatrixExtended matrixExtended(sum(globalPointIndices_.stencilSize()));

    // Calculate F
    dualFf_ = I + dualGradDf_.T();

    // Calculate Finv
    dualFinvf_ = inv(dualFf_);

    // Calculate J
    dualJf_ = det(dualFf_);

    // Interpolate the pressure to the dual faces
    dualPf_ = vfvc::interpolate
    (
        pointP_,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        debug
    );  
    
    // Calculate gradD at the primary mesh points
    pointGradD_ = vfvc::pGrad
    (
        pointD(),
        mesh()
    );
    
    // Calculate J at the primary points
    pointJ_ = det(I + pointGradD_.T());
    
    // Calculate cell P
    volP_ = vfvc::interpolate
    (
        pointP_,
        mesh()
    );

    // Store material tangent field for dual mesh faces
    Field<RectangularMatrix<scalar>> materialTangent
    (
        dualMechanicalPtr_().materialTangentFaceField()
    );

    // Calculate stress field at dual faces
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Calculate stress for primary cells
    mechanical().correct(sigma());

    // Global point index lists
    const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    const labelList& localToGlobalPointMap =
        globalPointIndices_.localToGlobalPointMap();

    // Coupled pressure and displacement correction
    Field<RectangularMatrix<scalar>> pointDPcorr(pointD().internalField().size(),RectangularMatrix<scalar>(4,1,0));

    // Newton-Raphson loop over momentum equation
    int iCorr = 0;
    scalar initResidualD = 0.0;
    scalar initResidualP = 0.0;
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solverPerf;
#else
    BlockSolverPerformance<vector> solverPerf;
#endif
    do
    {	
    	////// Update fields at the beginning of each outer iteration //////
        
        // Update gradD at dual faces
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
        
        // Update the pressure at the dual faces
        dualPf_ = vfvc::interpolate
        (
            pointP_,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            debug
        );

        // Update gradD at the primary mesh points
        pointGradD_ = vfvc::pGrad
        (
            pointD(),
            mesh()
        );
        
        // Update J at the primary points
    	pointJ_ = det(I + pointGradD_.T());
    	
	    // Update cell P
		volP_ = vfvc::interpolate
		(
		    pointP_,
		    mesh()
		);

        // Calculate stress at dual faces
        dualMechanicalPtr_().correct(dualSigmaf_);

        // Create the source vector for displacement-pressure implementation
        Field<RectangularMatrix<scalar>> sourceExtended(mesh().nPoints(), RectangularMatrix<scalar>(4,1,0));

        pointP_.correctBoundaryConditions();
        pointD().correctBoundaryConditions();
        
        ////// Assemble the source //////

        updateSource
        (
            sourceExtended,
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            pressureSmoothing_,
            zeta,
            debug
        );

        ////// Assemble the matrix //////

        // Assemble the matrix once per outer iteration
        matrixExtended.clear();

        // Update material tangent
        materialTangent = dualMechanicalPtr_().materialTangentFaceField();

        //Obtain undeformed surface vector field
        surfaceVectorField SfUndef = dualMesh().Sf();

        // Calculate geometric stiffness field for dual mesh faces
        Foam::Field<Foam::RectangularMatrix<Foam::scalar>> geometricStiffness
        (
            geometricStiffnessField
            (
                SfUndef,
                dualGradDf_
            )
        );

        // Calculate pBarSensitivity
        tensorField pBarSensitivity = -0.5*K*(
        	I + tr(pointGradD_)*I - pointGradD_.T() + cof(pointGradD_) +
        	pow(pointJ_,-2.0)*(I + tr(pointGradD_)*I - pointGradD_.T() + cof(pointGradD_))); 
        	
//        tensorField pBarSensitivity
//        (
//            pBarSensitivityField
//            (
//                pointGradD_
//            )
//        );
//        	
//        for (int i = 0; i < 60; i++)
//        {
//        	Info << pBarSensitivity[i] << endl;
//        }
        //Info << pBarSensitivity << endl;  
        
        // Add div(sigma) pressure and displacement coefficients
        vfvm::divSigma
        (
            matrixExtended,
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

        // Add laplacian coefficient to the pressure equation
        vfvm::laplacian
        (
            matrixExtended,
            compactStencil_,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            dualGradDf_,
            pressureSmoothing_,
            debug
        );

        // Add coefficients of pressure equation
        vfvm::Sp
        (
            matrixExtended,
            mesh(),
            dualMeshMap().dualCellToPoint(),
            pointVol_,
            pBarSensitivity,
            debug
        );           
        
//		    Info << endl << "Before enforcing DOFs: " << endl << endl;
//		    matrixExtended.print();
//		    Info << endl << "Print out the source: " << endl << endl;

//		    for (int i = 0; i < sourceExtended.size(); i++)
//		    {
//                      Info << "(" << i << ", 0) : " << sourceExtended[i] << endl;

//		    }
//		    Info << endl;

        // Enforce fixed DOF on the linear system for
        // the displacement
        sparseMatrixExtendedTools::enforceFixedDof
        (
            matrixExtended,
            sourceExtended,
            twoD_,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofValues_,
            fixedDofScale_
        );

//		    Info << endl << "After enforcing DOFs " << endl << endl;
//		    matrixExtended.print();
//		    Info << endl << "Print out the source: " << endl << endl;

//		    for (int i = 0; i < sourceExtended.size(); i++)
//		    {
//                      Info << "(" << i << ", 0) : " << sourceExtended[i] << endl;
//		    }

        ////// Solve the linear system //////
        
        if (debug)
        {
            Info<< "bool vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::evolve(): "
                << " solving linear system: start" << endl;
        }
        
        Info<< "    Solving" << endl;

        if (Switch(solidModelDict().lookup("usePETSc")))
        {
#ifdef USE_PETSC
            fileName optionsFile(solidModelDict().lookup("optionsFile"));

            // Solve for displacement and pressure correction

            solverPerf = sparseMatrixExtendedTools::solveLinearSystemPETSc
            (
                matrixExtended,
                sourceExtended,
                pointDPcorr,
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
            FatalErrorIn("vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::evolve()")
                << "PETSc not available. Please set the PETSC_DIR environment "
                << "variable and re-compile solids4foam" << abort(FatalError);
#endif
        }
        else
        {
            // Use Eigen SparseLU direct solver
            notImplemented("Not implemented with Eigen SparseLU direct solver")
        }

        if (debug)
        {
            Info<< "bool vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::evolve(): "
                << " solving linear system: end" << endl;
        }

        // Update point displacement field
        if (Switch(solidModelDict().lookup("lineSearch")))
        {
            notImplemented("Line search not implemented.")
        }
#ifdef OPENFOAM_NOT_EXTEND
        else if (mesh().relaxField(pointD().name()))
#else
        else if (mesh().solutionDict().relaxField(pointD().name()))
#endif
        {
            notImplemented("pointD or pointP relaxation not implemented.")
        }
        else
        {
            forAll(pointDPcorr, pointI)
            {
#ifdef OPENFOAM_NOT_EXTEND
                pointD().primitiveFieldRef()[pointI].component(0) += pointDPcorr[pointI](0,0);
                pointD().primitiveFieldRef()[pointI].component(1) += pointDPcorr[pointI](1,0);
                pointD().primitiveFieldRef()[pointI].component(2) += pointDPcorr[pointI](2,0);
                pointP_.primitiveFieldRef()[pointI] += pointDPcorr[pointI](3,0);
#else
                pointD().internalField()[pointI].component(0) += pointDPcorr[pointI](0,0);
                pointD().internalField()[pointI].component(1) += pointDPcorr[pointI](1,0);
                pointD().internalField()[pointI].component(2) += pointDPcorr[pointI](2,0);
                pointP_.internalField()[pointI] += pointDPcorr[pointI](3,0);
#endif
            }
        }
		
        pointD().correctBoundaryConditions();
        pointP_.correctBoundaryConditions();

        // Update point accelerations
        // Note: for NewmarkBeta, this needs to come before the pointU update
#ifdef OPENFOAM_NOT_EXTEND
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
#else
        pointA_.internalField() =
            vfvc::ddt
            (
                mesh().schemesDict().ddtScheme("ddt(pointU)"),
                mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
                pointU_
            );

        // Update point velocities
        pointU_.internalField() =
            vfvc::ddt
            (
                mesh().schemesDict().ddtScheme("ddt(pointD)"),
                mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
                pointD()
            );
#endif
    }
    while
    (
        !converged
        (
            iCorr,
            initResidualD,
            initResidualP,
#ifdef OPENFOAM_NOT_EXTEND
            cmptMax(solverPerf.nIterations()),
#else
            solverPerf.nIterations(),
#endif
            pointD(),
            pointP_,
            pointDPcorr
        ) && ++iCorr
    );
    
    ////// Update fields at the end of each time step //////

    // Update gradD at dual faces
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
    
    // Update the pressure at the dual faces
    dualPf_ = vfvc::interpolate
    (
        pointP_,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        debug
    );    

    // Update gradD at the primary mesh points
    pointGradD_ = vfvc::pGrad
    (
        pointD(),
        mesh()
    );
    
    // Update J at the primary points
	pointJ_ = det(I + pointGradD_.T());
    
    // Update cell P
    volP_ = vfvc::interpolate
    (
        pointP_,
        mesh()
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

    return true;
}

void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::setTraction
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
#ifdef OPENFOAM_NOT_EXTEND
    pointPatchVectorField& ptPatch = pointD().boundaryFieldRef()[patchID];
#else
    pointPatchVectorField& ptPatch = pointD().boundaryField()[patchID];
#endif

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
            "void Foam::vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::setTraction\n"
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

void vertexCentredNonLinGeomTotalLagPressureDisplacementSolid::writeFields(const Time& runTime)
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

#ifdef FOAMEXTEND
    pEpsilon.internalField() = symm(pGradD.internalField());
#else
    pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
#endif
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

#ifdef FOAMEXTEND
    pEpsilonEq.internalField() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#else
    pEpsilonEq.primitiveFieldRef() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#endif
    pEpsilonEq.write();

    Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << endl;

    solidModel::writeFields(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
