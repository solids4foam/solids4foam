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

#include "vertexCentredSegregatedLinGeomSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrix.H"
#include "symmTensor4thOrder.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "sparseMatrixTools.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredSegregatedLinGeomSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, vertexCentredSegregatedLinGeomSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void vertexCentredSegregatedLinGeomSolid::updateSource
(
    vectorField& source,
    const labelList& dualCellToPoint
)
{
    if (debug)
    {
        Info<< "void vertexCentredSegregatedLinGeomSolid::updateSource(...)"
            << endl;
    }

    // Reset to zero
    source = vector::zero;

    // The source vector is -F, where:
    // F = div(sigma) + rho*g - rho*d2dt2(D)

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Point density field
    const scalarField& pointRhoI = pointRho_.internalField();

    // Calculate the tractions on the dual faces
    surfaceVectorField dualTraction
    (
        (dualMesh().Sf()/dualMesh().magSf()) & dualSigmaf_
    );

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
    // const vectorField dualDivSigma = fvc::div(dualMesh().Sf() & dualSigmaf_);
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh().magSf());

    // Map dual cell field to primary mesh point field
    vectorField pointDivSigma(mesh().nPoints(), vector::zero);
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigma[pointID] = dualDivSigma[dualCellI];
    }

    // Add surface forces to source
    source -= pointDivSigma*pointVolI;

    // Add gravity body forces
    source -= pointRhoI*g().value()*pointVolI;

//     // Add transient term
//     source += vfvc::d2dt2
//     (
// #ifdef OPENFOAM_NOT_EXTEND
//         mesh().d2dt2Scheme("d2dt2(pointD)"),
// #else
//         mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
// #endif
//         pointD(),
//         pointU_,
//         pointA_,
//         pointRho_,
//         pointVol_,
//         int(bool(debug))
//     );
}


void vertexCentredSegregatedLinGeomSolid::setFixedDofs
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
                            "void vertexCentredSegregatedLinGeomSolid::setFixedDofs(...)"
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
                            "void vertexCentredSegregatedLinGeomSolid::setFixedDofs(...)"
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


void vertexCentredSegregatedLinGeomSolid::enforceTractionBoundaries
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
                    "void vertexCentredSegregatedLinGeomSolid::"
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

bool vertexCentredSegregatedLinGeomSolid::
vertexCentredSegregatedLinGeomSolid::converged
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
#ifdef OPENFOAM_NOT_EXTEND
    const scalar maxMagD = gMax(mag(pointD.primitiveField()));
#else
    const scalar maxMagD = gMax(mag(pointD.internalField()));
#endif

    // Check for convergence
    bool converged = false;
    if (residualNorm < solutionTol())
    {
        Info<< "    Converged" << endl;
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

        if (iCorr >= nCorr() - 1)
        {
            Warning
                << "Max iterations reached within the momentum loop"
                << endl;
        }
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredSegregatedLinGeomSolid::vertexCentredSegregatedLinGeomSolid
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
    impKf_(dualMechanicalPtr_().impKf()),
    steadyState_(false),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    twoDCorrector_(mesh()),
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
    globalPointIndices_(mesh())
#ifdef OPENFOAM_COM
    ,
    pointVolInterp_(pMesh(), mesh())
#endif
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
    scalarField& pointVolI = pointVol_.primitiveFieldRef();
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

vertexCentredSegregatedLinGeomSolid::~vertexCentredSegregatedLinGeomSolid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool vertexCentredSegregatedLinGeomSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

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
        impKf_.primitiveField(),
        debug
    );

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.0));
    Info<< "zeta: " << zeta << endl;

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
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solverPerf;
#else
    BlockSolverPerformance<vector> solverPerf;
#endif
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

        // Calculate stress at dual faces
        dualMechanicalPtr_().correct(dualSigmaf_);

        // Update the source vector
        source = vector::zero;
        pointD().correctBoundaryConditions();
        updateSource(source, dualMeshMap().dualCellToPoint());

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

#ifdef OPENFOAM_NOT_EXTEND
            pointD().primitiveFieldRef().replace
            (
                dirI, pointD().primitiveField().component(dirI) + pointDcorr
            );
#else
            pointD().internalField().replace
            (
                dirI, pointD().internalField().component(dirI) + pointDcorr
            );
#endif

            pointD().correctBoundaryConditions();
        }

//         // Update point accelerations
//         // Note: for NewmarkBeta, this needs to come before the pointU update
// #ifdef OPENFOAM_NOT_EXTEND
//         pointA_.primitiveFieldRef() =
//             vfvc::ddt
//             (
//                 mesh().ddtScheme("ddt(pointU)"),
//                 mesh().d2dt2Scheme("d2dt2(pointD)"),
//                 pointU_
//             );

//         // Update point velocities
//         pointU_.primitiveFieldRef() =
//             vfvc::ddt
//             (
//                 mesh().ddtScheme("ddt(pointD)"),
//                 mesh().d2dt2Scheme("d2dt2(pointD)"),
//                 pointD()
//             );
// #else
//         pointA_.internalField() =
//             vfvc::ddt
//             (
//                 mesh().schemesDict().ddtScheme("ddt(pointU)"),
//                 mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
//                 pointU_
//             );

//         // Update point velocities
//         pointU_.internalField() =
//             vfvc::ddt
//             (
//                 mesh().schemesDict().ddtScheme("ddt(pointD)"),
//                 mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
//                 pointD()
//             );
// #endif

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
            solverPerf.finalResidual()[vector::X],
#ifdef OPENFOAM_NOT_EXTEND
            cmptMax(solverPerf.nIterations()),
#else
            solverPerf.nIterations(),
#endif
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

#ifdef OPENFOAM_COM
    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());
#endif

    return true;
}

void vertexCentredSegregatedLinGeomSolid::setTraction
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
            "void Foam::vertexCentredSegregatedLinGeomSolid::setTraction\n"
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


void vertexCentredSegregatedLinGeomSolid::writeFields(const Time& runTime)
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
