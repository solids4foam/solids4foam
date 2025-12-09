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

#include "vertexCentredLinGeomSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#include "linearElasticMisesPlastic.H"
#include "compatibilityFunctions.H"
#include <vector>

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
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 1.0));

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
            boundaryFieldRef(dualTraction)[patchI] = vector::zero;
        }
    }

    // Calculate divergence of stress (force per unit volume) for the dual
    // cells
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

#ifdef OPENFOAM_NOT_EXTEND
    // Sum absolute forces in parallel
    pointConstraints::syncUntransformedData
    (
        mesh(), pointDivSigma, plusEqOp<vector>()
    );
#else
    if (Pstream::parRun())
    {
        notImplemented
        (
            "Running " + type() + " in parallel us currently only possible in "
            "OpenFOAM.com versions"
        );
    }
#endif

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
    symmTensorField& fixedDofDirections,
    vectorField& fixedDofDirectionsVec
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
                    fixedDofDirectionsVec[pointID] = vector::one;
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = disp;
                    fixedDofDirections[pointID] = symmTensor(I);
                    fixedDofDirectionsVec[pointID] = vector::one;
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

                    // If the point is not fully fixed then make sure the
                    // normal direction is fixed
                    if (mag(fixedDofDirections[pointID] - symmTensor(I)) > 0)
                    {
                        // If the directions are orthogonal we can add them
                        const symmTensor curDir = sqr(pointNormals[pI]);
                        if (mag(fixedDofDirections[pointID] & curDir) > 0)
                        {
                            FatalError
                                << "Point " << pointID << " is fixed in two "
                                << "directions: this is only implemented for "
                                << "Cartesian axis directions"
                                << abort(FatalError);
                        }

                        fixedDofDirections[pointID] += curDir;
                        fixedDofDirectionsVec[pointID] += pointNormals[pI];
                    }
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = normalDisp[pI]*pointNormals[pI];
                    fixedDofDirections[pointID] = sqr(pointNormals[pI]);
                    fixedDofDirectionsVec[pointID] = pointNormals[pI];
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

            // Map from primary mesh point field to second mesh face field
            // using the pointToDualFaces map
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
            boundaryFieldRef(dualTraction)[patchI] = dualFaceTraction;
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
            boundaryFieldRef(dualTraction)[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
        }
    }
}


#ifdef USE_PETSC

void vertexCentredLinGeomSolid::makeFixedDofRowsIS() const
{
    if (fixedDofRowsISPtr_ != nullptr)
    {
        FatalError
            << "Pointer already set" << exit(FatalError);
    }

    // Mark all fixed degrees of freedom
    List<boolList> mask
    (
        mesh().nPoints(), boolList(blockSize_, false)
    );

    const boolList& ownedByThisProc = globalPoints().ownedByThisProc();
    forAll(mask, pointI)
    {
        mask[pointI] = boolList(blockSize_, false);

        if (ownedByThisProc[pointI])
        {
            if (fixedDofs_[pointI])
            {
                for (label cmptI = 0; cmptI < blockSize_; ++cmptI)
                {
                    if (mag(fixedDofDirectionsVec_[pointI][cmptI]) > SMALL)
                    {
                        mask[pointI][cmptI] = true;
                    }
                }
            }
        }
    }

    std::vector<PetscInt> rows;
    const labelList& localToGlobalPointMap =
        globalPoints().localToGlobalPointMap();
    rows.reserve(localToGlobalPointMap.size()*blockSize_);

    forAll(localToGlobalPointMap, i)
    {
        if (ownedByThisProc[i])
        {
            // global block row (node id)
            const label gBlock = localToGlobalPointMap[i];

            for (PetscInt c = 0; c < blockSize_; ++c)
            {
                if (mask[i][c])
                {
                    // scalar global row
                    rows.push_back((PetscInt)gBlock*blockSize_ + c);
                }
            }
        }
    }

    ISCreateGeneral
    (
        PETSC_COMM_WORLD,
        (PetscInt)rows.size(),
        rows.data(),
        PETSC_COPY_VALUES,
        &fixedDofRowsISPtr_
    );
}

#endif // USE_PETSC


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

    const word predictorMethod
    (
        solidModelDict().lookupOrDefault<word>("predictorMethod", "linear")
    );

    if (predictorMethod == "linear")
    {
        // Assuming constant velocity
        pointD() = pointD().oldTime() + pointU_*runTime().deltaT();
    }
    else if (predictorMethod == "quadratic")
    {
        // Assuming constant acceleration
        pointD() =
            pointD().oldTime()
          + pointU_*runTime().deltaT()
          + 0.5*pointA_*pow(runTime().deltaT(), 2);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown predictorMethod = " << predictorMethod << ". Available "
            << "options are 'linear' and 'quadratic'" << exit(FatalError);
    }
}


bool vertexCentredLinGeomSolid::evolveSnes()
{
#ifdef USE_PETSC

    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update pointD boundary conditions
    pointD().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Map the pointD field to the SNES solution vector
        // Note: for point fields, the SNES solution vector may be smaller than
        // pointD.size() because points on processor boundaries may be owned by
        // other processors
        vectorField& pointDI = pointD();
        foamPetscSnesHelper::InsertFieldComponents<vector>
        (
            pointDI,
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD()
          ? makeList<label>({0,1})
          : makeList<label>({0,1,2})
        );
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the D field
    vectorField& pointDI = pointD();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        pointDI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    pointD().correctBoundaryConditions();

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 1.0));

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

    // Update point accelerations and velocities
    // Note: the acceleration needs to be updated before the
    // velocity for NewmarkBeta
    vectorField& pointAI = pointA_;
    vectorField& pointUI = pointU_;
    pointAI = vfvc::ddt(mesh(), pointU_);
    pointUI = vfvc::ddt(mesh(), pointD());

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
    // This is useful for visualisation but it is also needed when using
    // preCICE
    pointVolInterp_.interpolate(pointD(), D());
#endif

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

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

    if (solidModel::twoD())
    {
        solidModel::twoDCorrector().correctPoints(pointD());

        // Remove displacement in the empty directions
        forAll(mesh().geometricD(), dirI)
        {
            if (mesh().geometricD()[dirI] < 0)
            {
#ifdef OPENFOAM_NOT_EXTEND
                pointD().primitiveFieldRef().replace(dirI, 0.0);
#else
                pointD().internalField().replace(dirI, 0.0);
#endif
            }
        }
    }

    // Update the divergence of stress based on the latest pointD field
    updatePointDivSigma(pointD(), dualGradDf_, dualSigmaf_, pointDivSigma_);

    // Compute acceleration
#ifdef OPENFOAM_NOT_EXTEND
    pointA_ = pointDivSigma_/pointRho_ - dampingCoeff()*pointU_ + g();
#else
    pointA_ = pointDivSigma_/pointRho_ - dampingCoeff()*pointU_;

    if (mag(g().value()) > SMALL)
    {
        // foam-extend does not implement the addition of a uniform dimensioned
        // field to a geometric point field so we will do it manually
        vectorField& pointAI = pointA_;
        const vector gVec(g().value());
        forAll(pointAI, pointI)
        {
            pointAI[pointI] += gVec;
        }
        pointA_.correctBoundaryConditions();
    }
#endif

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
#ifdef USE_PETSC
    foamPetscSnesHelper
    (
        "pointD",
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solutionLocation::POINTS,
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
    ),
#endif
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
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    ),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    fixedDofs_(mesh().nPoints(), false),
    fixedDofValues_(fixedDofs_.size(), vector::zero),
    fixedDofDirections_(fixedDofs_.size(), symmTensor::zero),
    fixedDofDirectionsVec_(fixedDofs_.size(), vector::zero),
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
#ifdef USE_PETSC
    fixedDofRowsISPtr_(nullptr),
#endif
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
    )
#ifdef OPENFOAM_COM
    ,
    pointVolInterp_(pMesh(), mesh())
#endif
{
    if (solvePressure())
    {
        notImplemented("Not implemented when solvePressure is active");
    }

    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // pointD field must be defined
    pointDisRequired();

    // Set fixed degree of freedom list
    setFixedDofs
    (
        pointD(),
        fixedDofs_,
        fixedDofValues_,
        fixedDofDirections_,
        fixedDofDirectionsVec_
    );

    // Set point density field
    mechanical().volToPoint().interpolate(rho(), pointRho_);

    // Set the pointVol field
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

#ifdef OPENFOAM_NOT_EXTEND
    // Sum the shared point volumes to create the point global volumes
    pointConstraints::syncUntransformedData
    (
        mesh(), pointGlobalVol_, plusEqOp<scalar>()
    );
#else
    if (Pstream::parRun())
    {
        notImplemented
        (
            "Running " + type() + " in parallel us currently only possible in "
            "OpenFOAM.com versions"
        );
    }
#endif

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointU_.oldTime().storeOldTime();
    pointA_.storeOldTime();

    // Write fixed degree of freedom equation scale
    Info<< type() << "fixedDofScale: " << fixedDofScale_ << endl;

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
    if (fixedDofRowsISPtr_ != nullptr)
    {
        ISDestroy(&fixedDofRowsISPtr_);
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
        // For safety, we should use a time-step smaller than this e.g.
        // Abaqus uses stableTimeStep/sqrt(2): we will default to this value
        const scalar requiredDeltaT =
            1.0/
            gMax
            (
#ifdef OPENFOAM_NOT_EXTEND
                DimensionedField<scalar, Foam::surfaceMesh>
#else
                Field<scalar>
#endif
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


bool vertexCentredLinGeomSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    {
#ifdef OPENFOAM_COM
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Use "
            << solutionAlgorithmNames_.names()[solutionAlgorithm::PETSC_SNES]
            << " instead of "
            << solutionAlgorithmNames_.names()[solutionAlgorithm::IMPLICIT_COUPLED]
            << exit(FatalError);
#else
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Use "
            << solutionAlgorithmNames_[solutionAlgorithm::PETSC_SNES]
            << " instead of "
            << solutionAlgorithmNames_[solutionAlgorithm::IMPLICIT_COUPLED]
            << exit(FatalError);
#endif

        // Keep compiler happy
        return true;
    }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
#ifdef OPENFOAM_COM
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << solutionAlgorithmNames_.names()[IMPLICIT_SEGREGATED]
            << " is not implemented. The behaviour can be mimicked with "
            << solutionAlgorithmNames_.names()[solutionAlgorithm::PETSC_SNES]
            << exit(FatalError);
#else
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << solutionAlgorithmNames_[IMPLICIT_SEGREGATED]
            << " is not implemented. The behaviour can be mimicked with "
            << solutionAlgorithmNames_[solutionAlgorithm::PETSC_SNES]
            << exit(FatalError);
#endif

        // Keep compiler happy
        return true;
    }
    else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
    else
    {
#ifdef OPENFOAM_COM
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            << solutionAlgorithmNames_.names() << exit(FatalError);
#else
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            << solutionAlgorithmNames_.toc() << exit(FatalError);
#endif
    }

    // Keep compiler happy
    return true;
}


#ifdef USE_PETSC

label vertexCentredLinGeomSolid::initialiseJacobian(Mat& jac)
{
    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label vertexCentredLinGeomSolid::initialiseSolution(Vec& x)
{
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label vertexCentredLinGeomSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();

    // Extract pointD from x
    vectorField& pointDI = pointD();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        pointDI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the displacement boundary conditions
    pointD().correctBoundaryConditions();

    // Update the divergence of stress based on the latest pointD field
    // This also updates dualGradDf and dualSigmaf
    updatePointDivSigma(pointD(), dualGradDf_, dualSigmaf_, pointDivSigma_);

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Point density field
    const scalarField& pointRhoI = pointRho_.internalField();

    // The residual vector F calculated as:
    // F = div(sigma) + rho*g - rho*d2dt2(D)

    // Note: we integrate the residual here over the local point volumes as
    // opposed to the global point volumes but it does not matter: it does not
    // affect the Jacobian and only the residual from the processor which owns
    // the point is used.
    vectorField residual
    (
        pointDivSigma_*pointVolI
      + pointRhoI*g().value()*pointVolI
      - vfvc::d2dt2
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh.d2dt2Scheme("d2dt2(pointD)"),
#else
            mesh.schemesDict().d2dt2Scheme("d2dt2(pointD)"),
#endif
            pointD(),
            pointU_,
            pointA_,
            pointRho_,
            pointVol_,
            int(bool(debug))
        )
    );

    // Enforce fixed DOF by setting the residual to zero
    {
        const boolList& ownedByThisProc = globalPoints().ownedByThisProc();
        forAll(residual, pointI)
        {
            if (ownedByThisProc[pointI])
            {
                if (fixedDofs_[pointI])
                {
                    for (label cmptI = 0; cmptI < blockSize_; ++cmptI)
                    {
                        if (mag(fixedDofDirectionsVec_[pointI][cmptI]) > SMALL)
                        {
                            residual[pointI][cmptI] = 0.0;
                        }
                    }
                }
            }
        }
    }

    // Insert the residual into the PETSc residual
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    return 0;
}


label vertexCentredLinGeomSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();

    // Extract pointD from x
    vectorField& pointDI = pointD();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        pointDI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the displacement boundary conditions
    pointD().correctBoundaryConditions();

    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 1.0));
    const scalar zetaImplicit
    (
        solidModelDict().lookupOrDefault<scalar>("zetaImplicit", zeta)
    );

    // Calculate gradD at dual faces
    dualGradDf_ = vfvc::fGrad
    (
        pointD(),
        mesh,
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta,
        debug
    );

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf_);

    if (solidModelDict().lookupOrDefault<Switch>("approximateJacobian", false))
    {
        // Add laplacian term as a compact approximate linearisation of
        // div(sigma)
        vfvm::laplacian
        (
            jac,
            Switch(solidModelDict().lookup("compactImplicitStencil")),
            zetaImplicit,
            dualMesh(),
            blockSize_,     // nScalarEqns
            globalPoints().localToGlobalPointMap(),
            dualImpKf(),
            false           // flip sign
        );
    }
    else
    {
        // Calculate the material tangent
        List<mat66> materialTangent(mesh.nFaces());
        dualMechanicalPtr_().materialTangentFaceField(materialTangent);

        // Add linearisation of div(sigma) to jac
        vfvm::divSigma
        (
            jac,
            pointD(),
            mesh,
            dualMesh(),
            blockSize_,     // nScalarEqns
            globalPoints().localToGlobalPointMap(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            materialTangent,
            zetaImplicit,
            false           // flip sign
        );
    }

    // Lookup the d2dt2 scheme
#ifdef OPENFOAM_NOT_EXTEND
    ITstream& d2dt2Scheme = mesh.d2dt2Scheme("d2dt2(pointD)");
#else
    ITstream& d2dt2Scheme = mesh.schemesDict().d2dt2Scheme("d2dt2(pointD)");
#endif

    // Add d2dt2 coefficients to jac
    vfvm::d2dt2
    (
        jac,
        pointD(),
        pointRho_,
        pointVol_,
        d2dt2Scheme,
        blockSize_,     // nScalarEqns
        globalPoints().localToGlobalPointMap(),
        true           // flip sign
    );

    // Complete matrix assembly: this is required before we call
    // MatZeroRowsColumnsIS
    CHKERRQ(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
    CHKERRQ(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

    // Enforce fixed DOF
    // Zero all rows and columns of fixed DOFs and set -fixedDofScale_ on
    // the diagonal
    MatZeroRowsColumnsIS(jac, fixedDofRowsIS(), -fixedDofScale_, NULL, NULL);


    if (solvePressure())
    {
        notImplemented("solvePressure not implemented yet for formJacobian");
    }

    return 0;
}

#endif // USE_PETSC


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
    pointPatchVectorField& ptPatch = boundaryFieldRef(pointD())[patchID];

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

#ifdef OPENFOAM_NOT_EXTEND
    pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
#else
    pEpsilon.internalField() = symm(pGradD.internalField());
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

#ifdef OPENFOAM_NOT_EXTEND
    pEpsilonEq.primitiveFieldRef() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#else
    pEpsilonEq.internalField() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#endif
    pEpsilonEq.write();

    Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << endl;

    // To-do: add a general interface for calculating point stresses

    solidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
