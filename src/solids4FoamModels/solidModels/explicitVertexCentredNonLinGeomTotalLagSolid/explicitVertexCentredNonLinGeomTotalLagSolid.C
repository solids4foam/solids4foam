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

#include "explicitVertexCentredNonLinGeomTotalLagSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrix.H"
#include "symmTensor4thOrder.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#include "pointConstraints.H"
#include "minEdgeLength.H"
#include "fvc.H"
#include "fvcAverage.H"
#include "fvcLaplacian.H"
#include "sparseMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitVertexCentredNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable(solidModel, explicitVertexCentredNonLinGeomTotalLagSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void explicitVertexCentredNonLinGeomTotalLagSolid::updatePointDivSigma
(
    const pointVectorField& pointD,
    surfaceTensorField& dualGradDf,
    surfaceTensorField& dualFf,
    surfaceTensorField& dualFinvf,
    surfaceScalarField& dualJf,
    surfaceScalarField& dualPf,
    surfaceSymmTensorField& dualSigmaf,
    pointVectorField& pointDivSigma
)
{
    if (debug)
    {
        Info<< "void explicitVertexCentredNonLinGeomTotalLagSolid::"
            << "updatePointDivSigma(...): start" << endl;
    }

    // Calculate gradD at dual faces
    dualGradDf = vfvc::fGrad
    (
        pointD,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta_,
        debug
    );

    // Calculate F at dual faces
    if (useDeformationGradientInvolution_)
    {
        dualFf += dualGradDf.T() - dualGradDf.oldTime().T();
    }
    else
    {
        dualFf = I + dualGradDf.T();
    }

    // Calculate Finv
    dualFinvf = inv(dualFf);

    // Calculate the Jacobian of F
    dualJf = det(dualFf);

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Calculate pressure using the artificial compressibility approach
    // Note: gradU ~= (1/deltaT)*gradDeltaD
    //             ~= (1/deltaT)*(gradD - gradD.old)
    // dualPf is used for post-processing if the useArtificialCompressibility
    // switch is off
    dualPf -=
        dualKf_*dualJf*(dualFinvf && (dualGradDf - dualGradDf.oldTime()));

    if (useArtificialCompressibility_)
    {
        // Replace pressure component in stress
        dualSigmaf = dev(dualSigmaf) - dualPf_*I;
    }

    // Calculate the deformed unit normals
    const surfaceVectorField deformedSf(dualJf*dualFinvf.T() & dualMesh().Sf());
    const surfaceScalarField deformedMagSf(mag(deformedSf));
    const surfaceVectorField deformedN(deformedSf/deformedMagSf);

    // Calculate the tractions at the dual faces
    surfaceVectorField dualTraction(deformedN & dualSigmaf);

    // Enforce traction boundaries
    enforceTractionBoundaries
    (
        pointD,
        dualTraction,
        deformedN,
        mesh(),
        dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
#else
            dualTraction.boundaryField()[patchI] = vector::zero;
#endif
        }
    }

    // Calculate divergence of stress (force per unit volume) for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*deformedMagSf);

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
    const scalarField& pointVolI = pointVol_;
    forAll(pointDivSigmaI, pointI)
    {
        pointDivSigmaI[pointI] /= pointVolI[pointI];
    }

    if (debug)
    {
        Info<< "void explicitVertexCentredNonLinGeomTotalLagSolid::"
            << " updatePointDivSigma(...): end" << endl;
    }
}


void explicitVertexCentredNonLinGeomTotalLagSolid::enforceTractionBoundaries
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

            // Map from primary mesh point field to dual mesh face field using
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
                            const vector n =
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
                    "void explicitVertexCentredNonLinGeomTotalLagSolid::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;

            // Overwrite the dual patch face traction
#ifdef OPENFOAMESIORFOUNDATION
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
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#else
            dualTraction.boundaryField()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#endif
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitVertexCentredNonLinGeomTotalLagSolid::explicitVertexCentredNonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    useArtificialCompressibility_
    (
        solidModelDict().lookup("useArtificialCompressibility")
    ),
    useDeformationGradientInvolution_
    (
        solidModelDict().lookup("useDeformationGradientInvolution")
    ),
    zeta_(solidModelDict().lookupOrDefault<scalar>("zeta", 1.0)),
    dampingCoeff_(solidModelDict().lookup("dampingCoeff")),
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
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    dualKf_("dualKf", dualMechanicalPtr_->bulkModulus()),
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
            "pointRho",
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
    dualFf_
    (
        IOobject
        (
            "dualFf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("I", dimless, tensor(I)),
        "calculated"
    ),
    dualFinvf_
    (
        IOobject
        (
            "dualFinvf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("I", dimless, tensor(I)),
        "calculated"
    ),
    dualJf_
    (
        IOobject
        (
            "dualJf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedScalar("1", dimless, 1.0),
        "calculated"
    ),
    dualPf_
    (
        IOobject
        (
            "dualPf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
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
    nProcPerPoint_
    (
        IOobject
        (
            "nProcPerPoint",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("1", dimless, 1.0)
    ),
    globalPointIndices_(mesh())
#ifdef OPENFOAMESI
    ,
    pointVolInterp_(pMesh(), mesh())
#endif
{
    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // pointD field must be defined
    pointDisRequired();

    // Set point density field
    mechanical().volToPoint().interpolate(rho(), pointRho_);


    // Update nProcPerPoint field
    pointConstraints::syncUntransformedData
    (
        mesh(), nProcPerPoint_, plusEqOp<scalar>()
    );


    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
    scalarField& pointVolI = pointVol_;
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
    }

    // Sum the shared point volumes for consistency with serial runs
    pointConstraints::syncUntransformedData
    (
        mesh(), pointVol_, plusEqOp<scalar>()
    );

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointU_.oldTime().storeOldTime();
    pointA_.storeOldTime();
    dualGradDf_.oldTime();

    // Write flags
    Info<< "useArtificialCompressibility: " << useArtificialCompressibility_
        << nl << "zeta: " << zeta_ << endl;

    // Disable the writing of the unused fields
    D().writeOpt() = IOobject::NO_WRITE;
    D().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    DD().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    U().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

explicitVertexCentredNonLinGeomTotalLagSolid::
~explicitVertexCentredNonLinGeomTotalLagSolid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitVertexCentredNonLinGeomTotalLagSolid::setDeltaT(Time& runTime)
{
    if (time().timeIndex() == 0)
    {
        // Max wave speed in the domain
        scalar waveSpeed = 0.0;

        // if (useArtificialCompressibility_)
        // {
        //     // Calculate effective shear modulus
        //     // Assume impK == (4/3)*mu + bulkModulus
        //     const volScalarField mu
        //     (
        //         0.75*(mechanical().impK() - mechanical().bulkModulus())
        //     );

        //     // Use shear wave speed
        //     waveSpeed = max(Foam::sqrt(mu/mechanical().rho())).value();
        // }
        // else
        {
            // Use pressure wave speed
            waveSpeed = max
            (
                Foam::sqrt(mechanical().impK()/mechanical().rho())
            ).value();
        }

        // deltaT = cellWidth/waveVelocity
        // For safety, we should use a time-step smaller than this e.g. Abaqus
        // uses stableTimeStep/sqrt(2): we will default to this value
        const scalar requiredDeltaT = minEdgeLength(mesh())/waveSpeed;
        //const scalar requiredDeltaT = minEdgeLength(mesh(), pointD())/waveSpeed;

        // Lookup the desired Courant number
        const scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

        const scalar newDeltaT = maxCo*requiredDeltaT;

        physicsModel::printInfo() = false;

        // Info<< nl << "Setting deltaT = " << newDeltaT
        //     << ", maxCo = " << maxCo << endl;

        runTime.setDeltaT(newDeltaT);
    }

    // Update print info
    // physicsModel::printInfo() = bool
    // (
    //     runTime.timeIndex() % infoFrequency() == 0
    //  || mag(runTime.value() - runTime.endTime().value()) < SMALL
    // );
}


bool explicitVertexCentredNonLinGeomTotalLagSolid::evolve()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the solid momentum equation for pointD" << nl
            << "Simulation Time, Time Step, Clock Time, Max Stress" << endl;
    }

    physicsModel::printInfo() = bool
    (
        time().timeIndex() % infoFrequency() == 0
     || mag(time().value() - time().endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< time().value() << " " << time().deltaTValue()
            << " " << time().elapsedClockTime()
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

    // Remove displacement in the empty directions
    forAll(mesh().geometricD(), dirI)
    {
        if (mesh().geometricD()[dirI] < 0)
        {
            pointD().primitiveFieldRef().replace(dirI, 0.0);
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
        dualPf_,
        dualSigmaf_,
        pointDivSigma_
    );

    // Compute acceleration
    pointA_ = pointDivSigma_/pointRho_ - dampingCoeff_*pointU_ + g();

    // Check energies
    // To-do

    return true;
}


void explicitVertexCentredNonLinGeomTotalLagSolid::setTraction
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
#ifdef OPENFOAMESIORFOUNDATION
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
            "void Foam::explicitVertexCentredNonLinGeomTotalLagSolid::setTraction\n"
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


void explicitVertexCentredNonLinGeomTotalLagSolid::writeFields(const Time& runTime)
{
    Info<< nl << "Writing fields to " << runTime.timeName() << endl;

    // Write dualPf to a point field
    const volScalarField dualP(fvc::average(dualPf_));
    pointScalarField pointP
    (
        IOobject
        (
            "pointP",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    );
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map value from dual cell to point
        pointP[pointID] = dualP[dualCellI];
    }

    // Create the deformed normals point field
    Info<< "Writing dualGradDf etc" << endl;
    pointVol_.write();
    pointDivSigma_.write();
    dualGradDf_.write();
    dualFf_.write();
    dualFinvf_.write();
    dualFf_.write();
    const surfaceVectorField deformedSf(dualJf_*dualFinvf_.T() & dualMesh().Sf());
    const surfaceScalarField deformedMagSf(mag(deformedSf));
    const surfaceVectorField dualDeformedNormals(deformedSf/deformedMagSf);
    pointVectorField pointN
    (
        IOobject
        (
            "pointN",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimless, vector::zero)
    );
    const labelListList& pointToDualFaces = dualMeshMap().pointToDualFaces();
    forAll(pointN.boundaryField(), patchI)
    {
        if (pointN.boundaryField()[patchI].type() != "empty")
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            // Map from primary mesh point field to dual mesh face field using
            // the pointToDualFaces map
            for(int pI = 0; pI < meshPoints.size(); pI++)
            {
                const label pointID = meshPoints[pI];
                const labelList& curDualFaces = pointToDualFaces[pointID];

                forAll(curDualFaces, dfI)
                {
                    const label dualFaceID = curDualFaces[dfI];

                    if (!dualMesh().isInternalFace(dualFaceID))
                    {
                        // Check which patch this dual face belongs to
                        const label dualPatchID =
                            dualMesh().boundaryMesh().whichPatch(dualFaceID);

                        if (dualPatchID == patchI)
                        {
                            // Find local face index
                            const label localDualFaceID =
                                dualFaceID
                              - dualMesh().boundaryMesh()[dualPatchID].start();

                            // Dual face deformed unit normal
                            const vector n =
                                dualDeformedNormals.boundaryField()
                                [
                                    patchI
                                ][localDualFaceID];

                            pointN[pointID] = n;
                        }
                    }
                }
            }
        }
    }
    pointN.write();

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

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate
    mechanical().correct(sigma());

#ifdef OPENFOAMESI
    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());
#endif

//     // Calculate gradD at the primary points using least squares: this should
//     // be second-order accurate (... I think).
//     const pointTensorField pGradD(vfvc::pGrad(pointD(), mesh()));

//     // Calculate strain at the primary points based on pGradD
//     // Note: the symm operator is not defined for pointTensorFields so we will
//     // do it manually
//     // const pointSymmTensorField pEpsilon("pEpsilon", symm(pGradD));
//     pointSymmTensorField pEpsilon
//     (
//         IOobject
//         (
//             "pEpsilon",
//             runTime.timeName(),
//             runTime,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         pMesh(),
//         dimensionedSymmTensor("0", dimless, symmTensor::zero)
//     );

// #ifdef FOAMEXTEND
//     pEpsilon.internalField() = symm(pGradD.internalField());
// #else
//     pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
// #endif
//     pEpsilon.write();

//     // Equivalent strain at the points
//     pointScalarField pEpsilonEq
//     (
//         IOobject
//         (
//             "pEpsilonEq",
//             runTime.timeName(),
//             runTime,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         pMesh(),
//         dimensionedScalar("0", dimless, 0.0)
//     );

// #ifdef FOAMEXTEND
//     pEpsilonEq.internalField() =
//         sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
// #else
//     pEpsilonEq.primitiveFieldRef() =
//         sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
// #endif
//     pEpsilonEq.write();

//     Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << nl << endl;

    solidModel::writeFields(runTime);

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
